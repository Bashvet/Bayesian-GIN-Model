library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

gindata <- read.csv("GinJagsNig.csv",stringsAsFactors = T)

# Define the JAGS model
model_string <- "
model{
  for (i in 1:N) {
    # Infection status at two time points
    for (t in 1:2) {
      status[i, t] ~ dbern(p[t])
    }
    
    # Fecal Egg Count Before Treatment
    EggCountBefore[i] ~ dnegbin(k/(mean[i,1]*status[i,1]+k),k)
    
    # Use factors: Sex and Age, excluding reference category
    log(mean[i, 1]) <- beta0 +
                  beta_location[Location[i]] +
                  beta_breed[Breed[i]] +
                  beta_farm_type[Farm_type[i]] +
                  beta_mgt_system[Mgt_system[i]] +
                  beta_sex[Sex[i]] +  # Male is reference (Sex == 2); estimate for Female (Sex == 1)
                  beta_age[Age[i]]    # Adult is reference (Age == 1); estimate for Young (Young == 2)
    
    # Fecal Egg Count After Treatment
    EggCountAfter[i] ~ dnegbin(k2/(mean[i,2]*status[i,2]+k2),k2)
    
    # Mean after treatment
    mean[i,2] <- mean[i,1] * treatment_effect
    
    # Serum IgA Levels
    PreIgA[i] ~ dgamma(shape1[status[i,1]+1], rate1[status[i,1]+1])
    PostIgA[i] ~ dgamma(shape2[status[i,1]+1], rate2[status[i,1]+1])
  }
  
  # Priors
  for (t in 1:2) {
    p[t] ~ dbeta(10, 2)
  }
  treatment_effect ~ dbeta(1, 1)

  k ~ dgamma(0.001, 0.001)
  k2 ~ dgamma(0.001, 0.001)
  
  ### shape for uninfected
  #shape1[1] ~ dunif(0,100)
  #shape2[1] <- shape1[1]
  
  ## shape for infected
  #shape1[2] ~ dgamma(0.001, 0.001)
  #shape2[2] ~ dgamma(0.001, 0.001)
  shape1[1] <- shape1[2]
  shape1[2] ~ dgamma(0.001, 0.001)
  shape2[1] <- shape2[2]
  shape2[2] ~ dgamma(0.001, 0.001)
  
  for (j in 1:3) {
    unsortmeans[j] ~ dunif(0.001,100)
  }
  sortmeans <- sort(unsortmeans)

  # Incorporate IgA effect
  rate1[1] <- shape1[1] / sortmeans[1]
  rate1[2] <- shape1[2] / sortmeans[3]
  
  rate2[1] <- shape2[1] / sortmeans[1]
  rate2[2] <- shape2[2] / sortmeans[2]
  
  
  # Priors for regression coefficients
  beta0 ~ dnorm(0,0.001)
  
  for (j in 1:n_location) {
    beta_location[j] ~ dnorm(0, locprec)
  }
  locprec ~ dgamma(0.001,0.001)
  
  for (j in 1:n_breed) {
    beta_breed[j] ~ dnorm(0, breedprec)
  }
  breedprec ~ dgamma(0.001,0.001)
  
  beta_farm_type[1] <- 0
  for (j in 2:n_farm_type) {
    beta_farm_type[j] ~ dnorm(0, 0.001)
  }
  
  beta_mgt_system[1] <- 0
  for (j in 2:n_mgt_system) {
    beta_mgt_system[j] ~ dnorm(0, 0.001)
  }
  
  # Priors for Sex (only 1 level for comparison, i.e., Female)
  beta_sex[1] ~ dnorm(0, 0.001)
  beta_sex[2] <- 0
  
  # Priors for Age (only 1 level for comparison, i.e., Young)
  beta_age[1] ~ dnorm(0, 0.001)
  beta_age[2] <- 0
  
  #inits# beta0, beta_location, beta_breed, treatment_effect, p, status
  #monitor# beta0, beta_location, beta_breed, beta_sex, beta_age, beta_farm_type, beta_mgt_system, shape1, shape2, rate1,rate2, treatment_effect, k, k2, p
}"

# Define the number of unique levels for each factor
n_location <- length(unique(gindata$Location))
n_breed <- length(unique(gindata$Breed))
n_farm_type <- length(unique(gindata$Farm_type))
n_mgt_system <- length(unique(gindata$Mgt_system))
n_sex <- length(unique(gindata$Sex))
n_age <- length(unique(gindata$Age))

# Create the data list, now using the properly defined 'n_' variables
data_list <- list(
  N = nrow(gindata),
  EggCountBefore = gindata$EggCountBefore,
  EggCountAfter = gindata$EggCountAfter,
  PreIgA = gindata$PreIgA,
  PostIgA = gindata$PostIgA,
  Location = gindata$Location,
  Breed = gindata$Breed,
  Sex = gindata$Sex,
  Age = gindata$Age,
  Farm_type = gindata$Farm_type,
  Mgt_system = gindata$Mgt_system,
  n_location = n_location,
  n_breed = n_breed,
  n_farm_type = n_farm_type,
  n_mgt_system = n_mgt_system
)

# Define the initial values function
inits2 <- function() {
  list(
    beta0 = 1,
    beta_location = rep(0, n_location),
    beta_breed = rep(0, n_breed),
    #beta_sex = rep(0, 2),  # No need for the reference level, only one beta for comparison
    #beta_age = rep(0, 2),  # No need for the reference level, only one beta for comparison
    beta_farm_type = rep(0, n_farm_type),
    beta_mgt_system = rep(0, n_mgt_system),
    shape1 = 1,
    shape2 = 1,
    treatment_effect = 0.5,
    unsortmean = c(0,0.5,1),
    p = c(0.9, 0.9),  # Initialize P with 0.9
    status = matrix(1, nrow = data_list$N, ncol = 2)  # Initialize status to 1
  )
}

# Parameters to monitor
#params <- c("beta0", "beta_location", "beta_breed", 
#            "beta_sex", "beta_age", "beta_farm_type", 
#            "beta_mgt_system", "shape1", "shape2", "rate1","rate2",
#            "treatment_effect", "k", "k2", "p", "status")

# Run the model
results <- run.jags(
  model = model_string,
  data = data_list,
  inits = inits2,
  n.chains = 2,
  adapt = 1000,
  burnin = 5000,
  sample = 10000,
  thin = 10
)


# Display the results
print(results)
plot(results)
summary(results)

results$psrf

#save.image("Results.RData")

#####################################

setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
load("Results.RData")

