## Project title
Integration of Immunoglobulin A and Faecal Egg Counts to Evaluate Gastrointestinal Nematode Burden and Control Strategies in Sheep
# Project description
To assess the various host and management factors on the burden of GIN infections, the study employs a Bayesian Negative Binomial model implemented through the JAGS (Just Another Gibbs Sampler) framework using runjags (Denwood, 2016) in R version 4.4.1(R Core Team, 2024) and modeled to account for individual variation, infection dynamics, and treatment effects. The model aimed to quantify the effects of the different factors on post-treatment FEC of gastrointestinal nematode (GIN) infections in Nigerian indigenous sheep, with FEC and serum IgA levels before and after treatment incorporated as covariates to improve model precision. The model included several categorical risk factors: location, breed, farm type, management system, sex, and age. These were incorporated as fixed effects in the model, with appropriate reference categories for each factor. The model also accounts for the variability in host susceptibility and treatment efficacy across these factors. The primary outcome variable, FEC after treatment, was modeled using a Gamma-Poisson (negative binomial) regression due to the overdispersion of the count data. The model was hierarchically structured to account for individual-level variation in infection status and treatment response.
## Features
- Bayesian hierarchical modeling with JAGS.
- Analysis of treatment effects on parasite burden.
- Exploration of host traits like sex, age, and breed.
- Support for overdispersed fecal egg count data.
## Requirements
- R
- Packages: `runjags`, `rjags`, `ggplot2`
