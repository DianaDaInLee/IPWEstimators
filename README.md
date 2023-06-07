# IPWEstimators
This repository provides an R code that calculates and compares different IPW-estimators in response-adaptive experiment. Different estimators tested include: 

1) Simple mean
2) Hajek Estimator (Offer-Westort et al 2021) (also see https://mollyow.shinyapps.io/adaptive/)
3) Augmented IPW Estimator (Hadad et al 2021)
4) Adaptively-weighted Augmented IPW Estimator (Hadad et al 2021)

The R code runs a simulation of an adaptive-response trial and calculates the mean, variance, and confidence intervals of each treatment arm using the above 4 estimators. The function allows users to manipulate simulated experimental design (i.e., number of treatment arms, number of periods, floor rates, etc.) to compare estimates under different scenarios.
