# IPWEstimators
This repository provides an R code that calculates and compares different IPW-estimators and Thompson sampling methods in a response-adaptive experiment. In contrast to conventional experiment, adaptive design dynamically updates treatment assignment probabilities based on the performance of the prior rounds, hence allocating larger share of the sample in more promising arms. 

The R code provides a function that simulate an adaptive-response trial and calculates the mean, variance, and confidence intervals of each treatment arm using the above 4 estimators. The function allows users to manipulate simulated experimental design (i.e., treatment assignment probability algorithm, number of treatment arms, number of periods, floor rates, etc.) to compare estimates under different scenarios.

The simulation allows for three different ways to calculate treatment assignment probabilities in each round:

1) Static: assigns equal treatment assignment probability across all arms for all rounds (i.e., equivalent of the conventional static experiment)
2) Thompson sampling: assigns treatment assignment probability based on standard Thompson sampling (assuming posterior probability following beta distribution)
3) Thompson sampling with balancing weights: (see https://mollyow.shinyapps.io/adaptive/#4_Instability_and_balancing)

Different estimators tested include: 

1) Simple mean
2) Inverse-Probability Weighted (IPW) Estimator (Offer-Westort et al 2021) (also see https://mollyow.shinyapps.io/adaptive/)
3) Hajek Estimator (Offer-Westort et al 2021) 
4) Augmented IPW Estimator (Hadad et al 2021)
5) Adaptively-weighted Augmented IPW Estimator (Hadad et al 2021)


# Reference
Hadad, V., Hirshberg, D. A., Zhan, R., Wager, S., & Athey, S. (2021). [Confidence intervals for policy evaluation in adaptive experiments](https://doi.org/10.1073/pnas.2014602118). Proceedings of the national academy of sciences, 118(15), e2014602118.

Offer‐Westort, M., Coppock, A., & Green, D. P. (2021). [Adaptive experimental design: Prospects and applications in political science](https://doi.org/10.1111/ajps.12597). American Journal of Political Science, 65(4), 826-844.
