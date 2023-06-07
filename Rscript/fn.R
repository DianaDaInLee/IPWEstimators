impose_floor <- function(a, amin){
  new <- pmax(a, amin)
  total_slack <- sum(new) - 1
  individual_slack <- new - amin
  c <- total_slack / sum(individual_slack)
  new - c * individual_slack
}

draw_thompson_sampling <- function(batch, TYPE, K_ = K_s, S_ = S_s, t, s, W, YOBS, PMAT){
  if (TYPE == 'TS, BW'){
    PMAT <- PMAT[1:max(batch),]
    YOBS <- YOBS[1:max(batch)]
    W    <- W[1:max(batch)]
    N_    <- max(batch)
    
    # ipw       <- 1 / PMAT[cbind(1:max(batch), W[1:max(batch)])] # = 1 / pie
    # fit       <- lm_robust(YOBS[1:max(batch)] ~ as.factor(W[1:max(batch)]) -1, weights = ipw, se_type = 'HC2')
    # mu_hat    <- fit$coefficients
    # se_hat    <- fit$std.error
    
    ipw_score <- sapply(seq(K_), function(k) YOBS * (W == k) / PMAT[,k]) 
    mu_hat <- apply(ipw_score, 2, mean)
    se_hat <- sqrt(sapply(seq(K_), function(k) {(1 / (N_ * (N_-1))) * sum((ipw_score[,k] - mu_hat[k])^2)}))
    
    draws     <- replicate(S_, rnorm(K_, mean = mu_hat, sd = se_hat))
  }
  else{
    draws  <- replicate(S_, rbeta(K_, s + 1, t - s + 1)) 
  }
  return(draws)
}

sim_iter <- function(periods = periods_s, N = N_s, K = K_s, n_first = n_first_s, true_theta = true_theta_s, 
                     iter = iter_s, S = S_s, FLR = floor_s, type = NULL){
  
  sample_mean    <- sample_var <- ipw_mean <- haj_mean <- ipw_var <- haj_var <- aipw_mean <- aipw_var <- awaipw_mean <- awaipw_var <- posterior_final <- matrix(nrow = iter, ncol = K)
  posterior_full <- n_full <- NULL
  
  for (j in 1:iter){
    cat(paste(j, '...'))
    ymat <- mapply(rbinom, n = N, size = 1, prob = true_theta)       # hypothetical full potential outcomes matrix based on true mean (N X K)
    yobs <- w <- rep(NA, N)                                          # initialize vectors for observed outcomes (yobs) and treatment assignment (w) (N X 1)
    pmat <- mhat <- matrix(nrow = N, ncol = K)                       # initialize matrix to store treatment assignment probabilities (N row X K col)
    posterior <- matrix(nrow = periods + 1, ncol = K)
    nmat  <- matrix(nrow = periods, ncol = K)

    if ((N - n_first)%%(periods - 1) != 0){
      print("Error: Make sure to choose n_first such that the remaining sample size can be divisible evenly throughout the remaining periods.")
      stop()
    }
    n_oth <- (N - n_first) / (periods - 1)
        
    for (period in 1:periods){
      if (period == 1){
        # INITIAL PERIOD
        first_batch         <- seq(n_first)                                      # index for first batch
        pmat[first_batch, ] <- 1 / K                                             # randomly assign treatment in the first batch; observe outcomes
        posterior[period,]  <- 1 / K 
        mhat[first_batch, ] <- 0
        w[first_batch]      <- apply(pmat[first_batch, ], 1, function(x) sample(K, size = 1, prob = x)) # randomly assign with equal prob
        yobs[first_batch]   <- ymat[cbind(first_batch, w[first_batch])]          # get potential outcome for each unit in first batch under assigned K
        
        trials    <- sapply(seq(K), function(k) sum(w == k, na.rm = TRUE))       # number of assignments to each arm in the first batch (K X 1)
        successes <- sapply(seq(K), function(k) sum(yobs[w == k], na.rm = TRUE)) # number of successes under each arm in the first batch (K X 1)
        
        # Approximate Thompson Sampling probabilities
        draws    <- draw_thompson_sampling(first_batch, TYPE = type, t = trials, s = successes, W = w, YOBS = yobs, PMAT = pmat)
        argmax   <- apply(draws, 2, which.max)  # Check how many times w^th arm was maximal
        probs_ts <- table(cut(argmax, 0:K)) / S # Tally up the posterior probabilities
        probs_ts <- impose_floor(probs_ts, FLR)
        posterior[period+1,] <- probs_ts
        nmat[period,] <- trials
      }
      else{
        # SUBSEQUENT PERIODS
        next_batch           <- (sum(trials) + 1):((sum(trials) + 1) + n_oth - 1)
        mhat[next_batch,]    <- successes / trials
        if (type == 'Static'){
          pmat[next_batch, ] <- 1 / K                                                                     # static
          w[next_batch]      <- apply(pmat[next_batch, ], 1, function(x) sample(K, size = 1, prob = x))         
        }
        else{
          pmat[next_batch, ] <- rep(probs_ts, each = length(next_batch))                                  # thompson sampling
          w[next_batch]      <- sample(K, size = n_oth, prob = probs_ts, replace = TRUE)         
        }
        
        yobs[next_batch] <- ymat[cbind(next_batch, w[next_batch])]                                        # observe successes
        trials    <- sapply(seq(K), function(k) sum(w == k, na.rm = TRUE))                                # cumulative trials
        successes <- sapply(seq(K), function(k) sum(yobs[w == k], na.rm = TRUE))                          # cumulative successes
        
        draws     <- draw_thompson_sampling(batch = next_batch, TYPE = type, t = trials, s = successes, W = w, YOBS = yobs, PMAT = pmat)
        argmax    <- apply(draws, 2, which.max)
        probs_ts  <- impose_floor(table(cut(argmax, 0:K)) / S, FLR)
        posterior[period+1,] <- probs_ts
        nmat[period, ] <- trials
      }
    }
    
    # AFTER THE ENTIRE EXPERIMENT: WEIGHTING ESTIMATORS
    
    # 1) Naive Estimator
    sample_mean[j, ] <- successes / trials
    sample_var[j, ]  <- ((successes / trials) * (1 - successes / trials)) / trials
    
    # 2) IPW Estimator
    ipw           <- 1 / pmat[cbind(1:N, w)] # = 1 / pie (N X 1)
    ipw_score     <- sapply(seq(K), function(k) yobs * (w == k) / pmat[,k]) 
    ipw_mean[j, ] <- apply(ipw_score, 2, mean)
    # alternative: ipw_mean[j, ] <- sapply(seq(K), function(k) sum((yobs * ipw / N) * (w == k)))
    ipw_var[j, ]  <- sapply(seq(K), function(k) {(1 / (N * (N-1))) * sum((ipw_score[,k] - ipw_mean[j, k])^2)})
    
    # 3) Stabilized IPW Estimator (Hajek)
    # fit           <- lm_robust(yobs ~ as.factor(w) - 1, weights = ipw, se_type = 'HC2')
    # haj_mean[j, ] <- fit$coefficients
    # haj_var[j, ]  <- fit$std.error^2
    haj_score     <- sapply(seq(K), function(k) {(yobs * (w == k) / pmat[,k]) / sum(1 / N * (w == k) / pmat[,k]) })
    haj_mean[j, ] <- apply(haj_score, 2, mean)
    # alternative: haj_mean[j, ] <- sapply(seq(K), function(k) {sum((yobs * ipw / N) * (w == k) / (sum((ipw / N) * (w == k)))) })
    haj_var[j, ]  <- sapply(seq(K), function(k) {(1 / (N * (N-1))) * sum((haj_score[,k] - haj_mean[j, k])^2)})
    
    # 4) AIPW Estimator
    aipw_score     <- sapply(seq(K), function(k) yobs * (w == k) / pmat[,k] + (1 - (w == k) / pmat[,k]) * mhat[,k])
    aipw_mean[j,]  <- apply(aipw_score, 2, mean)
    aipw_var[j,]   <- sapply(seq(K), function(k) {(1 / (N * (N-1))) * sum((aipw_score[,k]- aipw_mean[j, k])^2)})
    
    # 5) Adaptively Weighted AIPW Estimator
    wts_lvdl        <- sqrt(pmat) 
    awaipw_mean[j,] <- sapply(seq(K), function(k) {sum(aipw_score[,k] * wts_lvdl[,k]) / sum(wts_lvdl[,k])})
    awaipw_var[j,]  <- sapply(seq(K), function(k) {sum(wts_lvdl[,k]^2 * (aipw_score[,k] - awaipw_mean[j,k])^2) / sum(wts_lvdl[,k])^2})
    
    posterior_final[j,] <- posterior[periods+1,]
    posterior_full <- rbind(posterior_full, as.data.frame(cbind(posterior, 0:periods)))
    n_full         <- rbind(n_full, as.data.frame(cbind(nmat, 1:periods)))
  }
  
  names(posterior_full) <- c(paste0('T', 1:K), 'period')
  names(n_full)         <- c(paste0('T', 1:K), 'period')
  
  results <- list(sample_mean      = sample_mean,
                  sample_var       = sample_var,
                  ipw_mean         = ipw_mean, 
                  ipw_var          = ipw_var, 
                  haj_mean         = haj_mean,
                  haj_var          = haj_var,
                  aipw_mean        = aipw_mean,
                  aipw_var         = aipw_var,
                  awaipw_mean      = awaipw_mean,
                  awaipw_var       = awaipw_var,
                  posterior_final  = posterior_final,
                  posterior_full   = posterior_full,
                  n_full           = n_full)
  return(results)
}

sim_iter_equalsize <- function(periods = periods_s, N = N_s, K = K_s, true_theta = true_theta_s, 
                     iter = iter_s, S = S_s, FLR = floor_s, type = NULL){
  
  sample_mean    <- sample_var <- ipw_mean <- haj_mean <- ipw_var <- haj_var <- aipw_mean <- aipw_var <- awaipw_mean <- awaipw_var <- posterior_final <- matrix(nrow = iter, ncol = K)
  posterior_full <- n_full <- comparison <- NULL
  
  for (j in 1:iter){
    cat(paste(j, '...'))
    ymat <- mapply(rbinom, n = N, size = 1, prob = true_theta)       # hypothetical full potential outcomes matrix based on true mean (N X K)
    yobs <- w <- rep(NA, N)                                          # initialize vectors for observed outcomes (yobs) and treatment assignment (w) (N X 1)
    pmat <- mhat <- matrix(nrow = N, ncol = K)                       # initialize matrix to store treatment assignment probabilities (N row X K col)
    posterior <- matrix(nrow = periods + 1, ncol = K)
    nmat  <- matrix(nrow = periods, ncol = K)
    
    n_oth <- floor(N / periods)
    n_first <- N - n_oth * (periods - 1)
  
    for (period in 1:periods){
      if (period == 1){
        # INITIAL PERIOD
        first_batch         <- seq(n_first)                                      # index for first batch
        pmat[first_batch, ] <- 1 / K                                             # randomly assign treatment in the first batch; observe outcomes
        posterior[period,]  <- 1 / K 
        mhat[first_batch, ] <- 0
        w[first_batch]      <- apply(pmat[first_batch, ], 1, function(x) sample(K, size = 1, prob = x)) # randomly assign with equal prob
        # for power analysis, ensure that every arm gets N_k=10
        if (type == 'Power, TS'){
          w[first_batch] <- as.numeric(randomizr::complete_ra(N = n_first, num_arms = K))
        }
        yobs[first_batch]   <- ymat[cbind(first_batch, w[first_batch])]          # get potential outcome for each unit in first batch under assigned K
      
        trials    <- sapply(seq(K), function(k) sum(w == k, na.rm = TRUE))       # number of assignments to each arm in the first batch (K X 1)
        successes <- sapply(seq(K), function(k) sum(yobs[w == k], na.rm = TRUE)) # number of successes under each arm in the first batch (K X 1)
    
        # Approximate Thompson Sampling probabilities
        draws    <- draw_thompson_sampling(batch = first_batch, K_ = K, S_ = S, TYPE = type, t = trials, s = successes, W = w, YOBS = yobs, PMAT = pmat)
        argmax   <- apply(draws, 2, which.max)  # Check how many times w^th arm was maximal
        probs_ts <- table(cut(argmax, 0:K)) / S # Tally up the posterior probabilities
        probs_ts <- impose_floor(probs_ts, FLR)
        posterior[period+1,] <- probs_ts
        nmat[period,] <- trials
      }
      else{
        # SUBSEQUENT PERIODS
        next_batch           <- (sum(trials) + 1):((sum(trials) + 1) + n_oth - 1)
        mhat[next_batch,]    <- rep(successes / trials, each = length(next_batch))
        if (type == 'Static'){
          pmat[next_batch, ] <- 1 / K                                                                     # static
          w[next_batch]      <- apply(pmat[next_batch, ], 1, function(x) sample(K, size = 1, prob = x))         
        }
        else{
          pmat[next_batch, ] <- rep(probs_ts, each = length(next_batch))                                  # thompson sampling
          w[next_batch]      <- sample(K, size = n_oth, prob = probs_ts, replace = TRUE)         
        }
        
        yobs[next_batch] <- ymat[cbind(next_batch, w[next_batch])]                                        # observe successes
        trials    <- sapply(seq(K), function(k) sum(w == k, na.rm = TRUE))                                # cumulative trials
        successes <- sapply(seq(K), function(k) sum(yobs[w == k], na.rm = TRUE))                          # cumulative successes
        
        draws     <- draw_thompson_sampling(batch = next_batch, K_ = K, S_ = S, TYPE = type, t = trials, s = successes, W = w, YOBS = yobs, PMAT = pmat)
        argmax    <- apply(draws, 2, which.max)
        probs_ts  <- impose_floor(table(cut(argmax, 0:K)) / S, FLR)
        posterior[period+1,] <- probs_ts
        nmat[period, ] <- trials
      }
    }
    
    # AFTER THE ENTIRE EXPERIMENT: WEIGHTING ESTIMATORS
    
    # 1) Naive Estimator
    sample_mean[j, ] <- successes / trials
    sample_var[j, ]  <- ((successes / trials) * (1 - successes / trials)) / trials
    
    # 2) IPW Estimator
    ipw           <- 1 / pmat[cbind(1:N, w)] # = 1 / pie (N X 1)
    ipw_score     <- sapply(seq(K), function(k) yobs * (w == k) / pmat[,k]) 
    ipw_mean[j, ] <- apply(ipw_score, 2, mean)
    # alternative: ipw_mean[j, ] <- sapply(seq(K), function(k) sum((yobs * ipw / N) * (w == k)))
    ipw_var[j, ]  <- sapply(seq(K), function(k) {(1 / (N * (N-1))) * sum((ipw_score[,k] - ipw_mean[j, k])^2)})
    
    # 3) Stabilized IPW Estimator (Hajek)
    # fit           <- lm_robust(yobs ~ as.factor(w) - 1, weights = ipw, se_type = 'HC2')
    # haj_mean[j, ] <- fit$coefficients
    # haj_var[j, ]  <- fit$std.error^2
    haj_score     <- sapply(seq(K), function(k) {(yobs * (w == k) / pmat[,k]) / sum(1 / N * (w == k) / pmat[,k]) })
    haj_mean[j, ] <- apply(haj_score, 2, mean)
    # alternative: haj_mean[j, ] <- sapply(seq(K), function(k) {sum((yobs * ipw / N) * (w == k) / (sum((ipw / N) * (w == k)))) })
    haj_var[j, ]  <- sapply(seq(K), function(k) {(1 / (N * (N-1))) * sum((haj_score[,k] - haj_mean[j, k])^2)})
    
    # 4) AIPW Estimator
    aipw_score     <- sapply(seq(K), function(k) yobs * (w == k) / pmat[,k] + (1 - (w == k) / pmat[,k]) * mhat[,k])
    aipw_mean[j,]  <- apply(aipw_score, 2, mean)
    aipw_var[j,]   <- sapply(seq(K), function(k) {(1 / (N * (N-1))) * sum((aipw_score[,k]- aipw_mean[j, k])^2)})
    
    # 5) Adaptively Weighted AIPW Estimator
    wts_lvdl        <- sqrt(pmat) 
    awaipw_mean[j,] <- sapply(seq(K), function(k) {sum(aipw_score[,k] * wts_lvdl[,k]) / sum(wts_lvdl[,k])})
    awaipw_var[j,]  <- sapply(seq(K), function(k) {sum(wts_lvdl[,k]^2 * (aipw_score[,k] - awaipw_mean[j,k])^2) / sum(wts_lvdl[,k])^2})
    
    # 6) Difference in Means
    for (ii in 2:K){
      dim  <- awaipw_mean[j,1] - awaipw_mean[j,ii]
      diff <- sapply(seq(K), function(k) {(aipw_score[,k] - awaipw_mean[j,k])})
      hsum <- sapply(seq(K), function(k) {sum(wts_lvdl[,k])})
      numerator <- hsum[ii] * wts_lvdl[1] * diff[,1] - hsum[1] * wts_lvdl[,ii] * diff[,ii]
      numerator <- sum(numerator^2)
      denominator <- hsum[1]^2 * hsum[ii]^2
      
      se <- sqrt(numerator / denominator)
      comparison <- bind_rows(comparison, data.frame(bestarm = 1, otherarm = ii, est = dim, se = se))
    }
    
    posterior_final[j,] <- posterior[periods+1,]
    posterior_full <- rbind(posterior_full, as.data.frame(cbind(posterior, 0:periods)))
    n_full         <- rbind(n_full, as.data.frame(cbind(nmat, 1:periods)))
  }
  
  names(posterior_full) <- c(paste0('T', 1:K), 'period')
  names(n_full)         <- c(paste0('T', 1:K), 'period')
  
  results <- list(sample_mean      = sample_mean,
                  sample_var       = sample_var,
                  ipw_mean         = ipw_mean, 
                  ipw_var          = ipw_var, 
                  haj_mean         = haj_mean,
                  haj_var          = haj_var,
                  aipw_mean        = aipw_mean,
                  aipw_var         = aipw_var,
                  awaipw_mean      = awaipw_mean,
                  awaipw_var       = awaipw_var,
                  posterior_final  = posterior_final,
                  posterior_full   = posterior_full,
                  n_full           = n_full,
                  comparison       = comparison)
  return(results)
}

plotit <- function(data, var){
  data$var <- data[[var]]
  xlabel <- 'Bias (True - Estimated)'
  if (var == 'ci_radius') xlabel <- 'Confidence Interval Radius (1.96 * S.E.)'
  p <- data %>%
    group_by(e, type, arm) %>%
    summarize( v = mean(var, na.rm = T),
               l = quantile(var, probs = 0.025, na.rm = T),
               h = quantile(var, probs = 0.975, na.rm = T)) %>%
    ggplot(aes(x = v, y = e, xmin = l, xmax = h)) +
    geom_point() + 
    geom_errorbarh(height = 0) + 
    facet_grid(type ~ arm, scales = 'free') +
    xlab(xlabel) + ylab('') +
    theme_bw() +
    theme(legend.position = 'none', 
          plot.caption = element_text(hjust = 0))  
  return(p)
}
