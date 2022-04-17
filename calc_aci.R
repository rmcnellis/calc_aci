# Fit light-limited and Cyt b6f-limited portions of the A/Ci curve to calculate vcmax and vqmax

calc_gammastar = function(temp) {
  
  gammastar25 = 4.332
  Hgm = 37830 # J mol-1
  R = 8.314        # J K-1 mol-1
  
  temp_k = 273.15 + temp
  
  gStar = gammastar25 * exp((Hgm / R) * (1 / 298.15 - 1 / temp_k))
  
  gStar
  
}

calc_km_pa = function(temp) {
  
  R = 8.314        
  O2 = 2.09476e5      
  Kc25 = 41.03
  Ko25 = 28210
  Hkc = 79430  
  Hko = 36380 
  
  temp_k = 273.15 + temp
  
  Kc_pa = Kc25 * exp(Hkc * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  Ko_pa = Ko25* exp(Hko * ((temp_k - 298.15) / (298.15 * R * temp_k)))
  
  O2_pa = O2 * (1e-6) 
  
  Km_pa = Kc_pa * (1 + O2_pa/Ko_pa)
  
  Km_pa 
  
}

calc_vmax <- function(data, # must include Photo, Ci, Rd, PAR, Tleaf
                      ci_trans,
                      nl = 0.75,
                      nc = 1,
                      Abs = 0.85,
                      beta = 0.52,
                      Kf = 0.05e09,
                      Kd = 0.55e09,
                      Kp1 = 14.5e09,
                      start_vmax = 30,
                      control_maxiter = 500,
                      control_minFactor = 1/10000){
  
  library(minpack.lm)

  a2 <- Abs * beta
  a1 <- Abs - a2
  phi1P_max <- Kp1 / (Kp1 + Kd + Kf)
  a1_phi1P_max <- a1 * phi1P_max
  
  data$C <- data$Ci
  data$Rd <- mean(data$Rd)
  data$Q <- mean(data$PAR)
  gammas <- mean(calc_gammastar(data$Tleaf))
  K <- calc_km_pa(data$Tleaf)
  
  vqmax_function <- Photo ~ ((((vqmax * Q) / ((vqmax / (a1 * phi1P_max)) + Q)) /
                                (1 - (nl / nc) + (3 + 7 * (gammas / C)) / ((4 + 8 * gammas / C) * nc))) /
                               (4 + 8 * gammas / C)) * (1 - gammas / C) - Rd
  
  vcmax_function <- Photo ~ vcmax * ((C - gammas) / (C + K))
  
  vqmax_data <- subset(data, C > ci_trans)
  vcmax_data <- subset(data, C < ci_trans & C > 0)
  
  vqmax_nls <- nlsLM(vqmax_function, data = data,
                    start = list(vqmax = start_vmax),
                    control = nls.control(maxiter = control_maxiter, minFactor = control_minFactor))
  
  vcmax_nls <- nlsLM(vcmax_function, data = data,
                     start = list(vcmax = start_vmax),
                     control = nls.control(maxiter = control_maxiter, minFactor = control_minFactor))
  
  output <- data.frame("vmax" = c(coef(vcmax_nls), coef(vqmax_nls)), 
                       "CI_2.5" = c(confint(vcmax_nls)[1], confint(vqmax_nls)[1]),
                       "CI_97.5" = c(confint(vcmax_nls)[2], confint(vqmax_nls)[2]))

  return(output)
}

