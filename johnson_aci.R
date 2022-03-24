# Test Johnson model with A/Ci curves

library(dplyr)
library(tidyr)
library(purrr)
library(plantecophys)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(minpack.lm)
library (nlme)

source("../johnsonberry_photosynthesis_R/johnson_model.R")

# value for b parameter in temperature response curve for respiration 
# (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)
b_tresp_R <- 0.1781075 
# value for c parameter in temperature response curve for respiration 
# (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)
c_tresp_R <- -0.00179152 

aci_c3 <- list.files(path = "lce_c3_aci/",
                     pattern = "*.csv",
                     full.names = T) %>%
  map_df(~read.csv(.))

# aci_c3 <- read.csv("lce_c3_aci/Auburn_Ivom_1.csv")
# id_c3 <- aci_c3$id[1]
slce <- read.csv('LCE_data.csv')
# adjust Rd to similar temperature as the A/Ci data
slce$Rd <- slce$Rd / exp((b_tresp_R * (slce$Tleaf_R - slce$Tleaf_photo)) + 
                           (c_tresp_R * (slce$Tleaf_R^2 - slce$Tleaf_photo^2))) 
slce_Rd <- subset(slce, select = c("aci_id", "Rd"))
names(slce_Rd) <- c("id", "Rd")

aci_c3 <- left_join(aci_c3, slce_Rd, by = "id")
aci_c3 <- subset(aci_c3, !is.na(aci_c3$Photo))

calc_gammastar = function(temp) {
  
  rat = 100000 / 101325
  
  gammastar25 = 4.332 * rat  # Pa
  Hgm = 37830 # J mol-1
  R = 8.314        # J K-1 mol-1
  O2 = 2.09476e5 # ppm
  O2_0 = O2 * 1e-6 * 101325
  O2_z = O2 * 1e-6 * 100000
  
  temp_k = 273.15 + temp
  
  gStar_pa = gammastar25 * exp((Hgm / R) * (1 / 298.15 - 1 / temp_k))
  
  gStar_pa
  
}

nl <- 0.75
nc <- 1
gammas <- mean(calc_gammastar(aci_c3$Tleaf))
Abs <- 0.85
beta <- 0.52
a2 <- Abs * beta
a1 <- Abs - a2
Kf <- 0.05e09
Kd <- 0.55e09
Kp1 <- 14.5e09
phi1P_max <- Kp1 / (Kp1 + Kd + Kf)
aci_c3$C <- aci_c3$Ci / 1e6
Rd <- mean(aci_c3$Rd / 1e6)
Q <- mean(aci_c3$Pari / 1e6)
ci_trans <- 150 / 1e6 # estimated Ci transition point

A_function <- Photo ~ ((((vqmax * Q) / ((vqmax / (a1 * phi1P_max)) + Q)) / 
                          (1 - (nl / nc) + (3 + 7 * (gammas / C)) / ((4 + 8 * gammas / C) * nc))) /
                         (4 + 8 * gammas / C)) * (1 - gammas / C) - Rd

# fit the enzyme limited portion of the curve (Vcmax and Vpmax)
aci_c3_c <- subset(aci_c3, C < ci_trans & C > 0)
fit_function <- aci_c3_c %>% group_by(id) %>% 
  group_map(~nlsLM(formula = A_function, data = .x, start = list(vqmax = 0.01), 
                   control = nls.control(maxiter = 500, minFactor = 1/10000)), .keep = TRUE)
# fit_function <- nlsLM(A_function, data = aci_c3_c,
#                       start = list(vqmax = 0.01),
#                       control = nls.control(maxiter = 500, minFactor = 1 / 10000))
fit_function


vqmax_model <- model(PAR = aci_c3$Pari, Temp = aci_c3$Tleaf, CO2 = aci_c3$Ci,
                     CB6F = (1.2/1e6), RUB = (27.8/1e6))
max(vqmax_model$Vqmax) # maximum Cyt b6f activity, mol e-1 m-2 s-1

ggplot() + geom_point(data = aci_c3, aes(x = Ci, y = Photo)) +
  geom_line(data = vqmax_model, aes(x = C * 1e6, y = An_a * 1e6), linetype = "dashed")


