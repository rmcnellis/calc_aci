# Test Johnson model with A/Ci curves

library(dplyr)
library(tidyr)
library(purrr)
library(plantecophys)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(minpack.lm)
library(nlme)

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
  
  gammastar25 = 4.332
  Hgm = 37830 # J mol-1
  R = 8.314        # J K-1 mol-1

  temp_k = 273.15 + temp
  
  gStar = gammastar25 * exp((Hgm / R) * (1 / 298.15 - 1 / temp_k))
  
  gStar
  
}

nl <- 0.75
nc <- 1
Abs <- 0.85
beta <- 0.52
a2 <- Abs * beta
a1 <- Abs - a2
Kf <- 0.05e09
Kd <- 0.55e09
Kp1 <- 14.5e09
phi1P_max <- Kp1 / (Kp1 + Kd + Kf)
a1_phi1P_max <- a1 * phi1P_max

# A_function <- Photo ~ ((((vqmax * Q) / ((vqmax / (a1 * phi1P_max)) + Q)) / 
#                           (1 - (nl / nc) + (3 + 7 * (gammas / C)) / ((4 + 8 * gammas / C) * nc))) /
#                          (4 + 8 * gammas / C)) * (1 - gammas / C) - Rd

A_function <- Photo ~ (((vqmax * Q) / ((vqmax / a1_phi1P_max) + Q)) / 4) - Rd

vqmax <- 360
Q <- 1000 
Rd <- 1
Photo <- (((vqmax * Q) / ((vqmax / a1_phi1P_max) + Q)) / 4) - Rd
Photo

# fit the enzyme limited portion of the curve (Vcmax and Vpmax)
id_list <- unique(aci_c3$id)
i = 2
for (i in 1:length(id_list)){
  temp_data <- subset(aci_c3, id == id_list[i] & Ci > 375)
  # temp_data$C <- temp_data$Ci / 1e6
  temp_data$Rd <- mean(temp_data$Rd)
  temp_data$Q <- mean(temp_data$Pari)
  gammas <- mean(calc_gammastar(temp_data$Tleaf))
  temp_nls <- nlsLM(A_function, data = temp_data,
                    start = list(vqmax = 0.01),
                    control = nls.control(maxiter = 500, minFactor = 1 / 10000))
  summary(temp_nls)
  
}
ggplot(subset(aci_c3, aci_c3$id == id_list[2]), aes(x = Ci, y = Photo)) + geom_point()
# fit_function <- aci_c3_c %>% group_by(id) %>% 
#   group_map(~nlsLM(formula = A_function, data = .x, start = list(vqmax = 0.01), 
#                    control = nls.control(maxiter = 500, minFactor = 1/10000)), .keep = TRUE)
# fit_function <- nlsLM(A_function, data = aci_c3_c,
#                       start = list(vqmax = 0.01),
#                       control = nls.control(maxiter = 500, minFactor = 1 / 10000))
fit_function


vqmax_model <- model(PAR = aci_c3$Pari, Temp = aci_c3$Tleaf, CO2 = aci_c3$Ci,
                     CB6F = (1.2/1e6), RUB = (27.8/1e6))
max(vqmax_model$Vqmax) # maximum Cyt b6f activity, mol e-1 m-2 s-1

ggplot() + geom_point(data = aci_c3, aes(x = Ci, y = Photo)) +
  geom_line(data = vqmax_model, aes(x = C * 1e6, y = An_a * 1e6), linetype = "dashed")


