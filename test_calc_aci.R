# Test A/Ci curve fitting

source("calc_aci.R")

test <- data.frame("Photo" = c(0.5, 8, 14, 20, 25, 27, 28, 30, 31, 31.5, 32, 32.8, 33, 33.8, 34, 34.6, 35, 35.3, 35.5, 36), 
                   "Ci" = seq(50, 1000, 50), 
                   "Rd" = rep(1, 20),
                   "PAR" = rep(1000, 20),
                   "Tleaf" = rep(25, 20))

ggplot(test, aes(x = Ci, y = Photo)) + geom_point()


calc_vmax(data = test, ci_trans = 400)

