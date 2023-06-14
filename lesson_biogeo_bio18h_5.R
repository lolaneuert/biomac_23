install.packages("sars")

library(sars) # species area relationships

data(galap) 
dim(galap)
head(galap)
?galap

arr_galap <- sar_power(data = galap, 
                       normaTest ="lillie", 
                       homoTest = "cor.fitted")
summary(arr_galap)
plot(arr_galap)
arr_galap$normaTest

gle_galap <- sar_loga(data = galap, normaTest ="lillie", homoTest = "cor.fitted")
summary(gle_galap)
plot(gle_galap)
gle_galap$normaTest
display_sars_models()
sars_models()

multi_galap <- sar_multi(data = galap, obj = c("power", "loga", "koba"))
summary(multi_galap$power)
summary(multi_galap$loga)
summary(multi_galap$koba)
plot(multi_galap)

sar_pred(arr_galap, area = 1000)

galap$t <- c(4, 1, 13, 16, 15, 2, 6, 4, 5, 11, 3, 9, 8, 10, 12, 7)
gdm(galap, model = "ATT2", mod_sel = TRUE)
gdm(galap, model = "all", mod_sel = TRUE)

?gdm
