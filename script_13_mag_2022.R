library(vegan)

data(BCI)
data(BCI.env)

dist(BCI.env[, c(1:2)])
geo_dist <- dist(BCI.env[, c("UTM.EW", "UTM.NS")])

library(betapart)
?decostand
BCI_pa <- decostand(BCI, method = "pa")

?beta.pair
dist <- beta.pair(BCI_pa)

dist$beta.sor

plot(geo_dist, dist$beta.sor, 
     xlab = "Geographical distance (m.)",
     ylab = "Beta-Diversity: Sorensen",
     main = "Distance decay of BCI")

mod <- glm(dist$beta.sor ~ geo_dist, family = "binomial")
summary(mod)

?predict
xs <- seq(0, 1000, 1)

ys <- predict(mod, list(geo_dist = xs), type = "response")

lines(xs, ys, col = "red", lwd = 3) 

??model.decay

mod_2 <- decay.model(dist$beta.sor, geo_dist)
mod_3 <- decay.model(dist$beta.sor, geo_dist, model.type = "pow")
summary(mod_2$model)
summary(mod_3$model)
plot.decay(mod_2)

# 1. Calculate distances
# 2. Plot points
# 3. Fit a model
# 4. Make a prediction
# 5. Plot the prediction

data(mite)
data(mite.xy)

# 1
mite_pa <- decostand(mite, method = "pa")
eco_dist <- beta.pair(mite_pa)
geo_dist <- dist(mite.xy)

# 2
plot(geo_dist, eco_dist$beta.sor,
     xlab = "Geographical distance",
     ylab = "Beta Diversity")

# 3
mod <- glm(eco_dist$beta.sor ~ geo_dist, family = "binomial")
summary(mod)

# 4
xs <- seq(0, 10, 0.1)
ys <- predict(mod, list(geo_dist = xs), type = "response")

# 5
lines(xs, ys, col = "red", lwd = 3)

forest <- readRDS("data/probabilistic_sample.rds")
str(forest)
comm <- forest$plot
env <- forest$env

library(tidyverse)
comm %>% 
  select(c("id", "site", "veg")) %>% 
  distinct(site, veg) %>% 
  dim()
comm <- comm[, c(9:ncol(comm))]
colSums(comm)

comm <- comm[, colSums(comm) > 0]
env %>% dim()

alt_dist <- dist(env$altitude)
comm_pa <- decostand(comm, method = "pa")
eco_dist <- beta.pair(comm_pa, index.family = "sorensen")

plot(alt_dist, eco_dist$beta.sor,
     xlab = "Elevational shift",
     ylab = "Beta Diversity: Sorensen")

mant <- mantel(eco_dist$beta.sor, alt_dist)
summary(mant)


### Species pool and Dark diversity

library(DarkDiv)
BCI_pa
pool <- beals(BCI_pa)

BCI_pa %>% dim()
pool %>% dim()
pool[1:5, 1:5]
BCI_pa[1:5, 1:5]

max(pool)
min(pool)

pool_2 <- DarkDiv(BCI_pa, method = "RawBeals")
pool_2 %>% class()
pool_2$AllProbs[1:5, 1:5]
pool_2$Pool[1:5, 1:5]
pool_2$Dark[1:5, 1:5]

dd <- rowSums(pool_2$Dark, na.rm = TRUE)
dd
hist(dd,
     xlab = "Total Dark diversity")
sr <- specnumber(BCI_pa)

par(mfrow = c(1, 3))
hist(sr,
     xlab = "Species richness")
hist(dd,
     xlab = "Total Dark diversity")
hist(comp,
     xlab = "Completeness")

comp <- log(sr / dd)
comp
?specnumber
sr
dd

sp_frec <- specnumber(BCI_pa, MARGIN = 2)
pr <- colSums(pool_2$Pool)

range_filling <- round((sp_frec / pr) * 100, 1)
range_filling

hist(range_filling,
     xlab = "Range filling index",
     ylab = "Species frequencies")

?sar_power

### Species-Area relationship
library(sars)

data(galap)

galap %>% dim()

head(galap)
?galap
arr <- sar_power(galap)
arr
summary(arr)

arr <- sar_power(galap,
                 normaTest = "lillie",
                 homoTest = "cor.fitted")
summary(arr)
arr$normaTest

plot(arr,
     xlab = "Area (km²)")

str(arr)
pr <- arr$calculated
arr$residuals
plot(galap$a, arr$residuals)
abline(h = 0)

log_mod <- sar_loga(galap)
summary(log_mod)

multi_mod <- sar_multi(galap, obj = c("power", "loga", "koba"))

par(mfrow = c(2, 2))
plot(multi_mod$power)
plot(multi_mod$loga)
plot(multi_mod$koba)


?gdm

t <- c(4, 1, 13, 16, 15, 2, 6, 4, 5, 11, 3, 9, 8, 10, 12, 7)
length(t)
data(galap)
galap$t <- c(4, 1, 13, 16, 15, 2, 6, 4, 5, 11, 3, 9, 8, 10, 12, 7)
galap <- cbind(galap, t)

gdm(galap, model = "ATT2", mod_sel = TRUE)

#### INDICATOR SPECIES ANALISYS

forest <- readRDS("data/probabilistic_sample.rds")
forest %>% class()
comm <- forest$plot %>% 
  select(!(1:8))
comm <- comm[, colSums(comm) > 0]

comm_pa <- decostand(comm, method = "pa")

library(indicspecies)
?multipatt

indi_sp <- multipatt(comm_pa, forest$plot$veg, func = "r.g")

indi_sp %>% summary

forest_ita <- readRDS("data/forest_sample.rds")

forest_ita %>% class()
