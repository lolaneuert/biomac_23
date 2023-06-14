install.packages(c("vegan", "Hmisc", "betapart"))

library(vegan)
library(Hmisc)
library(betapart)

# upload a dataset from a *.csv
dat <- read.csv("script/data_biomac_2.csv")
dim(dat)
names(dat)
head(names(dat), n = 10)
dat[1:5, 1:10]

# Hierarchical sampling design
# unique() to look at the different classes for each hierarchical category
unique(dat$Habitat)
unique(dat$Region)
unique(dat$Dataset)
unique(dat$Management)

comm <- dat[, -(1:5)] # to divide the dataset that is composed of 2 parts: species and environmental data, so we remove the first 5 columns to only have species
comm_pa <- decostand(comm, "pa") # now we can create presence/absence values

### Additive partitioning ###
?adipart
adipart(comm_pa) 

# Alpha: local diversity, mean species richness
# Gamma: regional diversity, total diversity
# Beta: missing species from the species pool, dark diversity

adi_mod <- adipart(comm_pa, dat[, c("Plot", "Habitat", "Region")])

adi_mod
# plot richness higher than expected
# missing species from the region: lower than expected

str(adi_mod)
adi_mod$statistic
obs <- adi_mod$statistic[c("alpha.1", "beta.1")]
expc <- adi_mod$oecosimu$means[c(1, 4)]
adi_mat <- matrix(data = c(obs, expc), ncol = 2)

barplot(adi_mat,
        space = 1,
        xlim = c(1, 5),
        ylab = "Number of species",
        names.arg = c("Observed", "Expected"),
        legend.text = c("Alpha plot", "Beta plot"),
        args.legend = list(x = "topright"))


# ...we go further across the hierarchical scheme
adipart(comm_pa, dat[, c(1, 5, 4)]) # plots (alpha1), management (alpha2) and dataset (gamma)
adipart(comm_pa, dat[, 1:4]) # plots (alpha1), habitat (alpha2), region (alpha3) and gamma (gamma)

### Multiplicative partitioning ###
?multipart
multipart(comm_pa, scales = 0)
multipart(comm_pa, dat[, c("Plot", "Habitat", "Region")], scales = 0)
multipart(comm_pa, dat[, c("Plot", "Habitat", "Region", "Dataset")], scales = 0)

mlt_mod <- multipart(comm_pa, dat[, 1:4])
str(mlt_mod)
obs <- mlt_mod$statistic[c(5, 6, 7)]
expc <- mlt_mod$oecosimu$means[c(5, 6, 7)]
multi_mat <- matrix(data = c(obs, expc), ncol = 2)

barplot(adi_mat,
        space = 1,
        xlim = c(1, 5),
        ylab = "Number of species",
        names.arg = c("Observed", "Expected"),
        legend.text = c("Alpha plot", "Beta plot"),
        args.legend = list(x = "topright"))

### Beta partitioning (Baselga framework) ###
### Load the dataset "Barrio Colorado Islands"

data(BCI)

### And check the structure

str(BCI)

### Exercise ###
### Calculate and plot species richness at plot-level for the entire dataset and
### convert matrix to presence/absence

sr <- specnumber(BCI)
hist(sr,
     breaks = 20,
     xlab = "Species Richness",
     ylab = "Number of plots",
     main = "")

BCI <- decostand(BCI, "pa") # transform the abundance matrix into a presence/absence matrix

### Calculate distance and similarity metrics among pairs of plots and multiple plots

vegdist(BCI[1:2, ], method = "jaccard") # Similarity index 
1 - vegdist(BCI[1:2, ], method = "jaccard") # Distance index
vegdist(BCI[2:3, ], method = "jaccard")
vegdist(BCI[1:3, ], method = "jaccard")
vegdist(BCI, method = "jaccard")
hist(vegdist(BCI, method = "jaccard"),
     xlab = "Jaccard's dissimilarity metric",
     ylab = "Number of pairs",
     main = "")

hist(vegdist(BCI, method = "bray"),
     xlab = "Sorensen's dissimilarity metric",
     ylab = "Number of pairs",
     main = "")

par(mfrow = c(1, 2))
hist(vegdist(BCI, method = "jaccard"),
     xlab = "Jaccard's dissimilarity metric",
     ylab = "Number of pairs",
     main = "")
hist(vegdist(BCI, method = "bray"),
     xlab = "Sorensen's dissimilarity metric",
     ylab = "Number of pairs",
     main = "")
dev.off()

### Partitioning beta diversity

?beta.multi # multiple site dissimilarity
beta.multi(BCI)
beta.multi(BCI[1:2, ]) # to compare certain rows
beta.multi(BCI[1:3, ])

?beta.pair # incidence-based pair-wise dissimilarities
beta.pair(BCI[1:3, ])
beta.pair(BCI)
boxplot(beta.pair(BCI))
boxplot(beta.pair(BCI),
        names = c("Turnover", "Nestedness", "Total beta"))

par(mfrow = c(1, 2))
boxplot(beta.pair(BCI),
        names = c("Turnover", "Nestedness", "Total beta"),
        main = "Sorensen")
boxplot(beta.pair(BCI, "jaccard"),
        names = c("Turnover", "Nestedness", "Total beta"),
        main = "Jaccard")
par(mfrow = c(1, 1))

### Distance decay

data(BCI.env)

spat_dist <- dist(BCI.env[, 1:2])
dissim_BCI <- beta.pair(BCI_pa)$beta.sor

plot(spat_dist, dissim_BCI$beta.sim, ylim = c(0.1, 0.6), xlim = c(0, max(spat_dist)))
pr <- decay.model(dissim_BCI$beta.sor, spat_dist, 
            model.type = "pow")
summary(pr$model)
summary(pr$data.y)
plot(pr)
bci_decay <- glm(dissim_BCI$beta.sor ~ spat_dist, family = "binomial")
plot(bci_decay)

summary(bci_decay)
bci_decay
xs <- seq(0, 1000, 10)
ys <- predict(bci_decay, list(spat_dist = xs), type = "response")
lines(xs, ys, col = "red", lwd = 3)



### EXERCISE
data(mite)
data(mite.xy)

mite_dist <- dist(mite.xy)

mite <- decostand(mite, "pa")
mite_dissim <- beta.pair(mite)$beta.sor

plot(mite_dist, mite_dissim, xlim = c(0, max(mite_dist)), las = 1,
     xlab = "Spatial distance (m.)", ylab = "Sorensen's dissimilarity")

mite_decay <- glm(mite_dissim ~ mite_dist, family = "binomial")
mite_decay
summary(mite_decay)

xs <- seq(0, 10, 0.1)
ys <- predict(mite_decay, list(mite_dist = xs), type = "response")

plot(mite_dist, mite_dissim, xlim = c(0, max(mite_dist)), las = 1,
     xlab = "Spatial distance (m.)", ylab = "Sorensen's dissimilarity")
lines(xs, ys, col = "red", lwd = 3)


forest <- readRDS("data/probabilistic_sample.rds")

comm_pa <- vegan::decostand(forest$plot[,c(9:ncol(forest$plot))], method = "pa")
eco_dist <- beta.pair(comm_pa)$turnover
hist(as.matrix(eco_dist))

plot(dist_for, eco_dist,  xlab = "Spatial distance (m.)", ylab = "Sorensen's dissimilarity")

forest$plot[, c(1:8)]

library(tidyverse)
sam <- forest$plot %>% dim()
  group_by(site) %>% 
  sample_n(1)
comm <- forest$plot[, 9:ncol(forest$plot)]
colSums(comm)

dist_for <- dist(forest$plot[, c(2:3)])
comm_c <- comm[, colSums(comm) > 0]
comm_pa <- vegan::decostand(comm_c, method = "pa")
eco_dist <- beta.pair(comm_pa)$beta.sor

plot(dist_for, eco_dist, 
     xlab = "Spatial distance (m.)", 
     ylab = "Sorensen's dissimilarity")

mite_decay <- glm(eco_dist ~ dist_for, family = "binomial")
mite_decay
summary(mite_decay)

xs <- seq(0, 10, 0.1)
ys <- predict(mite_decay, 
              list(dist_for = xs), type = "response")

plot(dist_for, eco_dist, 
     xlab = "Spatial distance (m.)", ylab = "Sorensen's dissimilarity")
lines(xs, ys, col = "red", lwd = 3)

sam
comm_pa <- forest$plot[, 9:ncol(forest$plot)] %>% 
  decostand(., method = "pa")
comm_pa <- comm_pa[, colSums(comm_pa) > 0]

adi_mod <- adipart(comm_pa, 
                   forest$plot[, c("id", "site", "veg")])
                   %>% 
                     ungroup())

ss <- read.csv("data/data_biomac_2.csv")
str(ss)
str(sam)
env <- forest$env
