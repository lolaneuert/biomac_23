install.packages("DarkDiv")

library(DarkDiv)
library(vegan)

dat <- read.csv("script/data_biomac_2.csv")
dim(dat)
dat[1:5, 1:10]

comm <- dat[, 6:ncol(dat)]

### DarkDiv package
?DarkDiv # dark diversity is the pool of species not present at the moment, but potentially occupying the area
# dark diversity concept is similar to beta diversity, both talk about and compare species that are not present,
# whilst in beta diversity the species are simply assumed to be absent in this locality, in dark diversity we talk about species that are not there but could be according to the parameters
dark_comm <- DarkDiv(comm, method = "RawBeals")
dark_comm_vegan <- vegan::beals(comm) # also works with the vegan package
str(dark_comm)
str(dark_comm_vegan)

dark_comm$AllProbs[1:5, 1:5]
dark_comm_vegan[1:5, 1:5]

max(dark_comm$AllProbs)
min(dark_comm$AllProbs)
mean(dark_comm$AllProbs)
median(dark_comm$AllProbs)

# very few species have an high probability of occurrence
# because few species occur frequently, and many species 
# occur rarely

quantile(dark_comm$AllProbs) 

### Dark Diversity s.s. 10.1111/geb.13203
rowSums(dark_comm$Dark) # NAs (i.e. missing values) should be accounted for when summing or doing other operations
dd <- rowSums(dark_comm$Dark, na.rm = T) # "na.rm = T" argument enables us to sum values by removing NAs
dd # our dark diversity vector (1 value per site/plot -> 19 values)
hist(dd,
     main = "Histogram of probabilistic dark diversity")
summary(dd)
str(dd)

boxplot(dd[1:10], dd[11:19]) # here we graphically compare dark diversity values between the first ten plots and last 9 plots
boxplot(dd[dat$Management == "F"], dd[dat$Management == "C"]) # here we do the same based on management
t.test(dd[dat$Management == "F"], dd[dat$Management == "C"]) # here we use a test to check for differences between the two groups
# it is not significant (above threshold of alpha), so the original hypothesis is accepted, the two box plots are similar/the means are similar enough


boxplot(dd[1:10], dd[11:19]) # here we graphically compare dark diversity values between the first ten plots and last 9 plots
boxplot(dd[dat$Region == "1"], dd[dat$Region == "2"]) # here we do the same based on management
t.test(dd[dat$Region == "1"], dd[dat$Region == "2"])

# for habitat, which is 4 levels, we need to use an anova instead of a t.test
dat$Habitat <- as.factor(dat$Habitat) # make habitat a factor
dd_matrix <- as.matrix(dd) # transform dd into a matrix format
boxplot(dd[dat$Habitat == "1"], dd[dat$Habitat == "2"], dd[dat$Habitat == "3"], dd[dat$Habitat == "4"])
anova <- (aov(dd_matrix~dat$Habitat))  # create an anova comparing the dd of different habitats
emmeans(anova, pairwise ~ Habitat) # see which of the habitats differ from each other 




### Community completeness -> ln(observed richness/dark diversity) (see 10.1007/s12224-013-9169-x)
# community completeness index, how much am I observing of the potential complete pool
# makes it easier to compare areas as they all show parts of a whole

sr <- specnumber(comm)
cc <- log(sr/dd) # community completeness index
cc
hist(cc)
summary(cc)

boxplot(cc[dat$Management == "F"], cc[dat$Management == "C"])
t.test(cc[dat$Management == "F"], cc[dat$Management == "C"])
# this time the p value is significant (below 0.05, therefore the alternative hypothesis is accepted, the means are not equal, the two plots are different from each other)

boxplot(cc[dat$Region == "1"], cc[dat$Region == "2"]) # here we do the same based on management
t.test(cc[dat$Region == "1"], cc[dat$Region == "2"])


cc_matrix <- as.matrix(cc)
boxplot(cc[dat$Habitat == "1"], cc[dat$Habitat == "2"], cc[dat$Habitat == "3"], cc[dat$Habitat == "4"])
anova_cc<- (aov(cc_matrix~dat$Habitat))  # create an anova comparing the dd of different habitats
TukeyHSD(anova_cc, pairwise ~ Habitat) # see which of the habitats differ from each other 



### Range filling -> realized/potential range size 10.1111/j.1461-0248.2004.00614.x
sf <- specnumber(comm, MARGIN = 2) # species frequencies -> realized range size
pr <- colSums(dark_comm$Pool) # potential range 
rf <- sf/pr
rf
head(sort(rf, decreasing = T), 20)
head(sort(rf), 20)
