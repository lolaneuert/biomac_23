### COMMUNITY MATRIX: vegan::dune CO-OCCURENCE DATA
###########################
library(vegan)
?vegan
# library(tidyverse)

data(dune)
data(dune.env)

# dimensions (rows, columns, cells)
# in R: data[row(s), column(s)]
dim(dune)
View(dune)

# co-occurrence species data
# each record/observation/vegetation plot (rows) corresponds
# to cover values of all the species, from 0 to max cover (9 in our case)
dune
dune[1:10, 1:10] 

# from co-occurence data to occurrence data using tidyverse
# each record/observation (rows) of occurrence data corresponds
# to cover values of a single species in a locality

library(tidyverse)
# ctrl + shift + M the shortcut for the pipe operator of tidyverse: %>% 

occurrence_dune <- dune %>% 
  rownames_to_column("id") %>% # create a new column using names of rows
  pivot_longer(!id, names_to = "sp", values_to = "cover") %>% # pivot table with species (present and absent) as single record for each plot
  filter(cover > 0) # exclude all absent species

occurrence_summary <- occurrence_dune %>% 
  group_by(sp) %>% 
  summarise(mean_cover = mean(cover),
            n = n())

### MULTIVARIATE ANALYSES: CLUSTERING AND ORDINATION
##########################

vegan::specnumber(dune) # built-in function of vegan library (called by :: to specify from which package) to calculate species richness

rowSums(dune) 
max(dune)

cover_level <- dune %>% 
  rownames_to_column("id") %>% 
  pivot_longer(!id, names_to = "sp", values_to = "cover") %>% 
  filter(cover > 0) %>% 
  pull(cover) %>% 
  unique() %>% 
  sort()

cover_level

# presence-absence matrix
dune_pa <- decostand(dune, method = "pa") # decostand with method pa transforms values into a presence/absence scale (0/1)
dune_pa[1:5, 1:7]

hist(colSums(dune_pa)) # ...many rare species

dune_pa %>% 
        colSums() %>% 
        hist()

dune_pa %>% 
        colSums() %>% 
        sort(decreasing = TRUE) 
        

# data transformation of species cover values
hist(as.matrix((dune))) # double-zero issue
hist(as.matrix((dune[dune > 0]))) # right skewed cover data

dune_log <- log1p(dune) # log transform
dune_sqrt <- sqrt(dune) # transform using square root

par(mfrow = c(2, 3))
hist(as.matrix((dune[dune > 0])),
     main = "Ordinal cover data",
     xlab = "Species cover values",
     ylab = "Frequency") # right skewed cover data
hist(as.matrix((dune_sqrt[dune_sqrt > 0])), 
     main = "Square root-transformed cover data",
     xlab = "Species cover values",
     ylab = "Frequency")
hist(as.matrix((dune_log[dune_log > 0])), 
     main = "Log-transformed cover data",
     xlab = "Species cover values",
     ylab = "Frequency")
car::qqPlot(dune[dune > 0])
car::qqPlot(dune_sqrt[dune_sqrt > 0])
car::qqPlot(dune_log[dune_log > 0])
dev.off()


### Ordination
# Principal component analysis --> reduces multidimensional data to principal components
### See detail on vegan::rda() at: 
# https://cran.r-project.org/web/packages/vegan/vegan.pdf

pca <- rda(dune_log)
pca

# ordination biplots
ordiplot(pca, display = "si")
ordiplot(pca, display = "sp", type = "text")


str(pca) # structure of pca object

pca$tot.chi
(pca$CA$eig[1] / pca$tot.chi) * 100 # explained variance by first axis

((pca$CA$eig[1] + pca$CA$eig[2]) / pca$tot.chi) * 100 # explained variance by first and second axes

summary(pca)


### Clustering
# plot partition according their species composition, using K-means clustering
set.seed(1221) # fix the random factor
clustering <- kmeans(dune_log, 3)
clustering
?kmeans

groups <- clustering$cluster # obtain group assignment

groups %>% 
        table()

# plot groups in the  multivariate space of the pca
ordiplot(pca, display = "si", type = "n")
for(i in 1:3) ordihull(pca, groups = groups, show.groups = i, col = i)
for(i in 1:3) ordispider(pca, groups = groups, show.groups = i, col = i)

# apply pca on environmental variables of our plot observations
head(dune.env)
str(dune.env)

par(mfrow = c(1, 3))

# Plot soil thickness
plot(pca, display = "sites", type = "n",
     main = "Soil thickness") 
with(dune.env, points(pca, disp = "si", cex = as.numeric(A1)))
with(dune.env, text(pca, disp = "si", labels = as.numeric(A1)))
ordisurf (pca, dune.env[, 'A1'], add = T, col = 'red', )

# Plot Moisture
plot(pca, display = "sites", type = "n", 
     main = "Soil Moisture")
with(dune.env, points(pca, disp = "si", pch = as.numeric(Moisture)))
with(dune.env, legend("topleft", levels(Moisture), pch = 1:4,
                      title = "Moisture"))
with(dune.env, ordispider(pca, Moisture, label = TRUE))

ordiplot(pca, display = "si", type = "n",
         main = "Species composition")
for(i in 1:3) ordihull(pca, groups = groups, show.groups = i, col = i)
for(i in 1:3) ordispider(pca, groups = groups, show.groups = i, col = i, label = T)
dev.off()
# Boxplot of species richness
boxplot_data <- dune %>% 
        decostand(., method = "pa") %>% 
        rowSums() %>% 
        as.data.frame() %>% 
        cbind(groups) %>% 
        rename(rich = ".") %>%
        mutate(groups = as.factor(groups)) 

boxplot(boxplot_data$rich ~ boxplot_data$groups)

vegan::specnumber(dune)

### INDICATOR SPECIES ANALYSIS
### Statistical Test for species co-occurrence patterns
###########################

# Before running the test, we take a look to species frequencies among groups
groups
groups == 1 # apply a condition

dune_pa <- decostand(dune, method = "pa")
group_1 <- dune_pa[groups == 1, ] # apply a condition to a dataset, on rows
# colSums(group_1) > 0 # new condition

group_1 <- group_1[, colSums(group_1) > 0] # apply a new condition to the dataset, on columns
# t(group_1) # transpose 
group_1 <- as.data.frame(t(group_1)) # transform the transpose of a matrix, to dataframe
group_1 <- sort(rowSums(group_1), decreasing = TRUE) # I sort species to be identify faster most frequent species
group_1 <- data.frame(row.names = c(1:length(group_1)), sp = names(group_1), fr = group_1) # I define a new dataframe with meaningful field and names

# Apply the same procedure to the other two groups (this could be done using a for() loop)
comm_pa <- decostand(dune, method = "pa")
group_2 <- comm_pa[groups == 2, ]
group_2 <- group_2[, colSums(group_2) > 0]
group_2 <- as.data.frame(t(group_2))
group_2 <- sort(rowSums(group_2), decreasing = TRUE) 
group_2 <- data.frame(row.names = c(1:length(group_2)), sp = names(group_2), fr = group_2)

# gr3...
comm_pa <- decostand(dune, method = "pa")
group_3 <- comm_pa[groups == 3, ]
group_3 <- group_3[, colSums(group_3) > 0]
group_3 <- as.data.frame(t(group_3))
group_3 <- sort(rowSums(group_3), decreasing = TRUE)
group_3 <- data.frame(row.names = c(1:length(group_3)), sp = names(group_3), fr = group_3)

### Merge the three dataframe using the "sp" field as link between dataframes
tot_species <- merge(x = group_1, y = group_2, by = "sp", all = TRUE)
tot_species <- merge(x = tot_species, y = group_3, by = "sp", all = TRUE)
tot_species[is.na(tot_species)] <- 0
names(tot_species) <- c("sp", "gr_1", "gr_2", "gr_3")


# install.packages("indicspecies")
library(indicspecies)

indi <- multipatt(dune_pa, groups) # see default parameters to be sure on the test you are using
# See tutorial at
# https://cran.r-project.org/web/packages/indicspecies/vignettes/indicspeciesTutorial.pdf

summary(indi)

### EXERCISE
# clustering (kmeans) of plots according the site chemical compositionn
# ordination plot (first and second pca axes) with site and "species"
# indicator species for each group

# Data set
data(varespec)
data(varechem)

# Main functions
kmeans()
rda()
multipatt()

# do it all for varespec
# presence-absence matrix
vare_pa <- decostand(varespec, method = "pa") # decostand with method pa transforms values into a presence/absence scale (0/1)

hist(as.matrix((varespec))) # double-zero issue
hist(as.matrix((varespec[varespec > 0]))) # right skewed cover data

vare_log <- log1p(varespec) # log transform
vare_sqrt <- sqrt(varespec) # transform using square root

pca_vare <- rda(vare_log)
pca_vare

# ordination biplots
ordiplot(pca_vare, display = "si", type= "text")
# ordiplot(pca_vare, display = "sp", type = "text")


str(pca_vare) # structure of pca object

pca$tot.chi
(pca$CA$eig[1] / pca$tot.chi) * 100 # explained variance by first axis

((pca$CA$eig[1] + pca$CA$eig[2]) / pca$tot.chi) * 100 # explained variance by first and second axes

summary(pca)

### Clustering
# plot partition according their species composition, using K-means clustering
set.seed(1221) # fix the random factor
clustering_vare <- kmeans(vare_log, 4)
clustering_vare
?kmeans

groups_vare <- clustering_vare$cluster # obtain group assignment

# plot groups in the  multivariate space of the pca
ordiplot(pca_vare, display = "si", type = "n")
for(i in 1:4) ordihull(pca_vare, groups = groups_vare, show.groups = i, col = i)
for(i in 1:4) ordispider(pca_vare, groups = groups_vare, show.groups = i, col = i)

indi_vare <- multipatt(vare_pa, groups_vare)
summary(indi_vare)

# do it all for varechem
# presence-absence matrix
chem_pa <- decostand(varechem, method = "pa") # decostand with method pa transforms values into a presence/absence scale (0/1)

hist(as.matrix((varechem))) # double-zero issue
hist(as.matrix((varechem[varechem > 0]))) # right skewed cover data

chem_log <- log1p(varechem) # log transform
chem_sqrt <- sqrt(varechem) # transform using square root

pca_chem <- rda(chem_log)
pca_chem

# ordination biplots
ordiplot(pca_chem, display = "si", type= "text")
# ordiplot(pca_chem, display = "sp", type = "text")


str(pca_chem) # structure of pca object

pca_chem$tot.chi
(pca_chem$CA$eig[1] / pca_chem$tot.chi) * 100 # explained variance by first axis

((pca_chem$CA$eig[1] + pca_chem$CA$eig[2]) / pca_chem$tot.chi) * 100 # explained variance by first and second axes

summary(pca_chem)

### Clustering
# plot partition according their species composition, using K-means clustering
set.seed(1221) # fix the random factor
clustering_chem <- kmeans(chem_log, 4)
clustering_chem
?kmeans

groups_chem <- clustering_chem$cluster # obtain group assignment

# plot groups in the  multivariate space of the pca
ordiplot(pca_chem, display = "si", type = "n")
for(i in 1:4) ordihull(pca_chem, groups = groups_chem, show.groups = i, col = i)
for(i in 1:4) ordispider(pca_chem, groups = groups_chem, show.groups = i, col = i)

indi_chem <- multipatt(chem_pa, groups_chem)
summary(indi_chem)
?multipatt
