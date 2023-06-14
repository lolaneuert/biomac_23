### CALCULATOR
###########################
3 + 3 # Simple sum

# assignment of a variable
s <- 3 + 3
s # print the variable

# assign & print
(s <- 5 + 6)

p <- 2 * s

# list of variables
calc <- list()
calc$s
calc$a <- 5
calc$p <- 2 * (3 + 3)

# several rules to call specific variables in a list
calc[1]
calc[1:3]
calc["s"]
calc[c(1, 3)]
calc$s

# variables can be different in the length - i.e, two or more values
calc$s[2] <- 7 + 23
calc$s

calc$result <- calc$s[1] + calc$s[2]
calc

### COMMUNITY MATRIX: vegan::dune CO-OCCURENCE DATA
###########################
# 
library(vegan)
library(tidyverse)
?vegan

# load community and environmental data (vegan package)
data(dune)
data(dune.env)

# dimensions (rows, columns, cells)
# in R: data[row(s), column(s)]
dim(dune) # two dimensions for data.frame objects
class(dune)
row.names(dune)
colnames(dune)
dune[5:7, 5:9]

head(dune)
tail(dune)
dune[3:5, 3:7] > 2
dune[1:5, 1:5] > 0

# we apply a condition, which can be true or false
rowSums(dune[1:5, 1:5] > 0) 

# species richness: summing all occuring species (cover > 0) per observations
rich <- rowSums(dune > 0)
rich
class(rich)

# tidyverse library is used for the data management
# code's grammar is a bit different from base R code
library(tidyverse)
# ctrl + shift + M the shortcut for the pipe operator of tidyverse: %>% 

# Using a pipe-line of tidyverse
a <- rich %>% 
  sort() %>% 
  as.data.frame() %>% 
  rename(richness = ".")

# Translating in base R, which is more complex to read
a <- sort(rich)
a <- as.data.frame(a)
colnames(a) <- "richness"
a


