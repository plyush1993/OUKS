##############################################################################################################################################################
# Table of contents
##############################################################################################################################################################

# Installation
# Base example
# Generate table for experiment with random order repeats and duplicates
# Generate table for experiment for batches with random order, repeats and duplicates
# Generate table for experiment for batches with random order, repeats and serial duplicates
# Generate table for experiment for batches with random order, serial duplicates and repeats
# "randomizr" package
# References

##############################################################################################################################################################
# Installation
##############################################################################################################################################################

# setup environment
setwd("D:/")

##############################################################################################################################################################
# Base example
##############################################################################################################################################################

set.seed(1234)
numbers <- seq(from = 1, to = 18, by = 1) # generate vector with length as total number of trials (samples)
sample(numbers, size = length(numbers), replace = F)

repeats <- rep(seq(from = 1, to = 18, by = 1),2) # generate repeats vector with length as total number of trials (samples)
sample(repeats, size = length(repeats), replace = F)

##############################################################################################################################################################
# Generate table for experiment with random order repeats and duplicates
##############################################################################################################################################################

n <- 1 # number of repeats (biological repeats)
t <- 2 # number of duplicates (technical repeats)
f_id <- 1 # first sample ID
l_id <- 124 # last sample ID
repeats <- rep(seq(from = f_id, to = l_id, by = 1),t) # generate repeats vector with length as number of samples*tech rep
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = n, replace = F)
df_repeats <- list()
df_repeats_df <- as.data.frame(1:length(repeats))

for (i in 1:n) {
 set.seed(seed[i])
 df_repeats[[i]] <- sample(repeats, size = length(repeats), replace = F)
 df_repeats_df <- cbind(df_repeats_df,df_repeats[[i]])
 }

vn <- c()
for (i in 1:n) {
  vn[i] <- paste("biorepeat", i)
}

vn <- c("id", vn)
colnames(df_repeats_df) <- vn
df_repeats_df

##############################################################################################################################################################
# Generate table for experiment for batches with random order, repeats and duplicates
##############################################################################################################################################################

library(tuple)
n <- 1 # number of repeats (biological repeats)
t <- 2 # number of duplicates (technical repeats)
b <- 10 # number of samples per batch
f_id <- 1 # first sample ID
l_id <- 124 # last sample ID
repeats <- rep(seq(from = f_id, to = l_id, by = 1),n) # generate repeats vector with length as  number of samples*tot rep
set.seed(1234) # seed for reproduce
repeats_random <- sample(repeats, size = length(repeats), replace = F)
m <- matrix(repeats_random, nrow = b)
cm <- c(m)
cm[which(tuplicated(cm, n+1), T)] <- NA
m2 <- matrix(cm, nrow = b)
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = ncol(m2), replace = F)
df_rep_batch <- list()
df_rep_batch_df <- as.data.frame(1:(length(m2[,1])*t))

for (i in 1:ncol(m2)) {
  set.seed(seed[i])
  df_rep_batch[[i]] <- sample(rep(m2[,i],t), size = length(m2[,i])*t, replace = F)
  df_rep_batch_df <- cbind(df_rep_batch_df,df_rep_batch[[i]])
}

vn <- c()
for (i in 1:(ncol(df_rep_batch_df)-1)) {
  vn[i] <- paste("batch", i)
}

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn
df_rep_batch_df

##############################################################################################################################################################
# Generate table for experiment for batches with random order, repeats and serial duplicates
##############################################################################################################################################################

library(tuple)
n <- 1 # number of repeats (biological repeats)
t <- 2 # number of duplicates (technical repeats)
b <- 10 # number of samples per batch
f_id <- 1 # first sample ID
l_id <- 124 # last sample ID
repeats <- rep(seq(from = f_id, to = l_id, by = 1),n) # generate repeats vector with length as  number of samples*tot rep
set.seed(1234) # seed for reproduce
repeats_random <- sample(repeats, size = length(repeats), replace = F)
m <- matrix(repeats_random, nrow = b)
cm <- c(m)
cm[which(tuplicated(cm, n+1), T)] <- NA
m2 <- matrix(cm, nrow = b)
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = ncol(m2), replace = F)
df_rep_batch <- list()
df_rep_batch_df <- as.data.frame(1:(length(m2[,1])*t))

for (i in 1:ncol(m2)) {
  set.seed(seed[i])
  df_rep_batch[[i]] <- sample(m2[,i], size = length(m2[,i]), replace = F)
  df_rep_batch[[i]] <- c(df_rep_batch[[i]],df_rep_batch[[i]])
  df_rep_batch_df <- cbind(df_rep_batch_df,df_rep_batch[[i]])
}

vn <- c()
for (i in 1:(ncol(df_rep_batch_df)-1)) {
  vn[i] <- paste("batch", i)
}

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn
df_rep_batch_df

##############################################################################################################################################################
# Generate table for experiment for batches with random order, serial duplicates and repeats
##############################################################################################################################################################

library(tuple)
n <- 1 # number of repeats (biological repeats)
t <- 1 # number of duplicates (technical repeats)
b <- 75 # number of samples per batch
f_id <- 1 # first sample ID
l_id <- 75 # last sample ID
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = n, replace = F)
repeats <- list()
for (i in 1:n) {
  set.seed(seed[i])
  repeats[[i]] <- paste(sample(f_id:l_id, size = l_id, replace = F),i)
}

m <- matrix(unlist(repeats), nrow = b)
cm <- c(m)
cm[which(tuplicated(cm, 2), T)] <- NA
m2 <- matrix(cm, nrow = b)
mm <- replicate(t, m2)
mm <- lapply(1:dim(mm)[3], function(y) as.matrix(mm[,,y]))
df_rep_batch_df <- do.call(rbind, mm)
df_rep_batch_df <- as.data.frame(cbind(1:nrow(df_rep_batch_df), df_rep_batch_df))

vn <- c()
for (i in 1:(ncol(df_rep_batch_df)-1)) {
  vn[i] <- paste("batch", i)
}

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn
df_rep_batch_df

##############################################################################################################################################################
# "randomizr" package
##############################################################################################################################################################

library(randomizr)

n <- 124 # total number of samples
set.seed(12) # seed for reproduce
factors <- c("CG", "TG") # list of factors/groups
prob <- c(0.5,0.5) # list of probabilities of assignment to factors/groups (sum is equal to 1)
ed <- simple_ra(N = n, conditions = factors, prob_each = prob)
table(ed)
ed

##############################################################################################################################################################
# References
##############################################################################################################################################################

# 1. Pezzatti, Julian, et al. "Implementation of liquid chromatography-high resolution mass spectrometry methods for untargeted metabolomic analyses of biological samples: A tutorial." Analytica Chimica Acta 1105 (2020): 28-44.
# 2. Dudzik, Danuta, et al. "Quality assurance procedures for mass spectrometry untargeted metabolomics. a review." Journal of pharmaceutical and biomedical analysis 147 (2018): 149-173.