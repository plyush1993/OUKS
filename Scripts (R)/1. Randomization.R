##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                              TABLE OF CONTENTS                           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Installation
# Base example
# Generate table for experiment with random order repeats and duplicates
# Generate table for experiment for batches with random order, repeats and duplicates
# Generate table for experiment for batches with random order, repeats and serial duplicates
# Generate table for experiment for batches with random order, serial duplicates and repeats
# Add specific sample type to table for experiment for batches
# Add parallel or sequential copy for factor
# Reshape, add names
# Generate table for experiment in tidy ready-to-use format
# "randomizr" package
# References

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Installation                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# setup environment
setwd("D:/")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Base example                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(1234)
numbers <- seq(from = 1, to = 18, by = 1) # generate vector with length as total number of trials (samples)
sample(numbers, size = length(numbers), replace = F)

repeats <- rep(seq(from = 1, to = 18, by = 1),2) # generate repeats vector with length as total number of trials (samples)
sample(repeats, size = length(repeats), replace = F)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Generate table for experiment with random order repeats and duplicates ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n <- 1 # number of repeats (biological repeats)
t <- 2 # number of duplicates (technical repeats)
f_id <- 1 # first sample ID
l_id <- 124 # last sample ID
repeats <- rep(seq(from = f_id, to = l_id, by = 1),t) # generate repeats vector with length as number of samples*tech rep
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = n, replace = F)
df_repeats <- list()
df_repeats_df <- as.data.frame(1:length(repeats))

for (i in 1:n) ({
 set.seed(seed[i])
 df_repeats[[i]] <- sample(repeats, size = length(repeats), replace = F)
 df_repeats_df <- cbind(df_repeats_df,df_repeats[[i]])
 })

vn <- c()
for (i in 1:n) ({
  vn[i] <- paste("biorepeat", i)
})

vn <- c("id", vn)
colnames(df_repeats_df) <- vn
df_repeats_df

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate table for experiment for batches with random order, repeats and duplicate----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

for (i in 1:ncol(m2)) ({
  set.seed(seed[i])
  df_rep_batch[[i]] <- sample(rep(m2[,i],t), size = length(m2[,i])*t, replace = F)
  df_rep_batch_df <- cbind(df_rep_batch_df,df_rep_batch[[i]])
})

vn <- c()
for (i in 1:(ncol(df_rep_batch_df)-1)) ({
  vn[i] <- paste("batch", i)
})

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn
df_rep_batch_df

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate table for experiment for batches with random order, repeats and serial duplicate----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
for (i in 1:(ncol(df_rep_batch_df)-1)) ({
  vn[i] <- paste("batch", i)
})

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn
df_rep_batch_df

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate table for experiment for batches with random order, serial duplicates and repeat----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tuple)
library(stringr)
n <- 2 # number of repeats (biological repeats)
t <- 3 # number of duplicates (technical repeats)
b <- 75 # number of samples per batch
f_id <- 1 # first sample ID
l_id <- 75 # last sample ID
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = n, replace = F)
repeats <- list()
for (i in 1:n) {
  set.seed(seed[i])
  repeats[[i]] <- paste0(sample(f_id:l_id, size = l_id, replace = F),"_",i)
}

m <- matrix(unlist(repeats), nrow = b)
cm <- c(m)
cm[which(tuplicated(cm, 2), T)] <- NA
cm <- as.data.frame(str_split(cm, "_"))[1,]
m2 <- matrix(cm, nrow = b)
mm <- replicate(t, m2)
mm <- lapply(1:dim(mm)[3], function(y) as.matrix(mm[,,y]))
df_rep_batch_df <- do.call(rbind, mm)
df_rep_batch_df <- as.data.frame(cbind(1:nrow(df_rep_batch_df), df_rep_batch_df))
df_rep_batch_df <- data.frame(lapply(df_rep_batch_df, unlist))

vn <- c()
for (i in 1:(ncol(df_rep_batch_df)-1)) {
  vn[i] <- paste("batch", i)
}

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn
df_rep_batch_df

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##        Add specific sample type to table for experiment for batches      ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note:
# it is possible to use it sequentially for several sample types (after each round just do: df_rep_batch_df <- df)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

type <- "QC" # name of sample type
freq <- 10 # frequency of specific sample in table (every "freq" element)
ind <- seq(1, nrow(df_rep_batch_df), freq)-1 # index for insert
add_sample <- as.data.frame(paste0(type, "_", 1:length(ind)))
add_s <- as.data.frame(cbind(ind, matrix(NA, nrow = nrow(add_sample), ncol = ncol(df_rep_batch_df)-1)))
add_s[,-1] <- add_sample 
colnames(add_s) <- colnames(df_rep_batch_df)
add_s[,-1] <- paste0(type, "-", 1:(length(ind)*ncol(add_s[,-1]))) # add unique names for specific sample type

df <- as.data.frame(rbind(df_rep_batch_df, add_s))
df <- data.frame(lapply(df, unlist))
df$id <- as.numeric(df$id)
df <- df[order(df$id, decreasing = F),]
df$id <- c(1:nrow(df))
df

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                 Add parallel or sequential copy for factor               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parallel (for polarity)
df_c <- as.data.frame(sapply(2:ncol(df), function(x) rep(df[,x], each = 2)))
df_c <- as.data.frame(cbind(ind = 1:nrow(df_c), pol = c("pos", "neg"), df_c)) # names of factor
colnames(df_c)[-c(1:2)] <- c(colnames(df))[-1] 

# Sequential (for mode)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note:
# it is possible to use it sequentially after "Parallel" (just do: df <- df_c) or instead "Parallel"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_c <- as.data.frame(sapply(2:ncol(df), function(x) rep(df[,x], times = 2)))
df_c <- as.data.frame(cbind(ind = 1:nrow(df_c), mode = rep(c("RP", "HILIC"), each = nrow(df_c)/2), df_c)) # names of factor
colnames(df_c)[-c(1:2)] <- c(colnames(df))[-1] 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             Reshape, add names                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(stringr)

# reshape
df_r <- stack(df, select = -c(1)) # in "select" add column which should be removed during reshape

# create and add names
df_r <- as.data.frame(paste0(1:nrow(df_r), "._", df_r[,2], "_", df_r[,1])) # add run order and batch index
colnames(df_r) <- "Sample"
df_r <- as.data.frame(paste0("Name:_", df_r$Sample)) # add at the beginning
colnames(df_r) <- "Sample"
df_r <- as.data.frame(paste0(Sample = df_r$Sample, "_.mzXML")) # add at the end
colnames(df_r) <- "Sample"

# create table with info
df_str <- as.data.frame(t(sapply(1:nrow(df_r), function(y) unlist(str_split(df_r[y,], "_")))))[,-c(1,5)]
colnames(df_str) <- c("Order", "Batch", "Name") 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##          Generate table for experiment in tidy ready-to-use format       ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tuple)
library(stringr)
n <- 5 # number of repeats (biological repeats)
t <- 3 # number of duplicates (technical repeats)
b <- 25 # number of samples per batch
f_id <- 1 # first sample ID
l_id <- 73 # last sample ID
set.seed(1234) # seed for reproduce
seed <-  sample(x = 1:12345, size = n, replace = F)
repeats <- list()
for (i in 1:n) {
  set.seed(seed[i])
  repeats[[i]] <- paste0(sample(f_id:l_id, size = l_id, replace = F),"_bio-",i) # add index biological repeats
}

m <- matrix(unlist(repeats), nrow = b)
cm <- c(m)
cm[which(tuplicated(cm, 2), T)] <- NA
m2 <- matrix(cm, nrow = b)
mm <- replicate(t, m2)
mm <- lapply(1:dim(mm)[3], function(y) as.matrix(mm[,,y]))
df_rep_batch_df <- do.call(rbind, mm)
df_rep_batch_df <- as.data.frame(cbind(1:nrow(df_rep_batch_df), df_rep_batch_df))
df_rep_batch_df <- data.frame(lapply(df_rep_batch_df, unlist))

vn <- c()
for (i in 1:(ncol(df_rep_batch_df)-1)) {
  vn[i] <- paste("batch", i)
}

vn <- c("id", vn)
colnames(df_rep_batch_df) <- vn

df_rep_batch_df[,-1] <- sapply(2:ncol(df_rep_batch_df), function(x) paste0(df_rep_batch_df[,x], "_tech-", rep(1:t, each = b))) # add index technical repeats

# add specific sample type
# if no specific sample type => do: df <- df_rep_batch_df
type <- "QC" # name of sample type
freq <- 10 # frequency of specific sample in table (every "freq" element)
ind <- seq(1, nrow(df_rep_batch_df), freq)-1 # index for insert
add_sample <- as.data.frame(paste0(type, "_", 1:length(ind)))
add_s <- as.data.frame(cbind(ind, matrix(NA, nrow = nrow(add_sample), ncol = ncol(df_rep_batch_df)-1)))
add_s[,-1] <- add_sample 
colnames(add_s) <- colnames(df_rep_batch_df)
add_s[,-1] <- paste0(type, "-", 1:(length(ind)*ncol(add_s[,-1]))) # add unique names for specific sample type

df <- as.data.frame(rbind(df_rep_batch_df, add_s))
df <- data.frame(lapply(df, unlist))
df$id <- as.numeric(df$id)
df <- df[order(df$id, decreasing = F),]
df$id <- c(1:nrow(df))

# add factor
# if no factor => do: df_comb <- df
df_c <- as.data.frame(sapply(2:ncol(df), function(x) rep(df[,x], each = 2)))
df_c <- as.data.frame(cbind(ind = 1:nrow(df_c), pol = c("pos", "neg"), df_c)) # names of factor
colnames(df_c)[-c(1:2)] <- c(colnames(df))[-1] 

df_comb <- df_c
df_comb[,-c(1:2)] <- sapply(3:ncol(df_comb), function(x) paste0(df_comb[,x], "_pol-", df_comb[,2])) # names of factor

# reshape and add names
df_long <- stack(df_comb, select = -c(1:2)) # in "select" add column which should be removed during reshape
df_long <- as.data.frame(paste0(1:nrow(df_long), "_order", "_", df_long[,2], "_", df_long[,1]))
colnames(df_long) <- "Sample"

# table with info
df_str <- sapply(1:nrow(df_long), function(y) unlist(str_split(df_long[y,], "_")))
ind_str <- as.numeric(sapply(df_str, length))
df_str_add_sample <- do.call(rbind.data.frame, df_str[ind_str == min(unique(ind_str))])[,-2]
colnames(df_str_add_sample) <- c("Order", "Batch", "Name") 
df_str_exp_sample <- do.call(rbind.data.frame, df_str[ind_str != min(unique(ind_str))])[,-2]
colnames(df_str_exp_sample) <- c("Order", "Batch", "Name", "Bio", "Tech")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                            "randomizr" package                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(randomizr)

n <- 124 # total number of samples
set.seed(12) # seed for reproduce
factors <- c("CG", "TG") # list of factors/groups
prob <- c(0.5,0.5) # list of probabilities of assignment to factors/groups (sum is equal to 1)
ed <- simple_ra(N = n, conditions = factors, prob_each = prob)
table(ed)
ed

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 References                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Pezzatti, Julian, et al. "Implementation of liquid chromatography-high resolution mass spectrometry methods for untargeted metabolomic analyses of biological samples: A tutorial." Analytica Chimica Acta 1105 (2020): 28-44.
# 2. Dudzik, Danuta, et al. "Quality assurance procedures for mass spectrometry untargeted metabolomics. a review." Journal of pharmaceutical and biomedical analysis 147 (2018): 149-173.