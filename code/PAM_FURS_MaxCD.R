# Laura Tibbs-Cortes
# Code for PAM, FURS, and MaxCD optimal training set methods implemented for inbreds.
# Code based on paper https://doi.org/10.1016/j.molp.2018.12.022
# and personal communication with Tingting Guo.

# load libraries
# If any libraries are not installed, first run: install.packages("libraryname")
library(data.table)
library(tidyverse)
library("cluster")
library(dendextend)

# Set your working directory to the "code" folder using setwd()

num.select <- 30 # set the number of individuals to select for the training set

# read in kinship matrix of the population that the training set will be selected from (in this case, BGEM):
K<- as.matrix(fread("../data/BGEM.kin.csv"), rownames="Taxa")

# PAM (partitioning around medoids) ---------------------------------------------------------------------

# In inbreds, PAM provides a unique solution for a given training set size in a given population

# group taxa into desired number of clusters and identify each cluster's medoid:
c.mat <- 1-cov2cor(K)
clust <- pam(x=c.mat,k=num.select,diss=TRUE) 

train.taxa.PAM <- clust$medoids # save training set 

train.taxa.PAM # view training set (unique solution)


# FURS (fast and unique representative subset selection) --------------------------------------------------------------------

FURS <- function(adajac_mat, trn_size,
                 thres = FALSE,thres_t=NULL){
  label_list <- colnames(adajac_mat)
  diag(adajac_mat) <- 0 # set diagonals to 0
  rep_set <- NULL # make empty training sample
  deg_cent <- sort(rowSums(adajac_mat),decreasing = TRUE) # sort row sums of adajac_mat
  cent_set <- data.frame(label=names(deg_cent),deg=deg_cent,active_status=TRUE) # each taxa's row sum
  cent_set$label <- as.character(cent_set$label)
  if(thres){ # by default, thres = F, so this is skipped
    if(is.null(thres_t)){ # by default, thres_t=null, so this is calculated
      thres_t <- median(cent_set$deg) # set median row sum as thres_t
    }
    cent_set <- subset(cent_set,deg>min(thres_t,median(cent_set$deg))) # remove lines with row sums below threshold
  }
  
  # Added by Laura Tibbs-Cortes May 26, 2021 to randomize starting point for training set among equally-good candidates
  # This enables identification of multiple training sets that are expected to be equally good
  cent_set <- as_tibble(cent_set) %>%
    group_by(deg) %>%
    slice(sample(1:n())) %>%
    arrange(-deg)
  
  while(length(rep_set)<trn_size){ # continue until rep_set reaches provided trn_size (=desired k)
    if(sum(cent_set$active_status==TRUE)==0){ # if there are no "active" lines, start over with all re-activated
      cent_set$active_status=TRUE
      # print("!!!")
    }
    select_label <- cent_set$label[cent_set$active_status==TRUE][1] # pull first line (max deg); but there are ties!
    rep_set <- c(rep_set,select_label)
    neigh <- label_list[which(adajac_mat[which(label_list==select_label),]==1)] # pull "neighbors" of the selected line: those with adajac_mat entries =1 with this line
    cent_set$active_status[cent_set$label %in% neigh] <- FALSE # set neightbors to inactive
    cent_set <- cent_set[-which(cent_set$label==select_label),] # remove the selected sample from the list of candidates
  }
  return(rep_set)
}

my_adajac_mat <- cov2cor(K)

# "tune" the matrix based on threshold 
tune=0.2 # see https://doi.org/10.1016/j.molp.2018.12.022 -- this is threshold to declare linked vs unlinked
my_adajac_mat[my_adajac_mat>tune] = 1
my_adajac_mat[my_adajac_mat<tune] = 0

train.taxa.FURS <- FURS(my_adajac_mat,num.select, FALSE) # get FURS optimal training set (will be different each run)

train.taxa.FURS #  view training set (non-unique solution)


# MaxCD (maximization of connectedness and diversity) -------------------------------------------------------------------

dist <- dist(K) # Calculate distance matrix from kinship
tree <- hclust(dist) # hierarchical cluster analysis on distance matrix

new.tree <- dendextend::shuffle(tree) # randomly shuffle the tree about its nodes. Enables identification of multiple training sets expected to be equally good.
tree.order <- new.tree$labels[new.tree$order] # get order of taxa in tree

# select every kth line to get up to num.select in training set:
train.taxa.MaxCD <- tree.order[floor(c(1:num.select)*(length(tree.order)/num.select))]
stopifnot(length(train.taxa.MaxCD)==num.select)

train.taxa.MaxCD #  view training set (non-unique solution)

