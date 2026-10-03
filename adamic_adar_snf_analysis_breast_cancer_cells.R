source("/home/ojinnagc/SNFtool/R/affinityMatrix.R")
source("/home/ojinnagc/SNFtool/R/chiDist2.R")
source("/home/ojinnagc/SNFtool/R/dist2.R")
source("/home/ojinnagc/SNFtool/R/SNF.R")
source("/home/ojinnagc/SNFtool/R/spectralClustering.r")
source("/home/ojinnagc/SNFtool/R/standardNormalization.R")
source("/home/ojinnagc/SNFtool/R/internal.R")
exists(".dominateset")
#getAnywhere(".dominateset")
#ls()
#readLines("gse240567_var.csv", n = 3)

gene <- read.csv(
  "breast_cancer_gene_var.csv",
  row.names = 1,
  check.names = FALSE
)

dim(gene)
head(colnames(gene))
head(rownames(gene))

meth <- read.csv(
  "breast_cancer_methyl_var.csv",
  row.names = 1,
  check.names = FALSE
)

dim(meth)
head(rownames(meth))
head(colnames(meth))

## To be safe explicitly convert column names to charcters/strings

colnames(gene) <- as.character(colnames(gene))
colnames(meth) <- as.character(colnames(meth))

# check
class(colnames(meth))
str(colnames(meth)[1:5])

# check identical patient ordering
identical(colnames(gene), colnames(meth))

# Then transpose

gene <- t(gene)
meth <- t(meth)

dim(gene)
dim(meth)

rownames(gene)[1:5]
rownames(meth)[1:5]

gene_norm <- standardNormalization(gene)
meth_norm <- standardNormalization(meth)

dim(gene_norm)
dim(meth_norm)

gene_norm[1:5,1:5]
meth_norm[1:5,1:5]

## Compute distance matrix
gene_dist <- dist2(gene_norm, gene_norm)
meth_dist <- dist2(meth_norm, meth_norm)

dim(gene_dist)
dim(meth_dist)

range(gene_dist)
range(meth_dist)
## Verify distance matrix properties; diagonal should be 0
gene_dist[1:5,1:5]
meth_dist[1:5,1:5]

saveRDS(gene_dist, "breast_cancer_gene_dist.rds")
saveRDS(meth_dist, "breast_cancer_meth_dist.rds")
list.files(pattern = "rds")

object.size(gene_dist)
object.size(meth_dist)

## Compute affinity matrix
gene_W <- affinityMatrix(
  gene_dist,
  K = 5,
  sigma = 0.5
)

meth_W <- affinityMatrix(
  meth_dist,
  K = 5,
  sigma = 0.5
)
saveRDS(gene_W, "breast_cancer_gene_W_affinity.rds")
saveRDS(meth_W, "breast_cancer_meth_W_affinity.rds")

dim(gene_W)
dim(meth_W)

range(gene_W)
range(meth_W)

gene_W[1:5,1:5]
meth_W[1:5,1:5]

## SNF implementation: here transition probability graph (Normalized affinity matrix) P, 
## Local KNN transition matrix S, and Diffusion are all computed.

# Step 1: Create wall
Wall <- list(gene_W, meth_W)

length(Wall)

# Step 2: Check patient names
#check_wall_names(Wall)

# Verify manually:
identical(
  rownames(gene_W),
  rownames(meth_W)
)

identical(
  colnames(gene_W),
  colnames(meth_W)
)

# Step 3: Run SNF 

W_fused <- SNF(
  Wall,
  K = 5,
  t = 30
)
# Step 4: Save

saveRDS(
  W_fused,
  "breast_cancer_W_fused.rds"
)

# Load the fused SNF network from disk (loaded back for association analysis 
# way down below to prevent starting from the very beginning to get w fused)

W_fused <- readRDS(
  "breast_cancer_W_fused.rds"
)
# Check dimensions
#dim(W_fused)

# View the first few patient IDs
head(rownames(W_fused))

# Step 5: Inspect result

dim(W_fused)

range(W_fused)

W_fused[1:5,1:5]

# then inspect dim(W_fused),range(W_fused) and W_fused[1:5,1:5] you just ran to 
# determine the number of patient cluster C using eigengap analysis rather than
# choosing a cluster number blindly

# eigengap analysis for determining optimum K
D <- diag(rowSums(W_fused))

L <- D - W_fused

D_inv_sqrt <- diag(
  1/sqrt(diag(D))
)

L_norm <- D_inv_sqrt %*% L %*% D_inv_sqrt

eig <- eigen(
  L_norm,
  symmetric = TRUE,
  only.values = TRUE
)

vals <- sort(eig$values)

vals[1:20]
# get the clusters

# Eigengaps
gaps <- diff(vals[1:20])

# Find maximum gap
which.max(gaps)
max(gaps)

# examine the non-trivial eigengaps

# Non-trivial eigengaps for K = 2 through K = 10
for (k in 2:10) {
  cat(
    "K =", k,
    " | eigengap =", vals[k + 1] - vals[k],
    "\n"
  )
}

# PLot the eigen
plot(
  1:20,
  vals[1:20],
  type = "b",
  pch = 19,
  main = "Eigengap Analysis – Standard SNF",
  xlab = "Eigenvalue Index",
  ylab = "Eigenvalue"
)

plot(
  2:20,
  vals[2:20],
  type = "b",
  pch = 19,
  main = "Non-trivial Eigengap Analysis – Standard SNF",
  xlab = "Eigenvalue Index",
  ylab = "Eigenvalue"
)
library(SNFtool)

estimate_snf <- estimateNumberOfClustersGivenGraph(
  W_fused,
  NUMC = 2:10
)

estimate_snf

## IMPORT THE SUBTYPE FILE

cancer_subtypes <- read.csv("breast_cancer_subtypes_179.csv")

head(cancer_subtypes)
dim(cancer_subtypes)
table(cancer_subtypes$hu_subtype)

#verify that patient order match

rownames(W_fused)[1:10]
cancer_subtypes$TAX[1:10]
identical(rownames(W_fused), cancer_subtypes$TAX)

# also check dimension
nrow(W_fused) == nrow(cancer_subtypes)

# Spectral clustering of the baseline SNF fused matrix
set.seed(123)
snf_cluster_k2 <- spectralClustering(W_fused, 2)

table(snf_cluster_k2)

set.seed(123)
snf_cluster_k3 <- spectralClustering(W_fused, 3)

table(snf_cluster_k3)

set.seed(123)
snf_cluster_k4 <- spectralClustering(W_fused, 4)

table(snf_cluster_k4)

set.seed(123)
snf_cluster_k5 <- spectralClustering(W_fused, 5)

table(snf_cluster_k5)

set.seed(123)
snf_cluster_k6 <- spectralClustering(W_fused, 6)

table(snf_cluster_k6)

# First compare K = 6 with the known subtypes
# We have six reference groups also, so make a simple cross tabulation

# Rows = clusters discovered by SNF
# Columns = Known breast-cancer subtypes
# Each number = how many patients of that subtype ended up in that SNF cluster

table(
  SNF_cluster = snf_cluster_k6,
  Subtype = cancer_subtypes$hu_subtype
)

## Now to cross-tabulation comparison for K = 2 to 5
table(
  SNF_cluster = snf_cluster_k2,
  Subtype = cancer_subtypes$hu_subtype
)
table(
  SNF_cluster = snf_cluster_k3,
  Subtype = cancer_subtypes$hu_subtype
)
table(
  SNF_cluster = snf_cluster_k4,
  Subtype = cancer_subtypes$hu_subtype
)
table(
  SNF_cluster = snf_cluster_k5,
  Subtype = cancer_subtypes$hu_subtype
)

## Calculating ARI and NMI for K = 2 to 6
# Make sure labeled are aligned - create subtype vector explicitly
# from the 179-patient file

subtype <- cancer_subtypes$hu_subtype

length(subtype)
table(subtype)

# Load packages for ARI and NMI

library(aricode)
# Calculate ARI and NMI for K = 2
# First check vectors for missing values

sum(is.na(snf_cluster_k2))
sum(is.na(subtype))

table(subtype, useNA = "ifany")
table(snf_cluster_k2, useNA = "ifany")

# Check data types
class(snf_cluster_k2)
class(subtype)

str(snf_cluster_k2)
str(subtype)

# Convert both to factors
snf_cluster_k2_factor <- factor(snf_cluster_k2)
subtype_factor <- factor(subtype)

class(snf_cluster_k2_factor)
class(subtype_factor)

# Convert the subtype labels to numerical categories to enable ARI and NMI run

subtype_num <- as.integer(factor(subtype))

table(subtype, subtype_num)

# Verify we still have 179 patients

length(subtype_num)
sum(is.na(subtype_num))

# Now run ARI/NMI

ARI_k2 <- ARI(snf_cluster_k2, subtype_num)
NMI_k2 <- NMI(snf_cluster_k2, subtype_num)

ARI_k2
NMI_k2

ARI_k3 <- ARI(snf_cluster_k3, subtype_num)
NMI_k3 <- NMI(snf_cluster_k3, subtype_num)

ARI_k3
NMI_k3

ARI_k4 <- ARI(snf_cluster_k4, subtype_num)
NMI_k4 <- NMI(snf_cluster_k4, subtype_num)

ARI_k4
NMI_k4

ARI_k5 <- ARI(snf_cluster_k5, subtype_num)
NMI_k5 <- NMI(snf_cluster_k5, subtype_num)

ARI_k5
NMI_k5

ARI_k6 <- ARI(snf_cluster_k6, subtype_num)
NMI_k6 <- NMI(snf_cluster_k6, subtype_num)

ARI_k6
NMI_k6

# Display as dataframe
data.frame(
  K = 2:6,
  ARI = c(ARI_k2, ARI_k3, ARI_k4, ARI_k5, ARI_k6),
  NMI = c(NMI_k2, NMI_k3, NMI_k4, NMI_k5, NMI_k6)
)

###############################################################################

### PROJECT 2: ADAMIC_ADAR

gene_W <- readRDS("breast_cancer_gene_W_affinity.rds")
meth_W <- readRDS("breast_cancer_meth_W_affinity.rds")

# Calculate Adamic Adar for gene
dim(gene_W)
class(gene_W)
gene_W[1:5, 1:5]

## Step 1 - Create the empty graph matrix

K <- 5

n <- nrow(gene_W)

adjacency <- matrix(
  0,
  nrow = n,
  ncol = n
)

for(i in 1:n)
{
  # Rank all patients from most similar to least similar
  neighbors <- order(
    gene_W[i, ],
    decreasing = TRUE
  )
  
  # Remove the patient itself
  neighbors <- neighbors[neighbors != i]
  
  # Keep only the K nearest neighbours
  neighbors <- neighbors[1:K]
  
  # Add edges to the graph
  adjacency[i, neighbors] <- 1
}

# Step 3 — Symmetrize the graph (i.e., make graph undirected)
adjacency_sym <- ifelse(
  adjacency + t(adjacency) > 0,
  1,
  0
)
# Step 4 — Verify that the graph is now symmetric
max(abs(adjacency_sym - t(adjacency_sym)))

# Step 5 — Let's also see how many edges were added
sum(adjacency)
sum(adjacency_sym)

## Now compute adamic adar from scratch

# Compute the degree of every patient
degree <- rowSums(adjacency_sym)
summary(degree)
table(degree)

## Now implement adamic adar for every patient
# Step 1 - Store the neighbors of every patient
neighbors_list <- lapply(
  1:n,
  function(i) which(adjacency_sym[i, ] == 1)
)
# Step 2 — Create an empty Adamic–Adar matrix
AA <- matrix(
  0,
  nrow = n,
  ncol = n
)
# Step 3 — Calculate Adamic–Adar for every patient pair
for(i in 1:(n-1))
{
  for(j in (i + 1):n)
  {
    common_neighbors <- intersect(
      neighbors_list[[i]],
      neighbors_list[[j]]
    )
    
    if(length(common_neighbors) > 0)
    {
      AA[i, j] <- sum(
        1 / log(degree[common_neighbors])
      )
    }
    
    AA[j, i] <- AA[i, j]
  }
}

# Calculate Adamic Adar for methylation

# Step 1 — Create the K-nearest-neighbor graph

K <- 5

n <- nrow(meth_W)

adjacency_meth <- matrix(
  0,
  nrow = n,
  ncol = n
)

for(i in 1:n)
{
  # Rank all patients from most similar to least similar
  neighbors <- order(
    meth_W[i, ],
    decreasing = TRUE
  )
  
  # Remove the patient itself
  neighbors <- neighbors[neighbors != i]
  
  # Keep only the top K neighbors
  neighbors <- neighbors[1:K]
  
  # Add edges
  adjacency_meth[i, neighbors] <- 1
}

# Step 2 — Make the graph undirected
adjacency_meth_sym <- ifelse(
  adjacency_meth + t(adjacency_meth) > 0,
  1,
  0
)
# Step 3 — Calculate patient degrees
degree_meth <- rowSums(adjacency_meth_sym)

# Step 4 — Store every patient's neighbors

neighbors_list_meth <- lapply(
  1:n,
  function(i) which(adjacency_meth_sym[i, ] == 1)
)
# Step 5 — Create the empty Adamic–Adar matrix
AA_meth <- matrix(
  0,
  nrow = n,
  ncol = n
)
# Step 6 — Calculate Adamic–Adar for every patient pair
for(i in 1:(n-1))
{
  for(j in (i + 1):n)
  {
    common_neighbors <- intersect(
      neighbors_list_meth[[i]],
      neighbors_list_meth[[j]]
    )
    
    if(length(common_neighbors) > 0)
    {
      AA_meth[i, j] <- sum(
        1 / log(degree_meth[common_neighbors])
      )
    }
    
    AA_meth[j, i] <- AA_meth[i, j]
  }
}

## CURRENT REFINEMENT ( W*ij = Wij (1 + lambdaAAij))

## Refinement equation for gene expression

# Normalize the Adamic-Adar matrices
AA_norm <- AA / max(AA)
AA_meth_norm <- AA_meth / max(AA_meth)


# Choose Lambda (later we can compare lambda = 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2)
lambda <- 1

# Refine gene expression affinity matrix

gene_refinement_factor <- 1 + lambda * AA_norm

gene_W_refined <- gene_W * gene_refinement_factor

# Refine methylation affinity matrix

meth_refinement_factor <- 1 + lambda * AA_meth_norm

meth_W_refined <- meth_W * meth_refinement_factor


## QUESTION B: DOES IT ACTUALLY IMPROVE SNF?

## Compute fused matrix with the refined affinities

W_fused_AA = SNF(
  list(gene_W_refined, meth_W_refined),
  K = 5,
  t = 30
)

dim(W_fused_AA)

# Eigengap analysis for AA-SNF

D_AA <- diag(rowSums(W_fused_AA))

L_AA <- D_AA - W_fused_AA

D_AA_inv_sqrt <- diag(
  1 / sqrt(diag(D_AA))
)

L_AA_norm <- D_AA_inv_sqrt %*%
  L_AA %*%
  D_AA_inv_sqrt

eig_AA <- eigen(
  L_AA_norm,
  symmetric = TRUE,
  only.values = TRUE
)

vals_AA <- sort(eig_AA$values)

vals_AA[1:20]

# Eigengaps
gaps_AA <- diff(vals_AA[1:20])

# Maximum eigengap
which.max(gaps_AA)
max(gaps_AA)

### Spectral clustering and ARI/NMI for for SNF/baseline and AA-SNF/baseline

### at Lambda = 1

# clusters, ARI/NMI

library(mclust)
library(aricode)

# ============================================================
# Spectral clustering: Baseline SNF vs AA-SNF (lambda = 1)
# K = 2 to 6
# ============================================================

results_lambda1 <- data.frame(
  K = 2:6,
  ARI_SNF = NA,
  NMI_SNF = NA,
  ARI_AA_SNF = NA,
  NMI_AA_SNF = NA
)

# ------------------------------------------------------------
# K = 2
# ------------------------------------------------------------

set.seed(123)
cluster_original_K2 <- spectralClustering(
  W_fused,
  K = 2
)

set.seed(123)
cluster_AA_K2 <- spectralClustering(
  W_fused_AA,
  K = 2
)

results_lambda1$ARI_SNF[1] <-
  adjustedRandIndex(cluster_original_K2, subtype_num)

results_lambda1$NMI_SNF[1] <-
  NMI(cluster_original_K2, subtype_num)

results_lambda1$ARI_AA_SNF[1] <-
  adjustedRandIndex(cluster_AA_K2, subtype_num)

results_lambda1$NMI_AA_SNF[1] <-
  NMI(cluster_AA_K2, subtype_num)


# ------------------------------------------------------------
# K = 3
# ------------------------------------------------------------

set.seed(123)
cluster_original_K3 <- spectralClustering(
  W_fused,
  K = 3
)

set.seed(123)
cluster_AA_K3 <- spectralClustering(
  W_fused_AA,
  K = 3
)

results_lambda1$ARI_SNF[2] <-
  adjustedRandIndex(cluster_original_K3, subtype_num)

results_lambda1$NMI_SNF[2] <-
  NMI(cluster_original_K3, subtype_num)

results_lambda1$ARI_AA_SNF[2] <-
  adjustedRandIndex(cluster_AA_K3, subtype_num)

results_lambda1$NMI_AA_SNF[2] <-
  NMI(cluster_AA_K3, subtype_num)


# ------------------------------------------------------------
# K = 4
# ------------------------------------------------------------

set.seed(123)
cluster_original_K4 <- spectralClustering(
  W_fused,
  K = 4
)

set.seed(123)
cluster_AA_K4 <- spectralClustering(
  W_fused_AA,
  K = 4
)

results_lambda1$ARI_SNF[3] <-
  adjustedRandIndex(cluster_original_K4, subtype_num)

results_lambda1$NMI_SNF[3] <-
  NMI(cluster_original_K4, subtype_num)

results_lambda1$ARI_AA_SNF[3] <-
  adjustedRandIndex(cluster_AA_K4, subtype_num)

results_lambda1$NMI_AA_SNF[3] <-
  NMI(cluster_AA_K4, subtype_num)


# ------------------------------------------------------------
# K = 5
# ------------------------------------------------------------

set.seed(123)
cluster_original_K5 <- spectralClustering(
  W_fused,
  K = 5
)

set.seed(123)
cluster_AA_K5 <- spectralClustering(
  W_fused_AA,
  K = 5
)

results_lambda1$ARI_SNF[4] <-
  adjustedRandIndex(cluster_original_K5, subtype_num)

results_lambda1$NMI_SNF[4] <-
  NMI(cluster_original_K5, subtype_num)

results_lambda1$ARI_AA_SNF[4] <-
  adjustedRandIndex(cluster_AA_K5, subtype_num)

results_lambda1$NMI_AA_SNF[4] <-
  NMI(cluster_AA_K5, subtype_num)


# ------------------------------------------------------------
# K = 6
# ------------------------------------------------------------

set.seed(123)
cluster_original_K6 <- spectralClustering(
  W_fused,
  K = 6
)

set.seed(123)
cluster_AA_K6 <- spectralClustering(
  W_fused_AA,
  K = 6
)

results_lambda1$ARI_SNF[5] <-
  adjustedRandIndex(cluster_original_K6, subtype_num)

results_lambda1$NMI_SNF[5] <-
  NMI(cluster_original_K6, subtype_num)

results_lambda1$ARI_AA_SNF[5] <-
  adjustedRandIndex(cluster_AA_K6, subtype_num)

results_lambda1$NMI_AA_SNF[5] <-
  NMI(cluster_AA_K6, subtype_num)


# ============================================================
# Final comparison
# ============================================================

results_lambda1

# ============================================================
# SNF vs AA-SNF comparison: K = 2 to 6
# Metrics:
#   1. Silhouette
#   2. Within-cluster similarity
#   3. Between-cluster similarity
#   4. Within/Between similarity ratio
#   5. Modularity
#   6. Normalized cluster entropy
# ============================================================

library(cluster)
library(igraph)


# ============================================================
# 1. Distance matrices
# ============================================================

D_original <- as.dist(1 - W_fused)

D_AA <- as.dist(1 - W_fused_AA)


# ============================================================
# 2. Cluster labels
# ============================================================

clusters_original <- list(
  K2 = snf_cluster_k2,
  K3 = snf_cluster_k3,
  K4 = snf_cluster_k4,
  K5 = snf_cluster_k5,
  K6 = snf_cluster_k6
)

clusters_AA <- list(
  K2 = cluster_AA_K2,
  K3 = cluster_AA_K3,
  K4 = cluster_AA_K4,
  K5 = cluster_AA_K5,
  K6 = cluster_AA_K6
)


# ============================================================
# 3. Within-cluster similarity
# ============================================================

within_similarity <- function(W, labels){
  
  clusts <- sort(unique(labels))
  vals <- c()
  
  for(cl in clusts){
    
    idx <- which(labels == cl)
    
    if(length(idx) > 1){
      
      subW <- W[idx, idx]
      
      vals <- c(
        vals,
        subW[upper.tri(subW)]
      )
    }
  }
  
  mean(vals)
}


# ============================================================
# 4. Between-cluster similarity
# ============================================================

between_similarity <- function(W, labels){
  
  vals <- c()
  clusts <- sort(unique(labels))
  
  for(i in 1:(length(clusts) - 1)){
    
    for(j in (i + 1):length(clusts)){
      
      idx1 <- which(labels == clusts[i])
      idx2 <- which(labels == clusts[j])
      
      vals <- c(
        vals,
        as.vector(W[idx1, idx2])
      )
    }
  }
  
  mean(vals)
}

# ============================================================
# 5. Normalized cluster entropy
# ============================================================

normalized_cluster_entropy <- function(labels){
  
  counts <- table(labels)
  proportions <- counts / sum(counts)
  
  entropy <- -sum(
    proportions * log(proportions)
  )
  
  normalized_entropy <- entropy / log(length(proportions))
  
  return(normalized_entropy)
}


# ============================================================
# 6. Create graph objects for modularity
# ============================================================

g_original <- graph_from_adjacency_matrix(
  W_fused,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

g_AA <- graph_from_adjacency_matrix(
  W_fused_AA,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)


# ============================================================
# 7. Calculate all metrics for K = 2 to 6
# ============================================================

results <- data.frame(
  K = 2:6,
  Silhouette_SNF = NA,
  Silhouette_AA_SNF = NA,
  Within_SNF = NA,
  Within_AA_SNF = NA,
  Between_SNF = NA,
  Between_AA_SNF = NA,
  Within_Between_SNF = NA,
  Within_Between_AA_SNF = NA,
  Modularity_SNF = NA,
  Modularity_AA_SNF = NA,
  Entropy_SNF = NA,
  Entropy_AA_SNF = NA
)


for(k in 2:6){
  
  # Get cluster labels
  labels_original <- clusters_original[[paste0("K", k)]]
  labels_AA <- clusters_AA[[paste0("K", k)]]
  
  
  # ----------------------------
  # Silhouette
  # ----------------------------
  
  sil_original <- silhouette(
    labels_original,
    D_original
  )
  
  sil_AA <- silhouette(
    labels_AA,
    D_AA
  )
  
  results$Silhouette_SNF[k - 1] <-
    mean(sil_original[, 3])
  
  results$Silhouette_AA_SNF[k - 1] <-
    mean(sil_AA[, 3])
  
  
  # ----------------------------
  # Within-cluster similarity
  # ----------------------------
  
  within_original <- within_similarity(
    W_fused,
    labels_original
  )
  
  within_AA <- within_similarity(
    W_fused_AA,
    labels_AA
  )
  
  results$Within_SNF[k - 1] <-
    within_original
  
  results$Within_AA_SNF[k - 1] <-
    within_AA
  
  
  # ----------------------------
  # Between-cluster similarity
  # ----------------------------
  
  between_original <- between_similarity(
    W_fused,
    labels_original
  )
  
  between_AA <- between_similarity(
    W_fused_AA,
    labels_AA
  )
  
  results$Between_SNF[k - 1] <-
    between_original
  
  results$Between_AA_SNF[k - 1] <-
    between_AA
  
  
  # ----------------------------
  # Within / Between ratio
  # ----------------------------
  
  results$Within_Between_SNF[k - 1] <-
    within_original / between_original
  
  results$Within_Between_AA_SNF[k - 1] <-
    within_AA / between_AA
  
  
  # ----------------------------
  # Modularity
  # ----------------------------
  
  results$Modularity_SNF[k - 1] <-
    modularity(
      g_original,
      membership = labels_original,
      weights = E(g_original)$weight
    )
  
  results$Modularity_AA_SNF[k - 1] <-
    modularity(
      g_AA,
      membership = labels_AA,
      weights = E(g_AA)$weight
    )
  
  
  # ----------------------------
  # Normalized cluster entropy
  # ----------------------------
  
  results$Entropy_SNF[k - 1] <-
    normalized_cluster_entropy(
      labels_original
    )
  
  results$Entropy_AA_SNF[k - 1] <-
    normalized_cluster_entropy(
      labels_AA
    )
}


# ============================================================
# Final results
# ============================================================

results

# AA-SNF contingency tables
# Comparing AA-SNF clusters with known breast-cancer subtypes

table(
  AA_SNF_cluster = cluster_AA_K2,
  Subtype = subtype
)

table(
  AA_SNF_cluster = cluster_AA_K3,
  Subtype = subtype
)

table(
  AA_SNF_cluster = cluster_AA_K4,
  Subtype = subtype
)

table(
  AA_SNF_cluster = cluster_AA_K5,
  Subtype = subtype
)

table(
  AA_SNF_cluster = cluster_AA_K6,
  Subtype = subtype
)

## focus on K =5 and K = 6

AA_table_K5 <- table(
  AA_SNF_cluster = cluster_AA_K5,
  Subtype = subtype
)

AA_table_K6 <- table(
  AA_SNF_cluster = cluster_AA_K6,
  Subtype = subtype
)

AA_table_K5
AA_table_K6

## Sensitivity plot of lambda

##install.packages("readODS")
library(readODS)

results <- read_ods("Adamic-Adar-ARI-NMI-Breast Cancer.ods")
nrow(results)
names(results)
unique(results$lambda)
unique(results$K)
str(results)


## Plot 1: NMI versus lambda
library(ggplot2)

ggplot(results, aes(x = lambda, y = NMI, group = K, color = factor(K))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Effect of Adamic–Adar Refinement on NMI",
    x = expression(lambda),
    y = "NMI",
    color = "Number of clusters (K)"
  ) +
  theme_minimal()

# Plot 2: ARI versus lambda

ggplot(results, aes(x = lambda, y = ARI, group = K, color = factor(K))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Effect of Adamic–Adar Refinement on ARI",
    x = expression(lambda),
    y = "ARI",
    color = "Number of clusters (K)"
  ) +
  theme_minimal()

## Plot 3: NMI and ARI versus lambda for gold-standard K = 6

# Create K = 6 subset:

results_K6 <- subset(results, K == 6)
results_K6

results_K6 <- subset(results, K == 6)

ggplot(results_K6, aes(x = lambda)) +
  geom_line(aes(y = NMI, color = "NMI", linetype = "NMI"), linewidth = 1) +
  geom_point(aes(y = NMI, color = "NMI", linetype = "NMI"), size = 2) +
  geom_line(aes(y = ARI, color = "ARI", linetype = "ARI"), linewidth = 1) +
  geom_point(aes(y = ARI, color = "ARI", linetype = "ARI"), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Effect of λ on Recovery of Six Breast-Cancer Subtypes",
    x = expression(lambda),
    y = "Clustering agreement",
    color = "Metric",
    linetype = "Metric"
  ) +
  theme_minimal()


## Histogram plot at lambda = 1.
# Make sure to check above that you have been working with lambda = 1 throughout.

## we are using only gene.No need to plot histogram for methylation; 
## we just need one omic plot to show what adamic adar refinement is doing)

# First verify that the refinement actually changed the network

# Compare refined and original gene-expression affinity matrices

sum(gene_W_refined > gene_W)
sum(gene_W_refined == gene_W)
sum(gene_W_refined < gene_W)

# Difference between refined and original affinity
difference <- gene_W_refined - gene_W

# check
range(difference)

# Keep only unique patient pairs

difference_pairs <- difference[upper.tri(difference)]
length(difference_pairs)

## Histogram 1
# Distribution of W* - W for all patient pairs

hist(
  difference_pairs,
  breaks = 50,
  main = "Change in Affinity After Adamic–Adar Refinement",
  xlab = "W^* - W",
  ylab = "Number of Patient Pairs"
)

## Histogram 2
# removes unchanged (zero) relationships and looks only at the relationships
# that actually increased

# Keep only patient pairs whose affinity increased
positive_changes <- difference_pairs[difference_pairs > 0]
length(positive_changes)

log_changes <- log10(positive_changes)

hist(
  log_changes,
  breaks = 50,
  main = "Distribution of Affinity Increases After Adamic–Adar Refinement",
  xlab = "log[10](W* - W)",
  ylab = "Number of Patient Pairs"
)
## Create the refinement factor 
# Refinement factor (this is already created way up above with lambda = 1)
# but just recreate here. it's ok as long as you're still working with lambda = 1.

refinement_factor <- 1 + lambda * AA_norm

# Keep only unique patient pairs
refinement_factor_pairs <- refinement_factor[upper.tri(refinement_factor)]
range(refinement_factor_pairs)

## Histogram 3
# Distribution of refinement factors

hist(
  refinement_factor_pairs,
  breaks = 50,
  main = "Distribution of Topology-Based Refinement Factors",
  xlab = expression(1 + lambda * tilde(AA)),
  ylab = "Number of Patient Pairs"
)

## Histogram 4 : the important subtype analysis
# Does Adamic-Adar preferentially strengthen relationships between patients
# belonging to the same molecular subtypes?

# Make sure subtype labels are aligned with the affinity matrix
identical(rownames(W_fused), cancer_subtypes$TAX)

# creat this again. Alrady done above, but just redo it
subtype <- cancer_subtypes$hu_subtype

length(subtype)
table(subtype)

# Get the patient-pair differences.
# We already have 'difference_pairs', but now we need to know which patient pair
# produced each differerence so that we can determine whether the two patients 
# have the same subtype.

# Create upper-triangle indices:

# Identify the patient pairs represented in difference_pairs
# Identify the patient pairs represented in difference_pairs
pair_indices <- which(
  upper.tri(difference),
  arr.ind = TRUE
)
head(pair_indices)

# Determine whether each pair has the same subtype

# Determine whether each patient pair has the same molecular subtype
same_subtype <- subtype[pair_indices[, 1]] ==
  subtype[pair_indices[, 2]]
table(same_subtype)

# Put everything into one data frame

# Create a data frame for the subtype comparison
histogram_data <- data.frame(
  difference = difference_pairs,
  same_subtype = same_subtype
)
head(histogram_data)

# The histogram plot - Affinity change: same subtype vs different subtype
# this plot will give two distributions overlaid:
# Same subtype: pairs where both patients belong to the same known molecular subtype.
# Different subtype: pairs where the patients belong to different known molecular subtypes.

# Keep only patient pairs whose affinity actually increased
histogram_data_positive <- histogram_data[
  histogram_data$difference > 0,
]

# Log10 transform the affinity increases
histogram_data_positive$log_difference <- 
  log10(histogram_data_positive$difference)

# Plot
ggplot(
  histogram_data_positive,
  aes(
    x = log_difference,
    fill = same_subtype
  )
) +
  geom_histogram(
    bins = 50,
    alpha = 0.6,
    position = "identity"
  ) +
  labs(
    title = "Affinity Increases by Breast-Cancer Subtype Relationship",
    x = "log[10](W* - W)",
    y = "Number of Patient Pairs",
    fill = "Same subtype"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "steelblue",
      "FALSE" = "orange"
    ),
    labels = c(
      "TRUE" = "Same subtype",
      "FALSE" = "Different subtype"
    )
  ) +
  theme_minimal()

# 2nd version of code for Histogram 4 - keep the x axis positive (preferred)
# since the change in affinites are all greater than zero.

# Keep only patient pairs whose affinity increased
histogram_data_positive <- histogram_data[
  histogram_data$difference > 0,
]

library(ggplot2)
library(scales)

ggplot(
  histogram_data_positive,
  aes(
    x = difference,
    fill = same_subtype
  )
) +
  geom_histogram(
    bins = 50,
    alpha = 0.6,
    position = "identity"
  ) +
  scale_x_log10(
    labels = label_scientific(digits = 1)
  ) +
  labs(
    title = "Affinity Increases by Breast-Cancer Subtype Relationship",
    x = "Increase in affinity (W* - W)",
    y = "Number of Patient Pairs",
    fill = "Same subtype"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "steelblue",
      "FALSE" = "orange"
    ),
    labels = c(
      "TRUE" = "Same subtype",
      "FALSE" = "Different subtype"
    )
  ) +
  theme_minimal()

## To answer the question: of all same-subtype pairs, what percentage were modified?
# and of all different-subtype pairs, what percentage were modified?

# Step 1: First, look at how many pairs you have in each category:
table(histogram_data$same_subtype)

# Step 2: Count modified vs. unmodified pairs
table(
  histogram_data$same_subtype,
  histogram_data$difference > 0
)
# Step 3: Calculate the proportions 
prop.table(
  table(
    histogram_data$same_subtype,
    histogram_data$difference > 0
  ),
  margin = 1
)
# Step 4: Covert the proportions to percentages
proportions <- prop.table(
  table(
    histogram_data$same_subtype,
    histogram_data$difference > 0
  ),
  margin = 1
)

proportions * 100

## making a clean table
modification_summary <- data.frame(
  Relationship = c(
    "Different subtype",
    "Same subtype"
  ),
  Total_pairs = c(
    12895,
    3036
  ),
  Modified_pairs = c(
    2013,
    1353
  ),
  Percent_modified = c(
    2013 / 12895 * 100,
    1353 / 3036 * 100
  )
)

modification_summary   # This is for gene expression

## GENE EXPRESSION: Bar plot of the percentage modification. 
#Put the actual percentages on top of the bars

library(ggplot2)

ggplot(
  modification_summary,
  aes(
    x = Relationship,
    y = Percent_modified,
    fill = Relationship
  )
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percent_modified)),
    vjust = -0.5,
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Proportion of Patient Pairs with Increased Gene Affinity",
    x = "Patient-pair relationship",
    y = "Pairs with increased affinity (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

## Methylation: Affinity Percentage increase

# Calculate how much each methylation affinity changed
meth_difference <- meth_W_refined - meth_W

# Check the dimensions
dim(meth_difference)

# Inspect the range of changes
range(meth_difference)

# Look at the first few rows and columns
meth_difference[1:5, 1:5]

# Extract only the upper triangle
# so that each patient pair is counted once
meth_difference_pairs <- meth_difference[
  upper.tri(meth_difference)
]

# Check number of unique patient pairs
length(meth_difference_pairs)

# Number of pairs with increased methylation affinity
sum(meth_difference_pairs > 0)

# Number of pairs with no increase
sum(meth_difference_pairs == 0)

# Check whether there are any negative changes
sum(meth_difference_pairs < 0)

# Create the methyaltion patient-pair data
cancer_subtypes$hu_subtype   # we already have this. No need to run it. 
# And the patient order was previously verified against the SNF matrices.

# Get the patient indices for each unique pair
pair_indices_meth <- which(
  upper.tri(meth_difference),
  arr.ind = TRUE
)

# Create the methylation pair-level data frame
meth_histogram_data <- data.frame(
  patient1 = rownames(meth_W)[pair_indices_meth[, 1]],
  patient2 = rownames(meth_W)[pair_indices_meth[, 2]],
  difference = meth_difference_pairs
)

# Add subtype information
meth_histogram_data$subtype1 <- cancer_subtypes$hu_subtype[
  pair_indices_meth[, 1]
]

meth_histogram_data$subtype2 <- cancer_subtypes$hu_subtype[
  pair_indices_meth[, 2]
]

# Determine whether the two patients have the same subtype
meth_histogram_data$same_subtype <-
  meth_histogram_data$subtype1 ==
  meth_histogram_data$subtype2

# Check the result
head(meth_histogram_data)
# We now have one row for every unique pair

# Check the number of same-subtype and different-subtype pairs
table(meth_histogram_data$same_subtype)

# Determine how many pairs were increased
table(
  meth_histogram_data$same_subtype,
  meth_histogram_data$difference > 0
)

# Calculate the percentage increased
prop.table(
  table(
    meth_histogram_data$same_subtype,
    meth_histogram_data$difference > 0
  ),
  margin = 1
)
# Create the methylation summary table
meth_modification_summary <- data.frame(
  Relationship = c(
    "Different subtype",
    "Same subtype"
  ),
  Total_pairs = c(
    sum(meth_histogram_data$same_subtype == FALSE),
    sum(meth_histogram_data$same_subtype == TRUE)
  ),
  Modified_pairs = c(
    sum(
      meth_histogram_data$same_subtype == FALSE &
        meth_histogram_data$difference > 0
    ),
    sum(
      meth_histogram_data$same_subtype == TRUE &
        meth_histogram_data$difference > 0
    )
  )
)

# Calculate percentage modified
meth_modification_summary$Percent_modified <-
  meth_modification_summary$Modified_pairs /
  meth_modification_summary$Total_pairs * 100

# Display the table
meth_modification_summary

# Methylation bar plot
library(ggplot2)

ggplot(
  meth_modification_summary,
  aes(
    x = Relationship,
    y = Percent_modified,
    fill = Relationship
  )
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percent_modified)),
    vjust = -0.5,
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Proportion of Patient Pairs with Increased Methylation Affinity",
    x = "Patient-pair relationship",
    y = "Pairs with increased affinity (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

