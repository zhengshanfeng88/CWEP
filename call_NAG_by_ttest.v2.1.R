#!/usr/bin/env Rscript
# author: zhengshanfeng and biotang, 2019.09
# version 2.1

#########################
### command arguments ###
#########################
opt <- commandArgs(TRUE)
# input
file_sample_list <- opt[1]
file_expr_data   <- opt[2]
# output
file_ttest_out   <- opt[3]

# example
# input
#file_sample_list <- ""
#file_expr_data   <- ""
# output
#file_ttest_out   <- ""

#################
### functions ###
#################
get_group_avg <- function(expr_data, group_idx_list){
  # calculate average of samples in one group
  avg <- data.frame()
  for(group in names(group_idx_list) ){
    idx <- group_idx_list[[group]]
    cal_group_avg <- function(x){
      y <- NA
      if(length(x) == 1) y <- x
      if(length(!is.na(x)) >= 1) y <- mean(x, na.rm = T)
      return(y)
    }
    value <- apply(expr_data[,idx], 1, cal_group_avg)
    #
    if(length(avg) >1 )  avg <- cbind(avg, value)
    if(length(avg) == 0) avg <- value
  }
  colnames(avg) <- names(group_idx_list)
  rownames(avg) <- rownames(expr_data)
  avg <- as.matrix(avg)
  return(avg)
}

#
get_group_sd <- function(expr_data, group_idx_list){
  # calculate standard deviation of samples in one group
  sd <- data.frame()
  for(group in names(group_idx_list) ){
    idx <- group_idx_list[[group]]
    cal_group_sd <- function(x){
      y <- NA
      if(length(x) == 1) y <- NA
      if(length(!is.na(x)) >= 1) y <- sd(x, na.rm = T)
      return(y)
    }
    value <- apply(expr_data[,idx], 1, cal_group_sd)
    #
    if(length(sd) >1 )  sd <- cbind(sd, value)
    if(length(sd) == 0) sd <- value
  }
  colnames(sd) <- names(group_idx_list)
  rownames(sd) <- rownames(sd)
  sd <- as.matrix(sd)
  return(sd)
}

#
filter_low_expr_of_group <- function(x, group_idx_list){
  # calculate average of samples in one group
  PASS <- FALSE
  for(group in names(group_idx_list) ){
    idx <- group_idx_list[[group]]
    if(mean(x[idx], na.rm = T) > 1 & length(which(x[idx] > 0) ) >= 2) PASS <- TRUE
  }
  return(PASS)
}
#######################
### read data       ###
#######################
# information of samples and experiment design
cat("reading information of experiment design\n")
sample_list <- read.table(file_sample_list, header = TRUE, stringsAsFactors = FALSE)
print(sample_list)

# read trait value
cat("reading data of expression or trait value\n")
expr_data <- read.table(file_expr_data, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

#######################
### check           ###
#######################
#
if(!"ids" %in% colnames(sample_list) ) stop("The column 'ids' of sample list is missing\n")
#if(!"reportName" %in% colnames(sample_list) ) stop("The column 'reportName' of sample list is missing\n")
if(!"group" %in% colnames(sample_list) ) stop("The column 'group' of sample list is missing\n")

#
if(any(!sample_list$ids %in% colnames(expr_data) ) ) stop("Sample is not existed in this data\n")

# the first column need to be unique
if(any(duplicated(expr_data[[1]]) ) ) stop("The first column need to be unique\n")
rownames(expr_data) <- expr_data[[1]]
head(expr_data)

#
if(!"pedigree" %in% colnames(sample_list) ) stop("The column 'pedigree' of sample list is missing\n")
if(any(!c("P", "H", "M") %in% sample_list$pedigree) ) stop("Pedigree information is missing\n")

#######################
### process data    ###
#######################
# select samplaes and change column names
cat("just select samples in experiment design table\n")
expr_data <- expr_data[,sample_list$ids]
# change column name if reportName existed
if("reportName" %in% colnames(sample_list) ) colnames(expr_data) <- sample_list$reportName
head(expr_data)

############################
### process sample list  ###
############################
sample_list$sample <- sample_list$ids
if("reportName" %in% colnames(sample_list) ) sample_list$sample <- sample_list$reportName
print(sample_list)

#############################################
### sample index and names of each group  ###
#############################################
group_name <- unique(sample_list$group)
group_idx_list <- list()
group_idx_name <- list()
for(group in group_name){
  idx <- which(sample_list$group == group)
  if(length(idx) > 0) group_idx_list[[group]] <- idx
  if(length(idx) > 0) group_idx_name[[group]] <- sample_list$sample[idx]
}

#
idx.P <- which(sample_list$pedigree == "P")
idx.M <- which(sample_list$pedigree == "M")
idx.H <- which(sample_list$pedigree == "H")

#
group_idx_list <- list(P = idx.P, M = idx.M, H = idx.H)
group_idx_name <- list(
  P = sample_list$sample[idx.P], 
  M = sample_list$sample[idx.M], 
  H = sample_list$sample[idx.H]
)

#############################################
### filter low expression of each group   ###
#############################################
PASS <- apply(expr_data, 1, function(x) filter_low_expr_of_group(x, group_idx_name) )
table(PASS)

################################
### sample size of groups    ###
################################
n.P <- apply(expr_data, 1, function(x) length(which(!is.na(x[idx.P]) ) ) )
n.M <- apply(expr_data, 1, function(x) length(which(!is.na(x[idx.M]) ) ) )
n.H <- apply(expr_data, 1, function(x) length(which(!is.na(x[idx.H]) ) ) )

################################
### calculate mean and sd    ###
################################
# calculate mean of group
cat("calculate mean of group (pedigree)\n")
mtx.avg <- get_group_avg(as.matrix(expr_data), group_idx_name)
head(mtx.avg)

# calculate sd of group
cat("calculate sd of group (pedigree)\n")
mtx.sd <- get_group_sd(as.matrix(expr_data), group_idx_name)
#cat("change colname\n")
colnames(mtx.sd) <- paste("sd.", names(group_idx_name), sep = "")
head(mtx.sd)

################################
### t-test for two samples   ###
################################
#
ttest_xy <- function(x, y){
  p <- NA
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if(length(x) >= 2 | length(y) >= 2){
    p <- t.test(x, y)$p.value
    if(mean(x) == mean(y) & sd(x) == sd(y)) p <- 1
  }
  return(p)
}

#
pvalue.PvsM <- apply(expr_data, 1, function(x) ttest_xy(x[idx.P], x[idx.M]))
pvalue.PvsH <- apply(expr_data, 1, function(x) ttest_xy(x[idx.P], x[idx.H]))
pvalue.MvsH <- apply(expr_data, 1, function(x) ttest_xy(x[idx.M], x[idx.H]))

#
pvalue.PvsM <- as.matrix(pvalue.PvsM)
pvalue.PvsH <- as.matrix(pvalue.PvsH)
pvalue.MvsH <- as.matrix(pvalue.MvsH)

#
colnames(pvalue.PvsM) <- c("pvalue.PvsM")
colnames(pvalue.PvsH) <- c("pvalue.PvsH")
colnames(pvalue.MvsH) <- c("pvalue.MvsH")

#
head(pvalue.PvsM)
head(pvalue.PvsH)
head(pvalue.MvsH)

################################
### Middle Parent value      ###
################################
#
cat("calculate MPV\n")
MPV <- ifelse((!is.na(mtx.avg[,c("M")] & !is.na(mtx.avg[,c("P")]) ) ), rowMeans(mtx.avg[,c("M", "P")]), NA)
MPV <- as.matrix(MPV)
colnames(MPV) <- c("MPV")
head(MPV)

#
estimate_MPV_var <- function(var.P, n.P, var.M, n.M) 0.25*var.P/n.P + 0.25*var.M/n.M
var.MPV <- ifelse(!is.na(MPV), estimate_MPV_var(mtx.sd[,c("sd.P")]^2, n.P, mtx.sd[,c("sd.M")]^2, n.M), NA)
sd.MPV  <- as.matrix(sqrt(var.MPV) )
colnames(sd.MPV) <- c("sd.MPV")
head(sd.MPV)

################################
### t-test for H vs MPV      ###
################################
#
calculate_t_statistic <- function(MPV, HV, var.MPV, var.HV, n.MPV, n.H){
  t_statistic <- abs(MPV - HV) / sqrt( (var.MPV/n.MPV + var.HV/n.H) )
}
#
estimate_df_of_t <- function(var.MPV, var.HV, n.MPV, n.H){
  k  <- (var.MPV/n.MPV) / ((var.MPV/n.MPV) + (var.HV/n.H))
  df.MPV <- n.MPV - 1
  df.H   <- n.H - 1
  df <- 1 / (k^2/df.MPV + (1-k)^2/df.H) 
}
#
cal_pvalue_of_ttest <- function(MPV, HV, var.MPV, var.HV, n.MPV, n.H){
  t_statistic <- calculate_t_statistic(MPV, HV, var.MPV, var.HV, n.MPV, n.H)
  df <- estimate_df_of_t(var.MPV, var.HV, n.MPV, n.H)
  p_value <- 2*pt(q = t_statistic, df = df, lower.tail = FALSE)
}

#
pvalue.HvsMPV <- ifelse(
  (
   !is.na(MPV) & !is.na(mtx.avg[,c("H")]) & 
   !is.na(sd.MPV) & !is.na(mtx.sd[,c("sd.H")]) & 
   sd.MPV+mtx.sd[,c("sd.H")] != 0 & 
   (n.P+n.M)/2 >= 2 & n.H >= 2 &
   !is.na(pvalue.MvsH) & !is.na(pvalue.PvsH)
  ),
  cal_pvalue_of_ttest(MPV, mtx.avg[,c("H")], sd.MPV^2, mtx.sd[,c("sd.H")]^2, (n.P+n.M)/2, n.H),
  NA
)
#pvalue.HvsMPV <- cal_pvalue_of_ttest(MPV, mtx.avg[,c("H")], sd.MPV^2, mtx.sd[,c("sd.H")]^2, (n.P+n.M)/2, n.H)
pvalue.HvsMPV <- as.matrix(pvalue.HvsMPV)
colnames(pvalue.HvsMPV) <- c("pvalue.HvsMPV")
head(pvalue.HvsMPV)

################################
### statistics of MPH        ###
################################
#
IndexMPH <- ifelse((!is.na(MPV) & !is.na(mtx.avg[,c("H")]) ), (mtx.avg[,c("H")] - MPV) / (MPV + 0.01), NA)
IndexMPH <- ifelse(is.infinite(IndexMPH), NA, IndexMPH)
IndexMPH <- as.matrix(IndexMPH)
colnames(IndexMPH) <- c("IndexMPH")
head(IndexMPH)

#
AdditiveEffect <- ifelse((!is.na(mtx.avg[,c("P")]) & !is.na(mtx.avg[,c("M")])), abs(mtx.avg[,c("P")] - mtx.avg[,c("M")])/2, NA)
AdditiveEffect <- as.matrix(AdditiveEffect)
colnames(AdditiveEffect) <- c("AdditiveEffect")
head(AdditiveEffect)

#
NonadditiveEffect <- ifelse((!is.na(mtx.avg[,c("H")]) & !is.na(MPV) ), mtx.avg[,c("H")] - MPV, NA)
NonadditiveEffect <- as.matrix(NonadditiveEffect)
colnames(NonadditiveEffect) <- c("NonadditiveEffect")
head(NonadditiveEffect)

#
DominanceDegree <- ifelse((!is.na(AdditiveEffect) & !is.na(NonadditiveEffect) ), abs(NonadditiveEffect / (AdditiveEffect + 0.01) ), NA)
DominanceDegree <- ifelse(is.infinite(DominanceDegree), NA, DominanceDegree)
DominanceDegree <- as.matrix(DominanceDegree)
colnames(DominanceDegree) <- c("DominanceDegree")
head(DominanceDegree)

################################
### fold change              ###
################################
#
logFC.HvsMPV <- ifelse((!is.na(MPV) & !is.na(mtx.avg[,c("H")]) ), log2( (mtx.avg[,c("H")] + 0.01) / (MPV + 0.01) ), NA)
logFC.HvsMPV <- ifelse(is.infinite(logFC.HvsMPV), NA, logFC.HvsMPV)
logFC.HvsMPV <- as.matrix(logFC.HvsMPV)
colnames(logFC.HvsMPV) <- c("logFC.HvsMPV")
head(logFC.HvsMPV)

#
logFC.PvsM <- ifelse((!is.na(mtx.avg[,c("P")]) & !is.na(mtx.avg[,c("M")]) ), log2( (mtx.avg[,c("P")] + 0.01) / (mtx.avg[,c("M")] + 0.01) ), NA)
logFC.PvsM <- ifelse(is.infinite(logFC.PvsM), NA, logFC.PvsM)
logFC.PvsM <- as.matrix(logFC.PvsM)
colnames(logFC.PvsM) <- c("logFC.PvsM")
head(logFC.PvsM)

#
logFC.PvsH <- ifelse((!is.na(mtx.avg[,c("P")]) & !is.na(mtx.avg[,c("H")]) ), log2( (mtx.avg[,c("P")] + 0.01) / (mtx.avg[,c("H")] + 0.01) ), NA)
logFC.PvsH <- ifelse(is.infinite(logFC.PvsH), NA, logFC.PvsH)
logFC.PvsH <- as.matrix(logFC.PvsH)
colnames(logFC.PvsH) <- c("logFC.PvsH")
head(logFC.PvsH)

#
logFC.MvsH <- ifelse((!is.na(mtx.avg[,c("M")]) & !is.na(mtx.avg[,c("H")]) ), log2( (mtx.avg[,c("M")] + 0.01) / (mtx.avg[,c("H")] + 0.01) ), NA)
logFC.MvsH <- ifelse(is.infinite(logFC.MvsH), NA, logFC.MvsH)
logFC.MvsH <- as.matrix(logFC.MvsH)
colnames(logFC.MvsH) <- c("logFC.MvsH")
head(logFC.MvsH)

################################
### write result             ###
################################
# create dir of output
if(!dir.exists(dirname(file_ttest_out) ) ) dir.create(dirname(file_ttest_out), recursive = T)

#
result <- data.frame(
  id = rownames(expr_data), 
  expr_data, 
  PASS,
  mtx.avg, 
  MPV, 
  mtx.sd, 
  sd.MPV, 
  pvalue.PvsM, 
  pvalue.PvsH, 
  pvalue.MvsH, 
  pvalue.HvsMPV, 
  logFC.PvsM,
  logFC.PvsH,
  logFC.MvsH,
  logFC.HvsMPV,
  IndexMPH, 
  AdditiveEffect,
  NonadditiveEffect,
  DominanceDegree,
  check.names = F
)
#
write.table(result, file_ttest_out, quote = F, sep = "\t", row.names = F)

