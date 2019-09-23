#!/usr/bin/env Rscript
# author: zhengshanfeng and biotang, 2019.09
# version 2.5

# load packages
library(ggplot2)
library(gridExtra)

#########################
### command arguments ###
#########################
opt <- commandArgs(TRUE)
# input
file_pedigree_conf <- opt[1]
file_pedigree_data <- opt[2]
# output
output_prefix      <- opt[3]

# example
# input
#file_pedigree_conf <- "example/Cross_ZS97AxMH63.conf.xls"
#file_pedigree_data <- "example/Cross_ZS97AxMH63.DEG_expr.xls"
# output
#output_prefix      <- "example/result/Cross_ZS97AxMH63"

# create output dir
if(!dir.exists(dirname(output_prefix) ) ) dir.create(dirname(output_prefix), recursive = T)

###########################
### global arguments    ###
###########################


###########################
### structure of output ###
###########################
file_expr_pattern  <- paste(output_prefix, ".hpep.xls", sep = "")
dir_summary        <- paste(dirname(output_prefix), "/summary/", sep = "") 

# output files under summary
# create dir
if(!dir.exists(dir_summary) ) dir.create(dir_summary)
# 
file_class_count   <- paste(dir_summary, basename(output_prefix), ".hpep_class_count.xls", sep = "")
file_diclass_count <- paste(dir_summary, basename(output_prefix), ".hpep_diClass_count.xls", sep = "")
file_super_count   <- paste(dir_summary, basename(output_prefix), ".hpep_superClass_count.xls", sep = "")
#
tiff_count_barplot <- paste(dir_summary, basename(output_prefix), ".hpep_count_barplot.tiff", sep = "")
tiff_rate_barplot  <- paste(dir_summary, basename(output_prefix), ".hpep_rate_barplot.tiff", sep = "")

###########################
### list of classCode   ###
###########################
classCode_list <- data.frame(
  ClassCode = c(
    "ClassNo.01",
    "ClassNo.02",
    "ClassNo.03",
    "ClassNo.04",
    "ClassNo.05",
    "ClassNo.06",
    "ClassNo.07",
    "ClassNo.08",
    "ClassNo.09",
    "ClassNo.10",
    "ClassNo.11",
    "ClassNo.12",
    "ClassNo.13",
    "ClassNo.14",
    "ClassNo.15",
    "ClassNo.16",
    "ClassNo.17",
    "ClassNo.18",
    "ClassNo.19",
    "ClassNo.20",
    "ClassNo.21",
    "ClassNo.22",
    "ClassNo.23",
    "ClassNo.24"
  ),  
  ExprPattern = c(
    # AHP
    "H>P>M",  # No.01
    "H>M>P",  # No.02
    "H>P>=M", # No.03
    "H>M>=P", # No.04
    # HP
    "H>=P>=M",# No.05
    "H>=M>=P",# No.06
    "H>=P>M", # No.07
    "P>=H>M", # No.08
    "H>=M>P", # No.09
    "M>=H>P", # No.10
    # MP
    "P>=H>=M",# No.11
    "P>H>M",  # No.12
    "M>H>P",  # No.13
    "M>=H>=P",# No.14
    # LP
    "M>H>=P", # No.15
    "M>P>=H", # No.16
    "P>H>=M", # No.17
    "P>M>=H", # No.18
    "M>=P>=H",# No.19
    "P>=M>=H",# No.20
    # BLP
    "M>=P>H", # No.21
    "P>=M>H", # No.22
    "M>P>H",  # No.23
    "P>M>H"   # No.24
  ),
  DiclassCode = c(
    "DiclassNo.0102", "DiclassNo.0102",
    "DiclassNo.0304", "DiclassNo.0304",
    "DiclassNo.0506", "DiclassNo.0506",
    "DiclassNo.0709", "DiclassNo.0810",
    "DiclassNo.0709", "DiclassNo.0810",
    "DiclassNo.1114", "DiclassNo.1213",
    "DiclassNo.1213", "DiclassNo.1114",
    "DiclassNo.1517", "DiclassNo.1618",
    "DiclassNo.1517", "DiclassNo.1618",
    "DiclassNo.1920", "DiclassNo.1920",
    "DiclassNo.2122", "DiclassNo.2122",
    "DiclassNo.2324", "DiclassNo.2324"
  ),
  DiclassPattern = c(
    "H>HP>LP",   "H>HP>LP",
    "H>HP>=LP",  "H>HP>=LP",
    "H>=HP>=LP", "H>=HP>=LP",
    "H>=HP>LP",  "HP>=H>LP",
    "H>=HP>LP",  "HP>=H>LP",
    "HP>=H>=LP", "HP>H>LP",
    "HP>H>LP",   "HP>=H>=LP",
    "HP>H>=LP",  "HP>LP>=H",
    "HP>H>=LP",  "HP>LP>=H",
    "HP>=LP>=H", "HP>=LP>=H",
    "HP>=LP>H",  "HP>=LP>H",
    "HP>LP>H",   "HP>LP>H"
  ),
  Superclass = c(
    "AHP", "AHP", "AHP", "AHP",
    "HP", "HP", "HP", "HP", "HP", "HP",
    "MP", "MP", "MP", "MP",
    "LP", "LP", "LP", "LP", "LP", "LP",
    "BLP", "BLP", "BLP", "BLP"
  ),
  stringsAsFactors = F,
  check.names = F
)
rownames(classCode_list)      <- classCode_list$ClassCode
classCode_list$ClassCode      <- factor(classCode_list$ClassCode, levels = as.character(classCode_list$ClassCode) )
classCode_list$ExprPattern    <- factor(classCode_list$ExprPattern, levels = as.character(classCode_list$ExprPattern) )
classCode_list$DiclassCode    <- factor(classCode_list$DiclassCode, levels = as.character(unique(classCode_list$DiclassCode) ) )
classCode_list$DiclassPattern <- factor(classCode_list$DiclassPattern, levels = as.character(unique(classCode_list$DiclassPattern) ) )
classCode_list$Superclass     <- factor(classCode_list$Superclass, levels = c("AHP", "HP", "MP", "LP", "BLP"))
classCode_list

#
diClass_list <- data.frame(
  DiclassCode = c(
    "DiclassNo.0102",
    "DiclassNo.0304",
    "DiclassNo.0506",
    "DiclassNo.0709",
    "DiclassNo.0810",
    "DiclassNo.1114",
    "DiclassNo.1213",
    "DiclassNo.1517",
    "DiclassNo.1618",
    "DiclassNo.1920",
    "DiclassNo.2122",
    "DiclassNo.2324"
  ),
  DiclassPattern = c(
    "H>HP>LP",
    "H>HP>=LP",
    "H>=HP>=LP",
    "H>=HP>LP", 
    "HP>=H>LP",
    "HP>=H>=LP",
    "HP>H>LP",
    "HP>H>=LP",
    "HP>LP>=H",
    "HP>=LP>=H",
    "HP>=LP>H",
    "HP>LP>H"
  ),
  Superclass = c(
    "AHP", "AHP",
    "HP", "HP", "HP",
    "MP", "MP",
    "LP", "LP", "LP",
    "BLP", "BLP"
  ),
  stringsAsFactors = F,
  check.names = F  
)
rownames(diClass_list)      <- diClass_list$DiclassCode
diClass_list$DiclassCode    <- factor(diClass_list$DiclassCode, levels = as.character(diClass_list$DiclassCode) )
diClass_list$DiclassPattern <- factor(diClass_list$DiclassPattern, levels = as.character(diClass_list$DiclassPattern) )
diClass_list$Superclass     <- factor(diClass_list$Superclass, levels = c("AHP", "HP", "MP", "LP", "BLP"))
diClass_list

#
superClass_list <- data.frame(
  SuperClass = c("AHP", "HP", "MP", "LP", "BLP")
)
superClass_list$SuperClass <- factor(superClass_list$SuperClass, levels = superClass_list$SuperClass)
superClass_list

#################
### functions ###
#################
get_avg <- function(expr, group_name, sample_list){
  # calculate average of samples in one group
  sample_name <- colnames(expr)
  avg <- data.frame()
  for(group in group_name){
    idx <- which(sample_list$group == group)
    samples <- sample_list$sample[idx]
    #jdx <- which(sample_name %in% sample_list$ids[idx])
    if(length(samples) == 1) value <- expr[,samples]
    if(length(samples) > 1 ) value <- rowMeans(expr[,samples])
    #
    if(length(avg) >1) avg <- cbind(avg, value)
    if(length(avg) == 0) avg <- value
  }
  colnames(avg) <- group_name
  rownames(avg) <- rownames(expr)
  avg <- as.matrix(avg)
  return(avg)
}

# class-mode
identify_express_pattern <- function(x){
  classCode <- NA
  # Criteria 1
  # AHP
  if(x["P"]> x["M"] && x["H"]> x["P"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.01" # H>P>M
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.02" # H>M>P
  if(x["P"]>=x["M"] && x["H"]> x["P"] && x["H"]> x["M"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.03" # H>P>=M
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["H"]> x["M"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.04" # H>M>=P(M!=P)
  # HP
  if(x["P"]>=x["M"] && x["H"]>=x["P"] && x["H"]> x["M"] && !x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.05" # H>=P>=M
  if(x["M"]>=x["P"] && x["H"]> x["P"] && x["H"]>=x["M"] && !x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.06" # H>=M>=P
  if(x["P"]> x["M"] && x["H"]>=x["P"] && x["H"]> x["M"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.07" # H>=P>M
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["H"]> x["M"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.08" # P>=H>M(P!=H)
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["H"]>=x["M"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.09" # H>=M>P
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.10" # M>=H>P(M!=H)
  # MP
  if(x["P"]> x["M"] && x["P"]>=x["H"] && x["H"]>=x["M"] &&  x["PvsM"] && !x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.11" # P>=H>=M
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.12" # P>H>M
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.13" # M>H>P
  if(x["M"]> x["P"] && x["H"]>=x["P"] && x["M"]>=x["H"] &&  x["PvsM"] && !x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.14" # M>=H>=P
  # LP
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["M"]> x["H"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.15" # M>H>=P(H!=P)
  if(x["M"]> x["P"] && x["P"]>=x["H"] && x["M"]> x["H"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.16" # M>P>=H
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.17" # P>H>=M(H!=M)
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["M"]>=x["H"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.18" # P>M>=H
  if(x["M"]>=x["P"] && x["P"]>=x["H"] && x["M"]> x["H"] && !x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.19" # M>=P>=H
  if(x["P"]>=x["M"] && x["P"]> x["H"] && x["M"]>=x["H"] && !x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.20" # P>=M>=H
  # BLP
  if(x["M"]> x["P"] && x["P"]> x["H"] && x["M"]> x["H"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.21" # M>=P>H(M!=P) 
  if(x["P"]>=x["M"] && x["P"]> x["H"] && x["M"]> x["H"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.22" # P>=M>H
  if(x["M"]> x["P"] && x["P"]> x["H"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.23" # M>P>H
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.24" # P>M>H
  # Criteria 2
  # adjust by H vs MPV
  if(!is.na(classCode) && classCode == "ClassNo.08" && x["H"] <= x["MPV"]) classCode <- "ClassNo.11" # P>=H>M (HP) change to P>=H>=M (MP)
  if(!is.na(classCode) && classCode == "ClassNo.10" && x["H"] <= x["MPV"]) classCode <- "ClassNo.14" # M>=H>P (HP) change to M>=H>=P (MP)
  if(!is.na(classCode) && classCode == "ClassNo.15" && x["H"] >= x["MPV"]) classCode <- "ClassNo.14" # M>H>=P (LP) change to M>=H>=P (MP)
  if(!is.na(classCode) && classCode == "ClassNo.17" && x["H"] >= x["MPV"]) classCode <- "ClassNo.11" # P>H>=M (LP) change to P>=H>=M (MP)
  #
  return(classCode)
}

## functions obtained from coolmap
cal_zScore <- function(x){
  M <- rowMeans(x, na.rm=TRUE)
  nsamples <- ncol(x)
  DF <- nsamples - 1L
  IsNA <- is.na(x)
  if(any(IsNA)) {
    mode(IsNA) <- "integer"
    DF <-  DF - rowSums(IsNA)
    DF[DF==0L] <- 1L
  }
  x <- x-M
  V <- rowSums(x^2L, na.rm=TRUE) / DF
  x <- x / sqrt(V+0.01)
  return(x)
}

#####################
### main pipeline ###
#####################
# read conf
cat("read configure file\n")
conf <- read.table(file_pedigree_conf, header = F, stringsAsFactors = F)
head(conf)
# global option from configure
maternal <- conf[conf[,1] == "Maternal", 2]
paternal <- conf[conf[,1] == "Paternal", 2]
hybrid   <- conf[conf[,1] == "Hybrid", 2]
sig.PvsM <- conf[conf[,1] == "PvsM", 2]
sig.PvsH <- conf[conf[,1] == "PvsH", 2]
sig.MvsH <- conf[conf[,1] == "MvsH", 2]

# read data
cat("read expr and significant data\n")
expr_sig <- read.table(file_pedigree_data, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_sig) <- expr_sig[,1]
head(expr_sig)

# check
cat("check configure and data colname\n")
if(any(!maternal %in% colnames(expr_sig)) ) stop("Sample name of maternal is not matched\n")
if(any(!paternal %in% colnames(expr_sig)) ) stop("Sample name of paternal is not matched\n")
if(any(!hybrid %in% colnames(expr_sig)) )   stop("Sample name of hybrid is not matched\n")
if(any(!sig.PvsM %in% colnames(expr_sig)) ) stop("Name of PvsM comparison is not matched\n")
if(any(!sig.PvsH %in% colnames(expr_sig)) ) stop("Name of PvsH comparison is not matched\n")
if(any(!sig.MvsH %in% colnames(expr_sig)) ) stop("Name of MvsH comparison is not matched\n")
cat("check is OK\n")

# construct sample list
cat("construct sample list\n")
sample_list <- data.frame(
  sample = c(maternal, paternal, hybrid), 
  group  = c(
    rep("M", length(maternal) ), 
    rep("P", length(paternal) ), 
    rep("H", length(hybrid) )
  ),
  check.names = F,
  stringsAsFactors = F
)
print(sample_list)

###############################################
# pre-process data
###############################################
# calculate mean of group
cat("calculate mean of group\n")
group_name <- unique(sample_list$group)
mtx.avg <- get_avg(as.matrix(expr_sig[,sample_list$sample] ), group_name, sample_list)
head(mtx.avg)

#
cat("calculate MPV\n")
MPV <- rowMeans(mtx.avg[,c("M", "P")])
MPV <- as.matrix(MPV)
colnames(MPV) <- c("MPV")
head(MPV)

#
cat("get a matrix of significance information\n")
mtx.sig <- as.matrix(expr_sig[,c(sig.PvsM, sig.PvsH, sig.MvsH)])
colnames(mtx.sig) <- c("PvsM", "PvsH", "MvsH")
head(mtx.sig)

# cbind in one
avg_sig_matrix <- cbind(mtx.avg, mtx.sig, MPV)
head(avg_sig_matrix)
rm(mtx.avg, mtx.sig, MPV)

###############################################
# expression pattern
###############################################
#
cat("identify expression pattern and add ClassCode\n")
classCodes <- apply(avg_sig_matrix, 1, identify_express_pattern)
expr_sig$ClassCode <- classCodes
head(expr_sig)

# add exprPattern and superclass information
cat("add exprPattern, diclass and superclass information\n")
expr_sig$ExprPattern    <- apply(as.matrix(classCodes), 1, function(x) ifelse(is.na(x), NA, as.character(classCode_list[x,]$ExprPattern) ) )
expr_sig$DiclassCode    <- apply(as.matrix(classCodes), 1, function(x) ifelse(is.na(x), NA, as.character(classCode_list[x,]$DiclassCode) ) )
expr_sig$DiclassPattern <- apply(as.matrix(classCodes), 1, function(x) ifelse(is.na(x), NA, as.character(classCode_list[x,]$DiclassPattern) ) )
expr_sig$SuperClass     <- apply(as.matrix(classCodes), 1, function(x) ifelse(is.na(x), NA, as.character(classCode_list[x,]$Superclass) ) )
head(expr_sig)

# order by classCode
cat("order by classCode\n")
expr_sig <- expr_sig[order(expr_sig$ClassCode),]
head(expr_sig)

# write out
cat("write out\n")
write.table(expr_sig, file_expr_pattern, sep = "\t", quote = F, row.names = F)

###############################################
# summary
###############################################
# total number of entries
totalNum <- nrow(expr_sig)

# count class
class_counts <- c()
for(classCode in rownames(classCode_list) ){
  count <- length( which(expr_sig$ClassCode == classCode) )
  if(length(class_counts) > 0)  class_counts <- c(class_counts, count)
  if(length(class_counts) == 0) class_counts <- c(count)
}
classCode_list$count <- class_counts
classCode_list$rate  <- class_counts / totalNum
classCode_list

# count diclass
diclass_counts <- c()
for(diclassCode in rownames(diClass_list) ){
  count <- length( which(expr_sig$DiclassCode == diclassCode) )
  if(length(diclass_counts) > 0)  diclass_counts <- c(diclass_counts, count)
  if(length(diclass_counts) == 0) diclass_counts <- c(count)
}
diClass_list$count <- diclass_counts
diClass_list$rate  <- diclass_counts / totalNum
diClass_list

# count superclass
superClass_counts <- c()
for(superClass in superClass_list$SuperClass){
  count <- length( which(expr_sig$SuperClass == superClass) )
  if(length(superClass_counts) > 0)  superClass_counts <- c(superClass_counts, count)
  if(length(superClass_counts) == 0) superClass_counts <- c(count)
}
superClass_list$count <- superClass_counts
superClass_list$rate  <- superClass_counts / totalNum
superClass_list

# draw count of class combined diclass and superclass
# ggplot of class
p1 <- ggplot(data = classCode_list, aes(x = ExprPattern, y = count) )
p1 <- p1 + geom_bar(stat = "identity")
p1 <- p1 + facet_wrap(~ Superclass, ncol = 5, scales = "free_x")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
p1 <- p1 + theme(axis.title = element_text(size = rel(2.5) ) )
p1 <- p1 + theme(axis.text = element_text(size = rel(1.6) ) )
p1 <- p1 + theme(strip.text = element_text(size = rel(2.2) ) )
p1 <- p1 + xlab("Expression Pattern of 24 classes") + ylab("Count")

# ggplot of diclass
p2 <- ggplot(data = diClass_list, aes(x = DiclassPattern, y = count) )
p2 <- p2 + geom_bar(stat = "identity")
p2 <- p2 + facet_wrap(~ Superclass, ncol = 5, scales = "free_x")
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
p2 <- p2 + theme(axis.title = element_text(size = rel(2.5) ) )
p2 <- p2 + theme(axis.text = element_text(size = rel(2.0) ) )
p2 <- p2 + theme(strip.text = element_text(size = rel(2.2) ) )
p2 <- p2 + xlab("Expression Pattern of 12 diclasses") + ylab("Count")
#p2 <- p2 + theme(axis.title.y = element_blank() )

# ggplot of superClass
p3 <- ggplot(data = superClass_list, aes(x = SuperClass, y = count) )
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + theme(axis.title = element_text(size = rel(2.5) ) )
p3 <- p3 + theme(axis.text = element_text(size = rel(2.2) ) )
p3 <- p3 + xlab("Superclass") + ylab("Count")

# open device
total_width  <- 12
total_height <- 18
tiff(filename = tiff_count_barplot,
     width = total_width,
     height = total_height,
     res = 1200,
     units = "in",
     compression = "lzw"
)

# width and height of each grid
grid_height  <- c(6, 6, 6)

# use gridExtra
gridExtra::grid.arrange(
  p1,  p2, p3,
  nrow = 3,
  heights = grid_height
)

# close device
dev.off()

# draw rate of class combined diclass and superclass
# ggplot of class
p1 <- ggplot(data = classCode_list, aes(x = ExprPattern, y = rate) )
p1 <- p1 + geom_bar(stat = "identity")
p1 <- p1 + facet_wrap(~ Superclass, ncol = 5, scales = "free_x")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
p1 <- p1 + theme(axis.title = element_text(size = rel(2.5) ) )
p1 <- p1 + theme(axis.text = element_text(size = rel(1.6) ) )
p1 <- p1 + theme(strip.text = element_text(size = rel(2.2) ) )
p1 <- p1 + xlab("Expression Pattern of 24 classes") + ylab("Rate")

# ggplot of diclass
p2 <- ggplot(data = diClass_list, aes(x = DiclassPattern, y = rate) )
p2 <- p2 + geom_bar(stat = "identity")
p2 <- p2 + facet_wrap(~ Superclass, ncol = 5, scales = "free_x")
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
p2 <- p2 + theme(axis.title = element_text(size = rel(2.5) ) )
p2 <- p2 + theme(axis.text = element_text(size = rel(2.0) ) )
p2 <- p2 + theme(strip.text = element_text(size = rel(2.2) ) )
p2 <- p2 + xlab("Expression Pattern of 12 diclasses") + ylab("Rate")
#p2 <- p2 + theme(axis.title.y = element_blank() )

# ggplot of superClass
p3 <- ggplot(data = superClass_list, aes(x = SuperClass, y = rate) )
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + theme(axis.title = element_text(size = rel(2.5) ) )
p3 <- p3 + theme(axis.text = element_text(size = rel(2.2) ) )
p3 <- p3 + xlab("Superclass") + ylab("Rate")

# open device
total_width  <- 12
total_height <- 18
tiff(filename = tiff_rate_barplot,
     width = total_width,
     height = total_height,
     res = 1200,
     units = "in",
     compression = "lzw"
)

# width and height of each grid
grid_height  <- c(6, 6, 6)

# use gridExtra
gridExtra::grid.arrange(
  p1,  p2, p3,
  nrow = 3,
  heights = grid_height
)

# close device
dev.off()

# write table of summary
# data frame of class
data_class_list <- data.frame(
  name  = paste(rownames(classCode_list), classCode_list$ExprPattern, sep = "."),
  count = classCode_list$count,
  check.names = F
)

# add number of missing data
numOfNA <- length(which(is.na(expr_sig$ClassCode) ) )
data_class_list <- rbind(data_class_list, data.frame(t(c(name = "Ambiguous", count = numOfNA))))

# add total number
data_class_list <- rbind(data_class_list, data.frame(t(c(name = "Total", count = totalNum))))

# add rate
data_class_list$rate <- sprintf("%.4f%%", 100*as.numeric(data_class_list$count)/totalNum)
data_class_list

# write out
write.table(data_class_list, file_class_count, quote = F, sep = "\t", row.names = F)

# data frame of diclass
data_diclass_list <- data.frame(
  name  = paste(diClass_list$DiclassCode, diClass_list$DiclassPattern, sep = "."),
  count = diClass_list$count,
  check.names = F
)
data_diclass_list <- rbind(data_diclass_list, data.frame(t(c(name = "Ambiguous", count = numOfNA))))
data_diclass_list$rate <- sprintf("%.4f%%", 100*as.numeric(data_diclass_list$count)/totalNum)
data_diclass_list

# write out
write.table(data_diclass_list, file_diclass_count, quote = F, sep = "\t", row.names = F)

#
superClass_list$rate <- sprintf("%.4f%%", 100*as.numeric(superClass_list$count)/totalNum)
superClass_list

# write out
write.table(superClass_list, file_super_count, quote = F, sep = "\t", row.names = F)

