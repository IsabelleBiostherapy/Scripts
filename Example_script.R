

---
  title: "Clinical Study Overview"
subtitle: "Isabelle Franco Moscardini - 05/09/2023"
output:
  html_document:
  toc: yes
df_print: paged
html_notebook:
  toc: yes
---
  
  

```{r setup, include=FALSE}

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load packages
pkgs <- c('pracma','writexl','readxl','lattice','ggplot2','ggrepel','dplyr',"tidyr", 'RColorBrewer', 'tmod', 'reshape2', 'rJava','tibble','survival','mgcv','mixOmics','ComplexHeatmap','circlize','Biobase','genefilter')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)

# Function to collapse probes into one single gene symbol signal
collapse.rows <- function(expr, probe.col, gene.col, data.table=F, method=c("maxMean", "minMean", "colMean", "colMedian")){
  if(length(grep('data.table', installed.packages())) == 0){
    install.packages('data.table')
    require(data.table)
  }else if(length(grep('data.table', search())) == 0){
    suppressPackageStartupMessages(require(data.table))
  }
  
  if (probe.col == "rownames"){
    expr <- data.table(expr, keep.rownames=T)
    setnames(expr, "rn", "rownames")
  }else{
    expr <- data.table(expr)
  }
  
  if(method=="maxMean" | method=="minMean"){
    expr[, rowmean := rowMeans(.SD[, !c(probe.col, gene.col), with=F])]
    if(method=="maxMean"){
      res <- expr[order(rowmean, decreasing=T)][, .SD[1], by=gene.col][, rowmean:=NULL]
    }
    else if(method=="minMean"){
      res <- expr[order(rowmean, decreasing=T)][, .SD[.N], by=gene.col][, rowmean:=NULL]
    }
  }
  else if(method=="colMean"){
    res <- expr[, lapply(.SD[, !c(probe.col), with=F], mean), by=gene.col]
  }
  else if(method=="colMedian"){
    res <- expr[, lapply(.SD[, !c(probe.col), with=F], median), by=gene.col]  
  }
  else stop("method must be 'maxMean', 'minMean', 'colMean' or 'colMedian'\n")
  
  if(!data.table){
    return(data.frame(res))
  }else{ return(res[]) }
}



```

### Section 1 data preprocessing

```{r}

## Filter 25% lowest variance
exprs <- read.delim("Data/Dummy_name.txt")
colnames(exprs) <- gsub(".sst.rma.gene.full.chp.Signal", "", colnames(exprs))
rownames(exprs) <- exprs$ID
exprs$ID <- NULL
eset <- ExpressionSet(assayData=as.matrix(exprs))
eset <- varFilter(eset, var.func=IQR, var.cutoff=0.25, filterByQuantile=TRUE)
exprs <- as.data.frame(exprs(eset))
dim(exprs) # 15985   159
head(exprs, 6)


# Load info samples
info <- read.delim("Data/Dummy_name1.txt")
info$File.Name <- gsub(".CEL", "", info$File.Name)
info$File.Name <- gsub("\\-", ".", info$File.Name)
info$TP <- gsub(".*_", "", info$File.Name) # Add time points
info$Group <- gsub(".*_", "", info$Condition.2) # Add group

# Same samples in expression and pheno
inter <- intersect(colnames(exprs), info$File.Name)
exprs <- exprs[, inter]
rownames(info) <- info$File.Name
info <- info[inter, ]
identical(colnames(exprs), rownames(info)) # TRUE
dim(exprs) # 15985   158

# Featuredata
feature <- read.csv("Data/Dummy_name2.csv")
colnames(feature)[1] <- "Probes"
feature$probeset_id <- NULL
exprs <- rownames_to_column(exprs, "Probes")
panel <- full_join(feature, exprs)
panel <- panel[complete.cases(panel$A02.01.04_V0), ]
dim(panel) # 15985   160
panel$Symbol <- gsub(" ", "", panel$Symbol)
head(panel, 6)

# Merge information of probes in a single gene symbol value 
panel <- collapse.rows(expr = panel, probe.col = "Probes", gene.col = "Symbol", method = "maxMean")
head(panel)
dim(panel) # 14752   160

# Create a table with specific Gene symbol for duplicated, excluding TSUnmmaped (but leaving DUX5, DUX3, DUX1)
info_duplicated <- panel[, 1:2]
info_duplicated$Probes <- ifelse(grepl("TSU", info_duplicated$Probes) & !grepl("DUX5|DUX3|DUX1", info_duplicated$Symbol), 0, info_duplicated$Probes) # Mark with 0 what should be removed
info_duplicated <- info_duplicated[info_duplicated$Probes != 0, ] # Remove what was marked to be removed
dim(info_duplicated) # 14688     3

# Generate the matrix expression with the gene symbols generate and use this one in the next steps
panel <- full_join(info_duplicated, panel[,c(2:160)])
dim(panel) #  14752   160
head(panel, 6)  
panel <- panel[complete.cases(panel$Symbol), ]
dim(panel) # 14688   161

## 
pheno <- as.data.frame(colnames(panel))
colnames(pheno) <- "SampleName"
pheno$Class <- gsub(".*_", "", pheno$SampleName)
dim(pheno) # 160   2

pheno <- pheno[3:160, ]
rownames(pheno) <- pheno$SampleName
rownames(panel) <- panel$Symbol
panel$Symbol <- NULL
panel$Probes <- NULL
head(panel, 6)

panel_names <- rownames_to_column(panel, "Genes")
head(panel_names, 6)
#write_xlsx(panel_names, "Interm/Expression_collapsed_MaxMean.xlsx")

```

### Section 2: generate the fold-changes for each patient

```{r}

## Fold-change
exp_V0 <- panel[,grepl("V0", colnames(panel))]
colnames(exp_V0) <- gsub("_V0", "", colnames(exp_V0))
colnames(exp_V0) <- gsub("^....", "", colnames(exp_V0))
exp_V0 <- as.data.frame(apply(exp_V0, 2, as.numeric))
rownames(exp_V0) <- rownames(panel)

exp_V2 <- panel[,grepl("V2", colnames(panel))]
colnames(exp_V2) <- gsub("_V2", "", colnames(exp_V2))
colnames(exp_V2) <- gsub("^....", "", colnames(exp_V2))
exp_V2 <- as.data.frame(apply(exp_V2, 2, as.numeric))
rownames(exp_V2) <- rownames(panel)

exp_V3 <- panel[,grepl("V3", colnames(panel))]
colnames(exp_V3) <- gsub("_V3", "", colnames(exp_V3))
colnames(exp_V3) <- gsub("^....", "", colnames(exp_V3))
exp_V3 <- as.data.frame(apply(exp_V3, 2, as.numeric))
rownames(exp_V3) <- rownames(panel)

exp_V5 <- panel[,grepl("V5", colnames(panel))]
colnames(exp_V5) <- gsub("_V5", "", colnames(exp_V5))
colnames(exp_V5) <- gsub("^....", "", colnames(exp_V5))
exp_V5 <- as.data.frame(apply(exp_V5, 2, as.numeric))
rownames(exp_V5) <- rownames(panel)

# V2
interV2V0 <- intersect(colnames(exp_V0), colnames(exp_V2))
exp_V0 <- exp_V0[, interV2V0]
exp_V2 <- exp_V2[, interV2V0]
identical(colnames(exp_V0), colnames(exp_V2)) # TRUE
foldv2v0 <- exp_V2-exp_V0

# V3
interV3V0 <- intersect(colnames(exp_V0), colnames(exp_V3))
exp_V0 <- exp_V0[, interV3V0]
exp_V3 <- exp_V3[, interV3V0]
identical(colnames(exp_V0), colnames(exp_V3)) # TRUE
foldv3v0 <- exp_V3-exp_V0

# V5
interV5V3 <- intersect(colnames(exp_V3), colnames(exp_V5))
exp_V3 <- exp_V3[, interV5V3]
exp_V5 <- exp_V5[, interV5V3]
identical(colnames(exp_V3), colnames(exp_V5)) # TRUE
foldv5v3 <- exp_V5-exp_V3


# Panel with all the FC together
colnames(foldv2v0) <- paste0(colnames(foldv2v0), "_V2V0")
colnames(foldv3v0) <- paste0(colnames(foldv3v0), "_V3V0")
colnames(foldv5v3) <- paste0(colnames(foldv5v3), "_V5V3")

foldv2v0 <- rownames_to_column(foldv2v0, "Symbol")
foldv3v0 <- rownames_to_column(foldv3v0, "Symbol")
foldv5v3 <- rownames_to_column(foldv5v3, "Symbol")

panel_fc <- Reduce(full_join, list(foldv2v0, foldv3v0, foldv5v3))

rownames(panel_fc) <- panel_fc$Symbol
panel_fc$Symbol <- NULL
head(panel_fc, 6)

```

### Section 3: Filtering groups, fold-Changes and performing analysis 

#### Section 3.1: Fold-Change V2V0 

```{r}
# prepare pheno
pheno <- as.data.frame(colnames(panel_fc))
pheno$ID <- gsub("_.*","",pheno$`colnames(panel_fc)`)
info$X <- gsub("^....","", info$File.Name)
info$X <- gsub("_.*","", info$X)
pheno$Group <- info$Group[match(pheno$ID, info$X)]
pheno$Visit <- gsub(".*_", "", pheno$`colnames(panel_fc)`)
head(pheno, 6)

## Filter by Fold Change
head(foldv2v0, 6)
fold <- foldv2v0

phenoV2 <- pheno
phenoV2 <- pheno[pheno$Visit == "V2V0" & pheno$Group == "BA", ]

rownames(fold) <- fold$Symbol
fold <- fold[, phenoV2$`colnames(panel_fc)`]
head(pheno, 6)

# fold[abs(fold) < 0.5] <- 0 # everything that has a absolute value of logFC less than 0.5 becomes 0
# fold <- fold[rowSums(fold) != 0, ] # Remove genes that are 0 in all samples
dim(fold) # 14677    17
head(fold, 6)

```


#### Section 3.2: Heatmap and Euclidean distances

```{r}
logfc_v2v0 <- rownames_to_column(fold, "Symbol")


# Select only columns for the heatmap
merged_heat <- logfc_v2v0[,c(1,7:24)]

head(merged_heat, 6)

```


```{r}
# Choose the colors of the heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "orange1")) 

merged_heat <- as.data.frame(merged_heat)
rownames(merged_heat) <- make.names(merged_heat$Symbol, unique = TRUE)
merged_heat$Symbol <- NULL

merged_heat <- merged_heat[order(merged_heat$Profile, decreasing = FALSE), ]

merged <- rownames_to_column(merged_heat, "Genes")
merged$Genes <- gsub("\\..*", "", merged$Genes)
merged <- merged[!duplicated(merged$Genes), ]
rownames(merged) <- merged$Genes
merged$Genes <- NULL

heat <- Heatmap(merged, col = col_fun, cluster_rows = FALSE, 
                row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 8))
plot(heat)


```


```{r}
# Extract Euclidean distances
distances <- dist(t(merged))
print(distances)

```



#### Section 5: Correlations with clinical data

```{r}
# Clinical inputs
id <- read.delim("Data/Dummy_name3.txt")
id$File.Name <- gsub(".CEL", "", id$File.Name)
id$File.Name <- gsub("\\-", ".", id$File.Name)
id$X <- gsub("^....", "", id$File.Name)


clinics <- read_excel("Data/Dummy_name4.xlsx", sheet = 4)
clinics[,3:6] <- apply(clinics[,3:6], 2, function(x) gsub(" .*", "", x))
clinics[,3:6] <- apply(clinics[,3:6], 2, as.numeric)

# Sum to score
clinics$score <- clinics$Dummy1 + clinics$Dummy2 +clinics$Dummy3 +clinics$Dummy4
head(clinics, 6)

clinicV0 <- as.data.frame(clinics[clinics$Visit == "Visit 0", ])
clinicV2 <- as.data.frame(clinics[clinics$Visit == "Visit 2", ])
interclinic <- intersect(clinicV0$`Patient Number`, clinicV2$`Patient Number`)
rownames(clinicV0) <- clinicV0$`Patient Number`
rownames(clinicV2) <- clinicV2$`Patient Number`
clinicV0 <- clinicV0[interclinic, ]
clinicV2 <- clinicV2[interclinic, ]
head(clinicV0, 6)

clinicV0 <- clinicV0[,3:7]
clinicV2 <- clinicV2[,3:7]
identical(rownames(clinicV0), rownames(clinicV2)) #TRUE
identical(colnames(clinicV0), colnames(clinicV2)) # TRUE
clinic_V2V0 <- clinicV2-clinicV0
clinic_V2V0$ID <- id$X[match(rownames(clinic_V2V0), gsub("_.*", "", id$Condition))]
clinic_V2V0$ID <- gsub("_.*", "_V2V0", clinic_V2V0$ID)
clinic_V2V0 <- clinic_V2V0[complete.cases(clinic_V2V0$ID), ]
rownames(clinic_V2V0) <- clinic_V2V0$ID
clinic_V2V0$ID <- NULL

clinic_V2V0$Volunteer <- NULL
head(clinic_V2V0, 6)

```

#### Section 5.1: Correlations V2V0 

```{r}
A <- logfc_v2v0
B <- t(clinic_V2V0)
head(A, 6)
head(B, 6)

rownames(A) <- A$Symbol
A$Symbol <- NULL

identical(colnames(A), colnames(B)) # TRUE

corr_clinic <- c() 
for(i in 1:5){ 
  print(rownames(B)[i]) 
  for(j in 1:nrow(A)){
    x   <- B[i,] 
    y   <- A[j,] 
    cor <- cor.test(as.numeric(x),as.numeric(y),method="spearman", exact = FALSE) 
    
    corr_clinic  <- rbind(corr_clinic ,cbind(rownames(A)[j],rownames(B)[i],cor$estimate,cor$p.value)) 
  }
  
}


colnames(corr_clinic) <- c("Gene", "Test", "R", "Pval")
corr_clinic <- as.data.frame(corr_clinic)
corr_clinic[,3:4] <- apply(corr_clinic[,3:4],2,as.numeric)

corr_clinic_filt <- corr_clinic[corr_clinic$Pval < 0.05, ]
corr_clinic_filt <- corr_clinic_filt[abs(corr_clinic_filt$R)  >  0.4, ]
head(corr_clinic_filt, 6)

corr_clinic_filt <- corr_clinic_filt[corr_clinic_filt$Test == "Score", ]
View(corr_clinic_filt)

```


```{r}
# Save by Symptom
for (i in unique(corr_clinic_filt$Test)) {
  tested <- corr_clinic_filt[corr_clinic_filt$Test == i, ]
  tested$Direction <- ifelse(tested$R > 0, 1, 2) 
  tested <- tested[,c(1,5,2:4)]
  tested$Gene <- gsub("\\.", "_", tested$Gene)
  writexl::write_xlsx(tested, paste0("Results/Correlations/Primary_endpoints/Filtered_spearman_corr_", i, ".xlsx"))
}


#write_xlsx(corr_clinic_filt, "Results/Correlations/corr_clinic_filt_primaryendpoints.xlsx")


```
