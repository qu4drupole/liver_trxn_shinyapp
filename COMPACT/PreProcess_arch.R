
library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)
library(tidyverse)
library(ggplot2)


#######################
#GSE20426
#######################

load("Data/GSE20426.RData")

#This can be a function
#vvvvvvvvvvvvvvvvvvvvvvvv
gset <- getGEO("GSE20426", GSEMatrix =TRUE, AnnotGPL=F)
gset <- gset[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
#^^^^^^^^^^^^^^^^^^^^^^^^

#make time variable
gset$Time <- factor(as.numeric(gset$`time (hours):ch1`))

#make trmt variable (in this case I need to combine treatments (surgery type and age))
gset$trmt1 <- character(length = nrow(pData(gset)))
gset$trmt1[gset$`treatment:ch1`=="none"] = "none"
gset$trmt1[gset$`treatment:ch1`=="partial hepatectomy"] = "PHx"
gset$trmt1[gset$`treatment:ch1`=="sham partial hepatectomy"] = "sham"

gset$trmt2 <- character(length = nrow(pData(gset)))
gset$trmt2[gset$`age:ch1` == levels(factor(gset$`age:ch1`))[1]] = "old"
gset$trmt2[gset$`age:ch1` == levels(factor(gset$`age:ch1`))[2]] = "young"

gset$Trmt = factor(sapply(1:nrow(pData(gset)), function(x) paste(gset$trmt1[x], gset$trmt2[x], sep = "_")))

#assign a title
experimentData(gset)@title <- "Hepatic gene expression during liver regeneration in response to partial hepatectomy: late time points (24h, 38h, 48h)"

save(gset,file="Data/GSE20426.RData")

#######################
#GSE100616
#######################

gset <- getGEO("GSE100616", GSEMatrix =TRUE, AnnotGPL=F)
gset <- gset[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

varLabels(gset)
#Weird labeling...
#make time variable
gset$Time <- factor(c(0,0,0,0,2,2,2,2,6,6,6,6))

#make trmt variable (in this case I need to combine treatments (surgery type and age))
gset$Trmt = factor(c('wt','wt','ko','ko','wt','wt','ko','ko','wt','wt','ko','ko'))

#assign a title
experimentData(gset)@title <- "Gene expression profiling of Pml wt and Pml KO mice liver with acetaminophen (apap) overdose (300mg/kg) i.p."

save(gset,file="Data/GSE100616.RData")

#######################
#GSE33785
#######################

gset <- getGEO("GSE33785", GSEMatrix =TRUE, AnnotGPL=T)
gset <- gset[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

varLabels(gset)
gset$Trmt <- factor(gset$`diet:ch1`)
gset$Time <- factor(sapply(pData(gset)[["timepoint:ch1"]],function(x) strsplit(x, "\\s+")[[1]][1], USE.NAMES = F))
experimentData(gset)@title <- "Gene expression in the chronic ethanol-treated rat liver during liver regeneration"
save(gset,file="Data/GSE33785.RData")

#####################################################################################################################


#Maybe the expression set should be object that can store other info like treatment names--instead of having to parse colnames

#Collapse expression into averages
#'g' is a gse object
eDcollapse <- function(g){
  #browser()
  #count unique trmt*time combinations
  cols = 0
  for(i in levels(g$Time)){
    cols = cols + length(unique(g$Trmt[g$Time==i]))
  }
  res.avg <- data.frame(matrix(nrow = nrow(fData(g)), ncol = cols))
  ci = 1
  #the double for loops can be replaced with "interaction" function
  for(i in levels(g$Time)){
    for(j in unique(g$Trmt[g$Time==i])){
      res.avg[,ci] <- rowMeans(exprs(g[,which(g$Time == i & g$Trmt == j)]))
      colnames(res.avg)[ci] = paste(j, as.character(i), sep = "//")
      ci = ci + 1
    }
  }
  rownames(res.avg) <- rownames(fData(g))
  return(res.avg)
}

res.test <- eDcollapse(gset)

eDcontrast <- function(eD, bl){
  #browser()
  #need to provide a warning if baseline isn't balanced with other Trmt*Time combinations
  bal.table <- table(strsplit2(names(eD),"//")[,2])
  if(mean(bal.table) == bal.table[as.character(bl)]){
    a <- eD[,which(strsplit2(names(eD),"//")[,2] == bl)]
    b <- eD[,which(strsplit2(names(eD),"//")[,2] != bl)]
    ct.matrix <- data.frame(matrix(nrow = nrow(eD), ncol = ncol(b)))
    colnames(ct.matrix) <- colnames(b)
    for(trmt.name in strsplit2(names(a),"//")[,1]){
      bl.data <- a[,grep(trmt.name, colnames(a))]
      trmt.key <- grep(trmt.name, strsplit2(names(b),"//")[,1])
      ct.matrix[, trmt.key] <- b[, trmt.key] - bl.data
    }
  }
  else{
    #The experiment is unbalanced, proceed with caution
    #hopefully the problem would be a generalized baseline
    if(bal.table[as.character(bl)] < mean(bal.table)){
      #find the differentiating trmt--this situation is likely caused by multilevel experiment design
      bl.names <- strsplit2(names(eD)[which(strsplit2(names(eD),"//")[,2]==bl)],"//")[,1]
      if(strsplit2(bl.names[1],"_")[,1] == strsplit2(bl.names[2],"_")[,1]){
        ct.key <- 2
      }
      else{ct.key <- 1}
      a <- eD[,which(strsplit2(names(eD),"//")[,2] == bl)]
      b <- eD[,which(strsplit2(names(eD),"//")[,2] != bl)]
      ct.matrix <- data.frame(matrix(nrow = nrow(eD), ncol = ncol(b)))
      colnames(ct.matrix) <- colnames(b)
      for(trmt.name in strsplit2(bl.names,"_")[,ct.key]){
        bl.data <- a[,grep(trmt.name, colnames(a))]
        trmt.key <- grep(trmt.name, strsplit2(names(b),"//")[,1])
        ct.matrix[, trmt.key] <- b[, trmt.key] - bl.data
      }
    }
    else{
      #you're fucked...too complicated
      #or maybe just drop treatments/columns that make things difficult
      #for now...just throw an exception
      }
  }
  rownames(ct.matrix) <- rownames(eD)
  return(ct.matrix)
}

res.contrast <- eDcontrast(res.test, 0)

eDdiscretize <- function(ctrst, minSet = 5, tn = 2){
  #browser()
  thresh <- log2(tn)
  res.unary <- data.frame(matrix(nrow = nrow(ctrst), ncol = ncol(ctrst)))
  colnames(res.unary) <- colnames(ctrst)
  rownames(res.unary) <- rownames(ctrst)
  res.unary <- apply(ctrst,2,function(x) {(abs(x)>thresh)*sign(x)})
  #make a data frame for treatment-specific patterns
  d.collapse <- data.frame(matrix(nrow = nrow(res.unary), ncol = length(unique(strsplit2(names(ctrst),"//")[,1]))))
  rownames(d.collapse) <- rownames(res.unary)
  colnames(d.collapse) <- unique(strsplit2(names(ctrst),"//")[,1])
  for(trmt in colnames(d.collapse)){
    a <- res.unary[,grep(trmt,colnames(res.unary))]
    a <- a[,order(strsplit2(colnames(a),"//")[,2])]
    d.collapse[,trmt] <- apply(a,1,paste,collapse=",")
  }
  
  #count the pattern sets, drop rare patterns
  for(pattern in unique(unlist(c(d.collapse)))){
    if(pattern %in% unlist(c(d.collapse)) & sum(d.collapse == pattern) < minSet){
      d.collapse <- d.collapse[-which(unname(apply(d.collapse,1,function(x) sum(grepl(pattern,x))>0))),]
    }
  }
  res.discrete = list(collapse = d.collapse, trmt = colnames(d.collapse), patterns = unique(unlist(c(d.collapse))))
  res.discrete$patterns <- res.discrete$patterns[order(res.discrete$patterns)]
  return(res.discrete)
}

res.discrete <- eDdiscretize(res.contrast, minSet = 15)

#making the COMPACT matrix
#the input to this has to be the output of eDdiscretize
makeCOMPACT <- function(res, t1, t2){
  #browser()
  dat <- subset(res$collapse, select = c(t1,t2))
  compact.m <- data.frame(matrix(nrow = length(res$patterns), ncol = length(res$patterns)))
  rownames(compact.m) <- res$patterns
  colnames(compact.m) <- res$patterns
  for(i in colnames(compact.m)){
    for(j in rownames(compact.m)){
      #create some object to store the coordinates of these probe IDs
      probeIDs <- intersect(rownames(dat)[dat[,t1] == i],rownames(dat)[dat[,t2] == j])
      compact.m[i,j] = length(probeIDs)
    }
  }
  center.index <- str_count(res$patterns[1],",")+1
  compact.m[paste(as.character(rep(0,center.index)),collapse = ','),paste(as.character(rep(0,center.index)),collapse = ',')] <- 0
  #eventually make the return a list with the COMPACT matrix and intersection IDs
  return(compact.m)
}

cm <- makeCOMPACT(res.discrete,"PHx_young","PHx_old")

#make a heat map of COMPACT
#...using heatmap.2
heatmap.2(as.matrix(cm), 
          trace = 'none', 
          Colv = NULL,
          Rowv = NULL,
          # dendrogram = 'row',
          # srtCol = 90, 
          # adjCol = c(0.8,0.6), 
          # labCol = c(d1lab,d2lab),
          # cexCol = 1, 
          # labRow = NA,
          # margins = c(7,1),
          keysize = 1)

#...using ggplot
hm.data <- cm %>% rownames_to_column('t2') %>% gather(rownames(cm), value, -t2)
#NEED TO DOUBLE CHECK THAT NAMING!
names(hm.data) <- c("Carb","High_fat","value")
a<-"Carb"
b<-"High_fat"
ggplot(hm.data, aes_string(a, b)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = value)) +
  scale_fill_gradient(low = "white", high = "red") 


