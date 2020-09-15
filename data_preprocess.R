library(GEOquery)
library(limma)
library(Biobase)
library(FNN)
library(gplots)
library(RColorBrewer)
library(factoextra)
library(sva)


gset <- getGEO("GSE63742", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

##VVVV I'll need a better way to do this... VVVV
gsms <- "000XXXXXXXXXXXXXXXXXXXXX111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#thankfully all of this is fast
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol"))


#Need to change these so GSEs can be readily compared and organized
#i.e. all GSEs will have a "trmt" variable
varLabels(gset)
gset$Time <- factor(sapply(pData(gset)[["time:ch1"]],function(x) strsplit(x, "\\s+")[[1]][1], USE.NAMES = F))
gset$Treatment <- factor(pData(gset)[["treatment:ch1"]])
a <- grep("///", fData(gset)$Gene.symbol)
fData(gset)$Gene.symbol[a] <- sapply(a, function(x) strsplit(fData(gset)$Gene.symbol[x],"///")[[1]][2])
experimentData(gset)@title <- "Expression data from rat liver tissue during liver regeneration after partial hepatectomy"
save(gset,file="Data/GSE63742.RData")
gset1 <- gset

###GSE33785
gset <- getGEO("GSE33785", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6247", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

gset$Treatment <- factor(gset$`diet:ch1`)
gset$Time <- factor(sapply(pData(gset)[["timepoint:ch1"]],function(x) strsplit(x, "\\s+")[[1]][1], USE.NAMES = F))
experimentData(gset)@title <- "Gene expression in the chronic ethanol-treated rat liver during liver regeneration"
#fData(gset)$Gene.symbol <- fData(gset)$`Gene symbol`
a <- grep("///", fData(gset)$Gene.symbol)
fData(gset)$Gene.symbol[a] <- sapply(a, function(x) strsplit(fData(gset)$Gene.symbol[x],"///")[[1]][2])
gset$ID <- gset$`animal id:ch1`
save(gset,file="Data/GSE33785.RData")
gset2 <- gset

###GSE97429
gset <- getGEO("GSE97429", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6247", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

unique(pData(gset)$`tissue:ch1`)
nrow(pData(gset))

gset$Treatment <- factor(gset$`treatment:ch1`)
gset$Time <- 0
gset$Time <- gset$Time + 24*(gset$`treatment:ch1` != "none")
gset$Time <- factor(gset$Time)
experimentData(gset)@title <- "Gene expression in the liver remnant is significantly affected by the size of partial hepatectomy - an experimental rat study"
fData(gset)$Gene.symbol <- sapply(fData(gset)$`gene_assignment`,function(x) {
  if(length(strsplit(x,"//")[[1]] > 1)){
    strsplit(x,"//")[[1]][2]}
  else{NA}
}
, USE.NAMES = F)
fData(gset)$Gene.symbol <- trimws(fData(gset)$Gene.symbol)

#there appears to be negative values in this dataset
head(exprs(gset)[,gset$Treatment=='none'])
#gset3 is RMA normalized (log2) and median transformed...whatever that last part means
hist(exprs(gset)[,1],breaks = 20)
hist(exprs(gset1)[,1],breaks = 20)
#The data is also centered...
median(exprs(gset1)[,1])
#7.09
#this is so sloppy but w/e...
hist(7+exprs(gset)[,1],breaks = 20)


exprs(gset) <- 7+exprs(gset)
save(gset,file="Data/GSE97429.RData")
gset3 <- gset

#GSE67022
gset <- getGEO("GSE67022", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6247", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

varLabels(gset)
gset$`time point:ch1`
a <- sapply(pData(gset)[["injected with:ch1"]],function(x) strsplit(x, "\\s+")[[1]][1], USE.NAMES = F)
for(i in grep(',',a)){
  a[i] <- substr(a[i],1,nchar(a[i])-1)
}
gset$Treatment <- factor(paste(a,gset$`molecule subtype:ch1`,sep="."))
gset$Time <- factor(gset$`time point:ch1`)
fData(gset)$Gene.symbol <- sapply(fData(gset)$`gene_assignment`,function(x) {
  if(length(strsplit(x,"//")[[1]] > 1)){
    strsplit(x,"//")[[1]][2]}
  else{NA}
}
, USE.NAMES = F)
fData(gset)$Gene.symbol <- trimws(fData(gset)$Gene.symbol)
experimentData(gset)@title <- "Regulation of Rat Hepatic Translation by mTOR"
save(gset,file="Data/GSE67022.RData")
gset4 <- gset

#for some reason gset4 can't be put in hm.data

grep(tT$Gene.symbol[8],fData(test1.cut)$Gene.symbol)
str(fData(test1.cut)$Gene.symbol)
str(fData(gset1)$Gene.symbol)
match(tT$Gene.symbol[8],fData(test1.cut)$Gene.symbol)
exprs(test1.cut)[13670,]

#ALSO treatment/level selection on gset4 appears to be inaccurate

###############################
#Making a combined variable
###############################
test.cut$dvar <- factor(match(paste(test.cut$Treatment,test.cut$Time),unique(paste(test.cut$Treatment,test.cut$Time))))
fData(test.cut)$logFC <- as.numeric(apply(exprs(test.cut),1,function(x) mean(x[test.cut$geo_accession[test.cut$dvar==2]]) - mean(x[test.cut$geo_accession[test.cut$dvar==1]])))



fData(test.cut)$`Gene symbol`[fData(test.cut)$`Gene symbol` == ""] <- NA
#fData(test.cut)$`Gene symbol`[which(fData(test.cut)$`Gene symbol` %in% tT$Gene.symbol)]
#tT2 <- subset(fData(test.cut)[which(fData(test.cut)$`Gene symbol` %in% tT$Gene.symbol),], select = c("Gene symbol","logFC"))

tT$logFC.d2 <- fData(test.cut)$logFC[match(tT$Gene.symbol,fData(test.cut)$`Gene symbol`)]

tT2 <- subset(tT, select = c("Gene.symbol","logFC","logFC.d2"))
p <- heatmap.2(as.matrix(tT2[,-1]), trace = 'none', dendrogram = 'row', labCol = c("dataset2","dataset1"),srtCol = 0, adjCol = c(0.5,0.6), cexCol = 3, labRow = NA)



###################################
#individual sample heat map
###################################

c1 <- which(gset1$Treatment == "partial hepatectomy" & gset1$Time == "0h")
c2 <- which(gset1$Treatment == "partial hepatectomy" & gset1$Time == "24h")
test1.cut <- gset1[,c(c1,c2)]

c1 <- which(gset2$Treatment == "Carbohydrate" & gset2$Time == 0)
c2 <- which(gset2$Treatment == "Carbohydrate" & gset2$Time == 24)
test2.cut <- gset2[,c(c1,c2)]
#c3 <- which(test2.cut$`animal id:ch1` %in% names(which(table(test2.cut$`animal id:ch1`)>1)))
if(any(table(test2.cut$ID)>1)){
  c3 <- which(test2.cut$ID %in% names(which(table(test2.cut$ID)>1)))
  test2.cut <- test2.cut[,c3]
}
test2.cut <- test2.cut[,c3]

#fData(test2.cut)$Gene.symbol <- fData(test2.cut)$`Gene symbol`
#match tT to test.cut
tT <- tT[tT$Gene.symbol != "",]
test1.sig <- exprs(test1.cut)[match(tT$Gene.symbol,fData(test1.cut)$Gene.symbol),]
test2.sig <- exprs(test2.cut)[match(tT$Gene.symbol,fData(test2.cut)$Gene.symbol),]

#combine and trim data
hm.data <- cbind(test1.sig,test2.sig)
a <- apply(hm.data,1,function(x) sum(!is.na(x)) == length(x))
hm.data <- hm.data[a,]
d1lab <- sapply(rownames(pData(test1.cut)), function(x) paste("d1",pData(test1.cut)[x,"Treatment"],pData(test1.cut)[x,"Time"],sep = "/"),USE.NAMES = F)
d2lab <- sapply(rownames(pData(test2.cut)), function(x) paste("d2",pData(test2.cut)[x,"Treatment"],pData(test2.cut)[x,"Time"],sep = "/"),USE.NAMES = F)
invisible(q <- heatmap.2(hm.data, 
                         trace = 'none', 
                         Colv = NULL,
                         dendrogram = 'row',
                         srtCol = 90, 
                         adjCol = c(0.8,0.6), 
                         labCol = c(d1lab,d2lab),
                         cexCol = 1, 
                         labRow = NA,
                         margins = c(7,1),
                         keysize = 1))

kl.data <- hm.data[p$rowInd,]
rownames(pData(test1.cut))[test1.cut$Time=="24h"]
a <- rowMeans(kl.data[,4:6])-rowMeans(kl.data[,1:3])
b <- rowMeans(kl.data[,10:12])-rowMeans(kl.data[,7:9])
c<-KL.divergence(a,b,100)

###################################
#individual sample heat map
###################################

rbind(subset(pData(test1.cut), select = c("Treatment","Time", "geo_accession")),
      subset(pData(test2.cut), select = c("Treatment","Time", "geo_accession")))

test2.cut.match <- exprs(test2.cut)[match(fData(test1.cut)$Gene.symbol,fData(test2.cut)$Gene.symbol),]
a <- rowMeans(test2.cut.match[,1:3])
b <- exprs(test2.cut.match)[,4:6]
#b <- rowMeans(exprs(test2.cut)[,4:6]) - rowMeans(exprs(test2.cut[,1:3]))
ma.data <- cbind(exprs(test1.cut),test2.cut.match)
par(mfrow=c(ceiling(ncol(ma.data)/3),3))
for (i in 1:ncol(ma.data)){
  plotMA(ma.data,i)
}

###################################
#Clustering
###################################

#basically need to start from the heat map data
#try clustering on treatment condition means
test1.cut$dvar <- factor(match(paste(test1.cut$Treatment,test1.cut$Time),unique(paste(test1.cut$Treatment,test1.cut$Time))))
test2.cut$dvar <- factor(match(paste(test2.cut$Treatment,test2.cut$Time),unique(paste(test2.cut$Treatment,test2.cut$Time))))
a <- rownames(pData(test1.cut))[test1.cut$dvar==1]
b <- rownames(pData(test1.cut))[test1.cut$dvar==2]
c <- rownames(pData(test2.cut))[test2.cut$dvar==1]
d <- rownames(pData(test2.cut))[test2.cut$dvar==2]

cl.data <- cbind(rowMedians(hm.data[,a]),rowMedians(hm.data[,b]),rowMedians(hm.data[,c]),rowMedians(hm.data[,d]))
#not doing any scaling...
#Trying k-means; not sure what nstart is doing...
wss <- sapply(2:15, function(x) {kmeans(cl.data,x,nstart = 5)$tot.withinss})
plot(2:15, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
#I think k=6 looks good
cl.test <- kmeans(cl.data,6,nstart = 5)
selcol2 <- brewer.pal(6,"Set1")
#make a new heat map, ordered the same as the "main" one, but add colors for the clusters
heatmap.2(cl.data[p$rowInd,], 
          trace = 'none', 
          Colv = NULL,
          Rowv = NULL,
          dendrogram = 'none',
          srtCol = 90, 
          adjCol = c(0.8,0.6), 
          cexCol = 1, 
          labRow = NA,
          margins = c(7,1),
          keysize = 1,
          RowSideColors = selcol2[cl.test$cluster[p$rowInd]])
#Trying with hclust
cl.fit <- hclust(dist(cl.data))
fviz_nbclust(cl.data, hcut, method = "wss")
cl.test <- cutree(cl.fit, k=5)
selcol <- brewer.pal(5,"Set1")
heatmap.2(cl.data[cl.fit$order,], 
          trace = 'none', 
          Colv = NULL,
          dendrogram = 'row',
          srtCol = 90, 
          adjCol = c(0.8,0.6), 
          #cexCol = 1, 
          #labRow = NA,
          lmat = rbind(c(0,0,4),c(3,1,2),c(5,0,0)),
          lhei = c(0.2,4,0.2),
          lwid = c(1,0.2,4),
          margins = c(2,20),
          key = F,
          RowSideColors = selcol[cl.test[cl.fit$order]],
          labRow = NA)
ggplot(data.frame(x=1:5,y=rep(1,5)),aes(x=x, y=y))+
  geom_bar(fill = selcol, stat = "identity", show.legend = F)+
  coord_equal()+
  xlab("Clusters")+
  geom_text(aes(label=1:5), size=10, nudge_y = -0.5)+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=20))
probeID <- rownames(hm.data)[cl.fit$order[cl.test==2]]
fData(test1.cut)[probeID,"Gene.symbol"]

#try clustering on fold change

###################################
#ComBat normalization
###################################

#need to match the datasets based on some common criteria--in this case gene symbols. Maybe also try Entrez/refseq/Ensembl/...

#collapse datasets
gset1.1 <- test1.cut[fData(test1.cut)$Gene.symbol != "",]
gset2.1 <- test2.cut[fData(test2.cut)$Gene.symbol != "",]
pData(gset1.1)$batch <- 1
pData(gset2.1)$batch <- 2
pheno_ <- rbind(subset(pData(gset1.1),select=c("Treatment", "Time", "batch")), subset(pData(gset2.1),select=c("Treatment","Time","batch")))
match.gl <- fData(gset2.1)$Gene.symbol[match(unique(fData(gset1.1)$Gene.symbol),fData(gset2.1)$Gene.symbol)]
match.gl <- match.gl[!is.na(match.gl)]
tT.trim <- tT[tT$Gene.symbol != "",]
exprs_1 <- cbind(exprs(gset1.1)[tT.trim$ID,],exprs(gset2.1)[match(tT.trim$Gene.symbol,fData(gset2.1)$Gene.symbol),])
a <- apply(exprs_1, 1, function(x) sum(!is.na(x)) == length(x))
exprs_1 <- exprs_1[a,]
match.gl <- match.gl[which(!(match.gl %in% tT.trim$Gene.symbol))]
exprs_2 <- cbind(exprs(gset1.1)[match(match.gl,fData(gset1.1)$Gene.symbol),],exprs(gset2.1)[match(match.gl,fData(gset2.1)$Gene.symbol),])
exprs_ <- rbind(exprs_1, exprs_2)
batch = pheno_$batch
modcombat = model.matrix(~1, data=pheno_)
combat_edata = ComBat(dat=exprs_, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)



###################################
#random stuff
###################################
d <- model.matrix(~ dvar + 0, test2.cut)
colnames(d) <- paste("G",levels(test2.cut$dvar),sep="")
fit <- lmFit(test2.cut, d)
cont.matrix <- makeContrasts(G2-G1, levels=d)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("ID","Gene.symbol","logFC"))

head(sort(table(fData(gset)$Gene.symbol), decreasing = T),20)
sum(table(fData(gset)$Gene.symbol)>2)
length(table(fData(gset)$Gene.symbol))
grep("Cdk",fData(gset)$Gene.symbol,value = T)

sum(exprs(gset)<0)

