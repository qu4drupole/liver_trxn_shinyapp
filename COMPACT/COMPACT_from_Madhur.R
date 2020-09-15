
cat("\014")
options(java.parameters = "-Xmx8000m")  #Required for xlsx package
#options(java.parameters = "- Xmx1024m")
library(limma)
library(pheatmap)
library(grid)
library(RColorBrewer)
library(xlsx)
library(magrittr)
library(vegan)
library(gplots)
library(ggfortify)
library(dplyr)

inputdir = "C:/Users/mxp212/Dropbox (SBG)/Madhur-Parihar/xenopus_retinoic_acid_madhur/Final paper sections/06 Significance Testing/"
outdir = "C:/Users/mxp212/Dropbox (SBG)/Madhur-Parihar/xenopus_retinoic_acid_madhur/Final paper sections/08 CTRL vs Treatment COMPACT/"
setwd(inputdir)

dbl = read.csv("Preprocessed_SignificanceFiltered_data.txt",stringsAsFactors=FALSE,as.is=T,row.names = 1)


##store experimental conditions as levels within a factor
batch    = factor(substr(strsplit2(colnames(dbl), "_")[, 1],1,1),levels=c("A","B","C","D","E","F"))
treatment= factor(strsplit2(colnames(dbl), "_")[, 2],levels=c("CTRL","DEAB","RA"))
time     = factor(strsplit2(colnames(dbl), "_")[, 3],levels=c("2hr","3.5hr","5hr","6.5hr"))

##I have no idea what this does, but maybe grouping the data by batch, adding new columns for new batch entries??##
sel =(batch=="") #select batch
#sel =(batch=="D"|batch=="E"|batch=="F") #select batch
if(sum(sel)==0){
  dat1 = dbl
}else{
  dat1 = dbl[,sel]  
}

batchX    = factor(substr(strsplit2(colnames(dat1), "_")[, 1],1,1)) #splits by "_" and then only takes the first element, 1,1 means start,stop
treatmentX= factor(strsplit2(colnames(dat1), "_")[, 2])
timeX     = factor(strsplit2(colnames(dat1), "_")[, 3])

### new attempt WORKING
i = interaction(treatmentX,timeX,sep="__")
# avg = apply(dat1,1,function(f) {
#   f %>% as.matrix %>% ave(i,FUN=mean) ###try with rowMeans
# })
# rownames(avg) = i
# avg2 = as.data.frame(t(avg[!duplicated(avg),]))

avg2 = sapply(split.default(dat1,i),function(x) rowMeans(x[,,drop=F])) #will probably fail on single columns

#####

ctrl = grep("CTRL", colnames(avg2),ignore.case = F)
ra   = grep("RA"  , colnames(avg2),ignore.case = F)
deab = grep("DEAB", colnames(avg2),ignore.case = F)

# avg.diff = avg2[,c(ctrl)] - avg2[,c(ra)] #RA-CTRL
# avg.diff1 = avg.diff*0
# th_hi_avg.diff = log2(1.3);    th_lo_avg.diff = -log2(1.3)
# avg.diff1[avg.diff >  th_hi_avg.diff] =  1
# avg.diff1[avg.diff <  th_lo_avg.diff] =  -1 
# avg.diff1$any = (rowSums(abs(avg.diff1))>0)
# avg2 = avg2[avg.diff1$any,]

avg3  = avg2[,c(ctrl,ra,deab)]
avg30 = avg2[,c(ctrl)]
avg3a = avg2[,c(ra)]
avg3b = avg2[,c(deab)]

##set ctrl as early timepoints - set as baseline - RA takes the place of later timepoints

c0_txt="CTRL"; cA_txt="RA"; cB_txt="DEAB"
## Generic COMPACT code starts with adjustments for timepoints
c0 = avg30[,grep(c0_txt, colnames(avg30),ignore.case = F)];rownames(c0) = rownames(avg30)
cA = avg3a[,grep(cA_txt, colnames(avg3a),ignore.case = F)];rownames(cA) = rownames(avg3a)
cB = avg3b[,grep(cB_txt, colnames(avg3b),ignore.case = F)];rownames(cB) = rownames(avg3b)


out_averaged = cbind(c0,cA,cB)
colnames(c0) = gsub(c0_txt,"0",colnames(c0))
colnames(cA) = gsub(cA_txt,"A",colnames(cA))
colnames(cB) = gsub(cB_txt,"B",colnames(cB))


##here we start to normalize the data by subtracting the baseline
## for my data, maybe use a common baseline 0 hours
c01 = cbind.data.frame(
  c0[,c(grep("_3.5hr",colnames(c0)),
        grep("_5hr"  ,colnames(c0)),
        grep("_6.5hr",colnames(c0)))] - c0[,grep("_2hr",colnames(c0))], stringsAsFactors = F)
rownames(c01) = rownames(avg30);colnames(c01) = c("0__3.5hr","0__5hr","0__6.5hr")

cA1 = cbind.data.frame(
  #A0.0 = cA[,grep("_2hr",colnames(cA))]  - c0[,grep("_2hr",colnames(c0))],
  A0.0 =0,
  cA[,c(grep("_3.5hr",colnames(cA)),
        grep("_5hr"  ,colnames(cA)),
        grep("_6.5hr",colnames(cA)))]- cA[,grep("_2hr",colnames(cA))], stringsAsFactors = F)
rownames(cA1) = rownames(avg3a);colnames(cA1) = c("A0.0","A__3.5hr","A__5hr","A__6.5hr")

cB1 = cbind.data.frame(
  #B0.0 = cB[,grep("_2hr",colnames(cB))]  - c0[,grep("_2hr",colnames(c0))],
  B0.0 = 0,
  cB[,c(grep("_3.5hr",colnames(cB)),
        grep("_5hr"  ,colnames(cB)),
        grep("_6.5hr",colnames(cB)))]- cB[,grep("_2hr",colnames(cB))], stringsAsFactors = F)
rownames(cB1) = rownames(avg3b);colnames(cB1) = c("B0.0", "B__3.5hr","B__5hr","B__6.5hr")


out_baseline_substracted = cbind.data.frame(c01,cA1[, !(colnames(cA1) %in% "A0.0")],cB1[, !(colnames(cB1) %in% "B0.0")])
out_baseline_substracted1 = out_baseline_substracted
colnames(out_baseline_substracted1)= gsub("^0_",paste0(c0_txt,"_"),colnames(out_baseline_substracted1))
colnames(out_baseline_substracted1)= gsub("^A_",paste0(cA_txt,"_"),colnames(out_baseline_substracted1))
colnames(out_baseline_substracted1)= gsub("^B_",paste0(cB_txt,"_"),colnames(out_baseline_substracted1))  #check again

#setting threshold for fold-change over the baselines
th_hi =  log2(2)
th_lo = -log2(2)
c011 = c01*0
cA11 = cA1*0
cB11 = cB1*0

#any value in the original dfs that passes the threshold in either direction is set as 1 or -1 respectively
c011[c01 >  th_hi] =  1
c011[c01 <  th_lo] =  -1 
cA11[cA1 >  th_hi] =  1
cA11[cA1 <  th_lo] =  -1 
cB11[cB1 >  th_hi] =  1
cB11[cB1 <  th_lo] =  -1 

c011$'0.all'= apply(c011[,!(colnames(c011) %in% "00.0")],1,paste,collapse="_") #No column 00.0 in c011 anyway
cA11$A.all  = apply(cA11[,!(colnames(cA11) %in% "A0.0")],1,paste,collapse="_")
cB11$B.all  = apply(cB11[,!(colnames(cB11) %in% "B0.0")],1,paste,collapse="_")

#finds all 'used' patterns that have genes
l1 = unique(c011$'0.all') ; l2 = unique(cA11$A.all); l3 = unique(cB11$B.all)
ll = unique(c(l1 , l2, l3))
ll = ll[order(ll)]

#just gives all possible patterns
ll.full = data.frame(ll.full = apply(expand.grid(c("-1","0","1"),c("-1","0","1"),c("-1","0","1")),1,paste0,collapse="_"))

patterndata = data.frame(name = rownames(c011),
                         pc01 = factor(c011$'0.all',levels=ll),
                         pcA1 = factor(cA11$A0.0),
                         pcA2 = factor(cA11$A.all,levels=ll),
                         pcB1 = factor(cB11$B0.0),
                         pcB2 = factor(cB11$B.all,levels=ll),
                         stringsAsFactors = F)

rownames(patterndata) = rownames(c011)
#generates mini table by counting #genes with each pattern?
pd_mini = lapply(names(patterndata)[c(2,4,6)], function(x){ #c(2,4,6) for "pc01","pcA2","pcB2"
  patterndata %>% group_by_(x) %>% summarise( Count= n(), genes = paste0(name,collapse = ",")) %>% merge(ll.full,by.x = x,by.y="ll.full",all= T)}) %>% as.data.frame

colnames(pd_mini) = c("Pattern.Control","Count.Control","Genes.Control",
                      "Pattern.RA"     ,"Count.RA"     ,"Genes.RA"     ,
                      "Pattern.DEAB"   ,"Count.DEAB"   ,"Genes.DEAB")
library(xlsx)
ggname = "C:/Users/mxp212/Dropbox (SBG)/Madhur-Parihar/xenopus_retinoic_acid_madhur/Final paper sections/08 CTRL vs Treatment COMPACT/Discrete Pattern Analysis1.xlsx"

##generates just gene counts for each pattern, whereas pd_mini has the gene names still
out0 = cbind(table(patterndata$pc01),  table(patterndata$pcA2),  table(patterndata$pcB2))
out1 = cbind(strsplit2(rownames(out0),"_"),out0)
#write.xlsx2(out1                  ,file=ggname, sheetName="patternwise count", row.names=T,quote=F)
#write.xlsx2(as.matrix(patterndata),file=ggname, sheetName="genewise pattern" , row.names=T,quote=F, append = T)


setwd(outdir)
## CTRL vs RA COMPACT
pd_0vsA = aggregate(name~pcA1+pc01+pcA2, data = patterndata,drop=T,na.action= NULL,FUN = function(z) {
  list1 = paste0(z,sep=",",collapse="");count=NROW(z);d = paste0(count,"---",list1)
})

COMPACT(pd_0vsA,c0_txt,cA_txt,"filename_0vsA.xlsx")


## CTRL vs DEAB COMPACT
pd_0vsB = aggregate(name~pcB1+pc01+pcB2, data = patterndata,drop=T,na.action= NULL,FUN = function(z) {
  list1 = paste0(z,sep=",",collapse="");  count=NROW(z);d = paste0(count,"---",list1)
})

COMPACT(pd_0vsB,c0_txt,cB_txt,"filename_0vsB.xlsx")


################
excel_ConvertToNumbers <- function(hhname) {
  ############################################
  ######Macro to convert text to numeric #####
  library(RDCOMClient)
  template.excel = "C:/Users/mxp212/Dropbox (SBG)/Madhur-Parihar/xenopus_retinoic_acid_madhur/Final paper sections/08 CTRL vs Treatment COMPACT/Macro template.xlsm"  
  target.excel   = hhname  
  # Open a specific workbook in Excel:
  excel1 <- COMCreate("Excel.Application")
  wb1 <- excel1$Workbooks()$Open(target.excel)
  excel2 <- COMCreate("Excel.Application")
  wb2 <- excel2$Workbooks()$Open(template.excel)
  # this line of code might be necessary if you want to see your spreadsheet:
  excel2[['Visible']] <- TRUE
  
  # Run the macro called "MyMacro":
  excel1$Run("ConvertTextToNumbers")
  
  # # Close the workbook and quit the app:
  # wb1$Save(Filename = template.excel)
  # wb2$Save(Filename = target.excel)
  # xlWbk$SaveAs(gsub("/", "\\\\", file.path(getwd(), "test.xls")))
  # wb1$Close(TRUE)
  # wb2$Close(TRUE)
  # excel1$Quit()
  
  # Release resources:
  rm(wb1, excel1,wb2, excel2)
  gc()
}

##########
COMPACT <- function(pd, row_txt,column_txt,ffname) {
  #browser()
  patterndata = patterndata[colnames(patterndata) %in% colnames(pd)]
  g = quads(patterndata)
  temp = strsplit2(pd$name,"---")
  pd$genecounts = as.numeric(temp[,1])
  pd$genelist = temp[,2]
  pd$genelist[nchar(pd$genelist)>32000] = substr(pd$genelist[nchar(pd$genelist)>32000],1,31000)
  pd$name = NULL
  pd = pd[order(pd$genecounts,decreasing = T),]
  #pd$serial = 1:dim(pd)[1]
  pd$genelist [is.na(pd$genelist)] = ""
  pd$genecounts [is.na(pd$genecounts)] = 0
  #x = table(patterndata[,c(2,3,4)])
  # x    = array(data = pd$genecounts[order(pd$serial)],
  #              dim = c(length(levels(pd_0vsA[,1])),length(levels(pd_0vsA[,2])),length(levels(pd_0vsA[,3]))),
  #              dimnames = list(levels(pd_0vsA[,1]),levels(pd_0vsA[,2]),levels(pd_0vsA[,3]))
  #              )
  # xxx  = array(data = pd$genelist[order(pd$serial)]  ,dim = dim(x))  # 3d matrix of gene names
  # 
  
  #x = table(patterndata$pconditionB1,patterndata$pconditionA1,patterndata$pconditionB2,deparse.level=2)
  x = table(patterndata[,c(3,2,4)],deparse.level=2) #Order: pA1,p01,pA2
  x[x==max(x)] = 0
  b = as.data.frame(x)
  b$serial = as.numeric(rownames(b))
  colnames(b)[1:3]= colnames(pd)[1:3]
  mrg = base::merge(b,pd,all.x=T)
  mrg$genelist [is.na(mrg$genelist)] = ""
  mrg$genecounts [is.na(mrg$genecounts)] = 0
  xxx  = array(data = mrg$genelist[order(mrg$serial)],dim = dim(x))  # 3d matrix of gene names
  
  
  k = rep(unlist(dimnames(x)[1]),each = dim(x)[2])
  #k = rep(unlist(pd$pc01), each = dim(x)[2])
  df = aperm(x,c(2,3,1))
  dim(df) = c(dim(x)[2],dim(x)[1]*dim(x)[3])  #flatten out 3d array to 3 column bind frames  ##alternate option arrayhelpers::array2df
  
  colnames(df) = paste0(column_txt,"_",k,"_",unlist(dimnames(x)[2]))
  rownames(df) = paste0(row_txt,"_",unlist(dimnames(x)[2]))
  df = addmargins(df)
  df[df==0]=""
  r1 = do.call("rbind",strsplit(rownames(as.data.frame.matrix(df)),"_"))
  c1 = do.call("rbind",strsplit(colnames(as.data.frame.matrix(df)),"_"))
  
  df1 = as.data.frame(cbind(matrix(data = "-",ncol= dim(r1)[2],nrow = dim(c1)[2]),t(c1)),stringsAsFactors = F)
  df2 = as.data.frame(cbind(r1,as.data.frame.matrix(df)),stringsAsFactors = F)
  names(df1) = names(df2)
  
  op = rbind(df1,df2)
  
  dfz = aperm(xxx,c(2,3,1))  # COMPACT of gene names
  dim(dfz) = c(dim(xxx)[2],dim(xxx)[1]*dim(xxx)[3])  #flatten out 3d array to 3 column bind frames  ##alternate option arrayhelpers::array2df
  
  r1 = do.call("rbind",strsplit(rownames(as.data.frame.matrix(df)),"_"))
  c1 = do.call("rbind",strsplit(colnames(as.data.frame.matrix(df)),"_"))
  r1 = r1[-nrow(r1),]
  c1 = c1[-nrow(c1),]
  df1 = as.data.frame(cbind(matrix(data = "-",ncol= dim(r1)[2],nrow = dim(c1)[2]),t(c1)),stringsAsFactors = F)
  df2 = as.data.frame(cbind(r1,as.data.frame.matrix(dfz)),stringsAsFactors = F)
  names(df1) = names(df2)
  
  opz = rbind(df1,df2)
  
  
  
  ##Save files using xlsx package
  library(xlsx)
  #ffname = "filename.xlsx"
  
  out_for_heatmaps = merge(patterndata,round(out_baseline_substracted1,digits = 2),by=0)[,-1]
  #out_for_heatmaps$pattern = paste0(column_txt,out_for_heatmaps$pcA1,"__",row_txt,"_",out_for_heatmaps$pcB1,"_",out_for_heatmaps$pcB2)
  
  print("Writing file 1/7")
  write.xlsx2(op,                     file=ffname, sheetName="count compact", row.names=F,quote=F)
  print("Writing file 2/7")
  write.xlsx2(opz,                    file=ffname, sheetName="name compact", row.names=F,append=T)
  print("Writing file 3/7")
  write.xlsx2(patterndata,            file=ffname, sheetName="genewise pattern", row.names=F,append=T)
  print("Writing file 4/7")
  write.xlsx2(pd,                     file=ffname, sheetName="patternwise genes", row.names=F,append=T)
  print("Writing file 5/7")
  write.xlsx2(round(out_averaged,2),  file=ffname, sheetName="averaged data", row.names=T,append=T)
  print("Writing file 6/7")
  write.xlsx2(out_for_heatmaps,       file=ffname, sheetName="baseline substracted data", row.names=T,append=T)
  print("Writing file 7/7")
  # write.xlsx2(round(avg3,2),        file=ffname, sheetName="all treatments averaged data", row.names=T,append=T)
  print("Done...")
  
  ###  Running a VBA macro to convert text into numeric in excel
  #excel_ConvertToNumbers(ffname)
}