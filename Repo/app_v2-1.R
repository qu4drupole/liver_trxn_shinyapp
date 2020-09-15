library(shiny)
library(GEOquery)
library(limma)
library(gplots)
library(FNN)
library(Biobase)
library(RColorBrewer)
library(factoextra)
library(sva)

#update version 1.1
#changed heatmap from averages to samples
#adjustment to account for paired samples

ui <- fluidPage(
  titlePanel("Liver Regeneration Meta-Analysis"),
  sidebarLayout(
    sidebarPanel(
      tags$head(tags$style(type="text/css", "
             #loadmessage {
                           position: fixed;
                           top: 0px;
                           left: 0px;
                           width: 100%;
                           padding: 5px 0px 5px 0px;
                           text-align: center;
                           font-weight: bold;
                           font-size: 100%;
                           color: #000000;
                           background-color: #f3f936;
                           z-index: 105;
                           }
                           ")),
      radioButtons("aType","Analysis type",
                   choices = c("Discover" = "discover",
                               "Use gene list" = "geneList")),
      conditionalPanel(condition = "input.aType == 'geneList'", 
                       textAreaInput("inGenes","Paste gene list:", height = 100)),
      ###DATASET 1####
      h4("Dataset 1"),
      fluidRow(
        column(4,selectInput("data1","Choose data:",
                             c("null" = "null",
                               "GSE63742" = "GSE63742.RData",
                               "GSE33785" = "GSE33785.RData",
                               "GSE97429" = "GSE97429.RData",
                               "GSE67022" = "GSE67022.RData"))),
        column(2, br(),
               downloadButton("downloadData1","Download"))
      ),
      p(textOutput("Title1")),
      conditionalPanel(condition = "input.aType == 'discover'", 
                       numericInput("numGenes", "Number of genes in top table:", 250, width = 100)),
      fluidRow(
        em("CONTRASTS:"),
        br(),
        column(4,
               selectInput("trmt1a","Treatment:",c("null" = "null")),
               em("versus"),
               selectInput("trmt1b","Treatment:",c("null" = "null"))
        ),
        column(4,
               selectInput("time1a","Time:",c("null" = "null")),
               br(),
               selectInput("time1b","Time:",c("null" = "null"))),
        column(2,
               br(),br(),br(),br(),br(),br(),
               actionButton("set1","Prepare Dataset 1"))
      ),
      ###DATASET 2####
      h4("Dataset 2"),
      fluidRow(
        column(4,selectInput("data2","Choose data:",
                             c("null" = "null",
                               "GSE63742" = "GSE63742.RData",
                               "GSE33785" = "GSE33785.RData",
                               "GSE97429" = "GSE97429.RData",
                               "GSE67022" = "GSE67022.Rdata"))),
        column(2, br(),
               downloadButton("downloadData2","Download"))
      ),
      p(textOutput("Title2")),
      fluidRow(
        em("CONTRASTS:"),
        br(),
        column(4,
               selectInput("trmt2a","Treatment:",c("null" = "null")),
               em("versus"),
               selectInput("trmt2b","Treatment:",c("null" = "null"))
        ),
        column(4,
               selectInput("time2a","Time:",c("null" = "null")),
               br(),
               selectInput("time2b","Time:",c("null" = "null"))),
        column(2,
               br(),br(),br(),br(),br(),br(),
               actionButton("set2","Prepare Dataset 2"))
      ),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Working...",id="loadmessage")),
      conditionalPanel(condition = "input.aType == 'discover'", actionButton("GO", "Discover!")),
      conditionalPanel(condition = "input.aType == 'geneList'", 
                       actionButton("GO.gl", "Analyze!"))
      ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Main", 
                           textOutput("divergence"),
                           textOutput("test"),
                           plotOutput("heatmap")),
                  tabPanel("MA plots", 
                           fluidRow(column(3,tableOutput("samples1")),
                                    column(1),
                                    column(3,tableOutput("samples2"))),
                           plotOutput("MA")),
                  tabPanel("Clustering", 
                           fluidRow(column(3,plotOutput("elbows",height = 300)), 
                                    column(2,selectInput("nclusters", "Number of clusters:", c("null" = "null", 1:10)))),
                           #textOutput("test2"),
                           plotOutput("clusterKey",height = 100),
                           plotOutput("clusts", height = 600),
                           fluidRow(column(2,selectInput("clusterN", "Select cluster",
                                                         c("null" = "null"))),
                                    column(2,tableOutput("geneList")))),
                  tabPanel("Help", plotOutput("help")))
    )
  ))

server <- function(input, output, session) {
  rv <- reactiveValues(
    d1 = 1,
    cut1 = 1,
    tT = 1,
    d2 = 1,
    cut2 = 1,
    eData = 1
  )
  
  # output$test <- renderText({
  #   strsplit(input$inGenes, "\n")[[1]][3]
  # })
  
  #fp <- reactive({
  #  paste("Data/mRNA/",input$data1,sep="")})
  
  observeEvent(input$data1, {
    if(input$data1 != "null"){
      load(paste("Data/",input$data1,sep=""))
      rv$d1 <- gset
      output$Title1 <- renderText({paste("Title:",experimentData(rv$d1)@title)})
      updateSelectInput(session, "trmt1a",
                        choices = c("null", levels(rv$d1$Treatment)))
      updateSelectInput(session, "trmt1b",
                        choices = c("null",levels(rv$d1$Treatment)))
    }
  })
  
  observeEvent(input$trmt1a, {
    if(input$trmt1a != "null"){
      updateSelectInput(session, "time1a",
                        choices = c("null",
                                    levels(factor((rv$d1$Time[rv$d1$Treatment == input$trmt1a])))))}
  })
  
  observeEvent(input$trmt1b, {
    if(input$trmt1b != "null"){
      updateSelectInput(session, "time1b",
                        choices = c("null",
                                    levels(factor((rv$d1$Time[rv$d1$Treatment == input$trmt1b])))))}
  })
  
  observeEvent(input$set1,{
    if(input$trmt1a != "null" & input$time1a != "null" & input$trmt1b != "null" & input$time1b != "null"){
      if(input$aType == "discover"){
        c1 <- which(rv$d1$Treatment == input$trmt1a & rv$d1$Time == input$time1a)
        c2 <- which(rv$d1$Treatment == input$trmt1b & rv$d1$Time == input$time1b)
        rv$cut1 <- rv$d1[, c(c1,c2)]
        if(any(table(rv$cut1$ID) > 1)){
          c3 <- which(rv$cut1$ID %in% names(which(table(rv$cut1$ID)>1)))
          rv$cut1 <- rv$cut1[,c3]
        }
        rv$cut1$dvar <- factor(match(paste(rv$cut1$Treatment,rv$cut1$Time),unique(paste(rv$cut1$Treatment,rv$cut1$Time))))
        design <- model.matrix(~ dvar + 0, rv$cut1)
        colnames(design) <- paste("G",levels(rv$cut1$dvar),sep="")
        fit <- lmFit(rv$cut1, design)
        cont.matrix <- makeContrasts(G2-G1, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        rv$tT <- topTable(fit2, adjust="fdr", sort.by="B", number=input$numGenes)
        #Need to include error handling for available feature data
        #....maybe make the fvarLablels a checklist
        rv$tT <- subset(rv$tT, select=c("ID","Gene.symbol","logFC"))
        rv$tT <- rv$tT[rv$tT$Gene.symbol != "",]
        rv$tT$ID <- as.character(rv$tT$ID)}
      #This just skips the differential expression step
      if(input$aType == "geneList"){
        c1 <- which(rv$d1$Treatment == input$trmt1a & rv$d1$Time == input$time1a)
        c2 <- which(rv$d1$Treatment == input$trmt1b & rv$d1$Time == input$time1b)
        rv$cut1 <- rv$d1[, c(c1,c2)]
        if(any(table(rv$cut1$ID) >1 )){
          c3 <- which(rv$cut1$ID %in% names(which(table(rv$cut1$ID)>1)))
          rv$cut1 <- rv$cut1[,c3]
        }
        rv$cut1$dvar <- factor(match(paste(rv$cut1$Treatment,rv$cut1$Time),unique(paste(rv$cut1$Treatment,rv$cut1$Time))))
      }
    }
    })
  
  observeEvent(input$data2, {
    if(input$data2 != "null"){
      load(paste("Data/",input$data2,sep=""))
      rv$d2 <- gset
      output$Title2 <- renderText({paste("Title:",experimentData(rv$d2)@title)})
      updateSelectInput(session, "trmt2a",
                        choices = c("null",levels(rv$d2$Treatment)))
      updateSelectInput(session, "trmt2b",
                        choices = c("null",levels(rv$d2$Treatment)))
    }
  })
  
  observeEvent(input$trmt2a, {
    if(input$trmt2a != "null"){
      updateSelectInput(session, "time2a",
                        choices = c("null",
                                    levels(factor((rv$d2$Time[rv$d2$Treatment == input$trmt2a])))))}
  })
  
  observeEvent(input$trmt2b, {
    if(input$trmt2b != "null"){
      updateSelectInput(session, "time2b",
                        choices = c("null",
                                    levels(factor((rv$d2$Time[rv$d2$Treatment == input$trmt2b])))))}
  })
  
  observeEvent(input$set2,{
    if(input$trmt2a != "null" & input$time2a != "null" & input$trmt2b != "null" & input$time2b != "null"){
      c1 <- which(rv$d2$Treatment == input$trmt2a & rv$d2$Time == input$time2a)
      c2 <- which(rv$d2$Treatment == input$trmt2b & rv$d2$Time == input$time2b)
      rv$cut2 <- rv$d2[, c(c1,c2)]
      if(any(table(rv$cut2$ID) > 1)){
        c3 <- which(rv$cut2$ID %in% names(which(table(rv$cut2$ID)>1)))
        rv$cut2 <- rv$cut2[,c3]
      }
      rv$cut2$dvar <- factor(match(paste(rv$cut2$Treatment,rv$cut2$Time),unique(paste(rv$cut2$Treatment,rv$cut2$Time))))
    }
  })
  
  cv <- reactiveValues(
    hm.data = 1,
    cl.data = 1,
    cl.fit = 1,
    cl.test =1
  )
  
  observeEvent(input$GO,{
#ComBat normalization
    gset1.1 <- rv$cut1[fData(rv$cut1)$Gene.symbol != "",]
    gset2.1 <- rv$cut2[fData(rv$cut2)$Gene.symbol != "",]
    pData(gset1.1)$batch <- 1
    pData(gset2.1)$batch <- 2
    pheno_ <- rbind(subset(pData(gset1.1),select=c("Treatment", "Time", "dvar", "batch")), subset(pData(gset2.1),select=c("Treatment","Time", "dvar", "batch")))
    match.gl <- fData(gset2.1)$Gene.symbol[match(unique(fData(gset1.1)$Gene.symbol),fData(gset2.1)$Gene.symbol)]
    match.gl <- match.gl[!is.na(match.gl)]
    exprs_1 <- cbind(exprs(gset1.1)[rv$tT$ID,],exprs(gset2.1)[match(rv$tT$Gene.symbol,fData(gset2.1)$Gene.symbol),])
    a <- apply(exprs_1, 1, function(x) sum(!is.na(x)) == length(x))
    exprs_1 <- exprs_1[a,]
    match.gl <- match.gl[which(!(match.gl %in% rv$tT$Gene.symbol))]
    exprs_2 <- cbind(exprs(gset1.1)[match(match.gl,fData(gset1.1)$Gene.symbol),],exprs(gset2.1)[match(match.gl,fData(gset2.1)$Gene.symbol),])
    exprs_ <- rbind(exprs_1, exprs_2)
    batch = pheno_$batch
    modcombat = model.matrix(~1, data=pheno_)
    rv$eData = ComBat(dat=exprs_, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
#Making the main heatmap
    cv$hm.data <- exprs_1
    d1lab <- sapply(rownames(pData(rv$cut1)), function(x) paste("d1",pData(rv$cut1)[x,"Treatment"],pData(rv$cut1)[x,"Time"],sep = "."),USE.NAMES = F)
    d2lab <- sapply(rownames(pData(rv$cut2)), function(x) paste("d2",pData(rv$cut2)[x,"Treatment"],pData(rv$cut2)[x,"Time"],sep = "."),USE.NAMES = F)
    output$heatmap <- renderPlot({heatmap.2(cv$hm.data,
                                             trace = 'none',
                                             Colv = NULL,
                                             dendrogram = 'row',
                                             srtCol = 90,
                                             adjCol = c(1,0.6),
                                             labCol = c(d1lab,d2lab),
                                             cexCol = 1.8,
                                             labRow = NA,
                                             margins = c(20,1),
                                             keysize = 0.75)}, height = 1000)
    output$divergence <- renderText({
      kl.data <- cv$hm.data
      a <- rowMeans(kl.data[,rownames(pData(rv$cut1))[rv$cut1$dvar==2]])-rowMeans(kl.data[,rownames(pData(rv$cut1))[rv$cut1$dvar==1]])
      b <- rowMeans(kl.data[,rownames(pData(rv$cut2))[rv$cut2$dvar==2]])-rowMeans(kl.data[,rownames(pData(rv$cut2))[rv$cut2$dvar==1]])
      KL.entropy <- KL.divergence(a,b,50)
      KL.index <- which(KL.entropy != Inf)[1]
      KL.value <- round(KL.entropy[KL.index],4)
      paste("KL divergence =",KL.value,"with NN =",KL.index)})
#making the MA tab
    output$samples1 <- renderTable({subset(pData(rv$cut1), select = c("Treatment","Time", "geo_accession"))})
    output$samples2 <- renderTable({subset(pData(rv$cut2), select = c("Treatment","Time", "geo_accession"))})
    output$MA <- renderPlot({
      par(mfrow=c(ceiling(ncol(rv$eData)/3),3), mar = c(1,1,1,1))
      for (i in 1:ncol(rv$eData)){
        plotMA(rv$eData,i)
      }}, height = 800)
#making the clustering tab
    output$elbows <- renderPlot({
      a <- rownames(pData(rv$cut1))[rv$cut1$dvar == 1]
      b <- rownames(pData(rv$cut1))[rv$cut1$dvar == 2]
      c <- rownames(pData(rv$cut2))[rv$cut2$dvar == 1]
      d <- rownames(pData(rv$cut2))[rv$cut2$dvar == 2]
      cv$cl.data <- cbind(rowMedians(cv$hm.data[,a]),rowMedians(cv$hm.data[,b]),rowMedians(cv$hm.data[,c]),rowMedians(cv$hm.data[,d]))
      plot(fviz_nbclust(cv$cl.data, hcut, method = "wss"))
    }, height = 300)
    observeEvent(input$nclusters,{
      if(input$nclusters != "null"){
        updateSelectInput(session, "clusterN",
                          choices = c("null",1:as.numeric(input$nclusters)))
        cv$cl.fit <- hclust(dist(cv$cl.data))
        cv$cl.test <- cutree(cv$cl.fit, k=as.numeric(input$nclusters))
        selcol <- brewer.pal(as.numeric(input$nclusters),"Set1")
        output$clusts <- renderPlot({
          heatmap.2(cv$cl.data[cv$cl.fit$order,], 
                    trace = 'none', 
                    Colv = NULL,
                    dendrogram = 'row',
                    srtCol = 90, 
                    adjCol = c(0.8,0.6), 
                    cexCol = 1, 
                    labRow = NA,
                    lmat = rbind(c(0,0,4),c(3,1,2),c(5,0,0)),
                    lhei = c(0.2,4,0.2),
                    lwid = c(1,0.2,4),
                    margins = c(2,20),
                    key = F,
                    RowSideColors = selcol[cv$cl.test[cv$cl.fit$order]])
          }, height = 600)
        output$clusterKey <- renderPlot({
          ggplot(data.frame(x=1:as.numeric(input$nclusters),y=rep(1,as.numeric(input$nclusters))),aes(x=x,y=y))+
            geom_bar(fill = selcol, stat = "identity", show.legend = F)+
            coord_equal()+
            xlab("Clusters")+
            geom_text(aes(label=1:as.numeric(input$nclusters)), size = 5, nudge_y = -0.5)+
            theme(panel.background = element_rect(fill = "white"),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.x = element_text(size=20))}, height = 100)
      }
    })  
    output$geneList <- renderTable({
      if(input$clusterN != "null"){
        probeID <- rownames(cv$hm.data)[cv$cl.fit$order[cv$cl.test==as.numeric(input$clusterN)]]
        fData(rv$cut1)[probeID,"Gene.symbol"]
      }
    })
  })
  
#running the gene list analysis
  observeEvent(input$GO.gl,{
    sig1 <- exprs(rv$cut1)[match(strsplit(input$inGenes, "\n")[[1]],fData(rv$cut1)$Gene.symbol),]
    sig2 <- exprs(rv$cut2)[match(strsplit(input$inGenes, "\n")[[1]],fData(rv$cut2)$Gene.symbol),]
    cv$hm.data <- cbind(sig1, sig2)
    a <- apply(cv$hm.data,1,function(x) sum(!is.na(x)) == length(x))
    cv$hm.data <- cv$hm.data[a,]
    d1lab <- sapply(rownames(pData(rv$cut1)), function(x) paste("d1",pData(rv$cut1)[x,"Treatment"],pData(rv$cut1)[x,"Time"],sep = "."),USE.NAMES = F)
    d2lab <- sapply(rownames(pData(rv$cut2)), function(x) paste("d2",pData(rv$cut2)[x,"Treatment"],pData(rv$cut2)[x,"Time"],sep = "."),USE.NAMES = F)
    output$heatmap <- renderPlot({heatmap.2(cv$hm.data,
                                            trace = 'none',
                                            Colv = NULL,
                                            dendrogram = 'row',
                                            srtCol = 90,
                                            adjCol = c(1,0.6),
                                            labCol = c(d1lab,d2lab),
                                            cexCol = 1.8,
                                            labRow = NA,
                                            margins = c(20,1),
                                            keysize = 0.75)}, height = 1000)
  })
  
  #still need to finish this
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file){
      write.csv(isolate(rv$d1),file)
    })
}

shinyApp(ui = ui, server = server)