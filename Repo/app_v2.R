library(shiny)
library(GEOquery)
library(limma)
library(gplots)
library(FNN)
library(Biobase)
library(RColorBrewer)
library(factoextra)

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
        hr(),
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
        actionButton("GO", "Analyze!")
      ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Main", 
                           textOutput("divergence"),
                           tableOutput("test"),
                           plotOutput("heatmap")),
                  tabPanel("MA plots", 
                           fluidRow(column(3,tableOutput("samples1")),
                                    column(1),
                                    column(3,tableOutput("samples2"))),
                           plotOutput("MA")),
                  tabPanel("Clustering", 
                           fluidRow(column(3,plotOutput("elbows")), 
                                    column(2,selectInput("nclusters", "Number of clusters:", c("null" = "null",1:10)))),
                           plotOutput("clusts"),
                           tableOutput("clusters"),
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
    cut2 = 1
  )
  
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
      c1 <- which(rv$d1$Treatment == input$trmt1a & rv$d1$Time == input$time1a)
      c2 <- which(rv$d1$Treatment == input$trmt1b & rv$d1$Time == input$time1b)
      rv$cut1 <- rv$d1[, c(c1,c2)]
      if(any(table(rv$cut1$ID>1))){
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
      rv$tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
      #Need to include error handling for available feature data
      #....maybe make the fvarLablels a checklist
      rv$tT <- subset(rv$tT, select=c("ID","Gene.symbol","logFC"))
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
      if(any(table(rv$cut2$ID>1))){
        c3 <- which(rv$cut2$ID %in% names(which(table(rv$cut2$ID)>1)))
        rv$cut2 <- rv$cut2[,c3]
      }
      rv$cut2$dvar <- factor(match(paste(rv$cut2$Treatment,rv$cut2$Time),unique(paste(rv$cut2$Treatment,rv$cut2$Time))))
    }
  })
  
  observeEvent(input$GO,{
    sig1 <- exprs(rv$cut1)[match(rv$tT$Gene.symbol,fData(rv$cut1)$Gene.symbol),]
    sig2 <- exprs(rv$cut2)[match(rv$tT$Gene.symbol,fData(rv$cut2)$Gene.symbol),]
    hm.data <- cbind(sig1, sig2)
    a <- apply(hm.data,1,function(x) sum(!is.na(x)) == length(x))
    hm.data <- hm.data[a,]
    d1lab <- sapply(rownames(pData(rv$cut1)), function(x) paste("d1",pData(rv$cut1)[x,"Treatment"],pData(rv$cut1)[x,"Time"],sep = "."),USE.NAMES = F)
    d2lab <- sapply(rownames(pData(rv$cut2)), function(x) paste("d2",pData(rv$cut2)[x,"Treatment"],pData(rv$cut2)[x,"Time"],sep = "."),USE.NAMES = F)
    # invisible(p <- heatmap.2(hm.data,
    #                trace = 'none',
    #                Colv = NULL,
    #                dendrogram = 'row',
    #                srtCol = 90,
    #                adjCol = c(1,0.6),
    #                labCol = c(d1lab,d2lab),
    #                cexCol = 1.8,
    #                labRow = NA,
    #                margins = c(20,1),
    #                keysize = 0.75))
    kl.data <- hm.data[p$rowInd,]
    a <- rowMeans(kl.data[,rownames(pData(rv$cut1))[rv$cut1$dvar==2]])-rowMeans(kl.data[,rownames(pData(rv$cut1))[rv$cut1$dvar==1]])
    b <- rowMeans(kl.data[,rownames(pData(rv$cut2))[rv$cut2$dvar==2]])-rowMeans(kl.data[,rownames(pData(rv$cut2))[rv$cut2$dvar==1]])
    KL.entropy <- KL.divergence(a,b,50)
    KL.index <- which(KL.entropy != Inf)[1]
    KL.value <- round(KL.entropy[KL.index],4)
    output$heatmap <- renderPlot({heatmap.2(hm.data,
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
    output$divergence <- renderText({paste("KL divergence =",KL.value,"with NN =",KL.index)})
    #making the MA tab
    output$samples1 <- renderTable({subset(pData(rv$cut1), select = c("Treatment","Time", "geo_accession"))})
    output$samples2 <- renderTable({subset(pData(rv$cut2), select = c("Treatment","Time", "geo_accession"))})
    output$MA <- renderPlot({
      cut2.match <- exprs(rv$cut2)[match(fData(rv$cut1)$Gene.symbol,fData(rv$cut2)$Gene.symbol),]
      ma.data <- cbind(exprs(rv$cut1),cut2.match)
      par(mfrow=c(ceiling(ncol(ma.data)/3),3), mar = c(1,1,1,1))
      for (i in 1:ncol(ma.data)){
        plotMA(ma.data,i)
      }}, height = 800)
    output$elbows <- renderPlot({
      a <- rownames(pData(rv$cut1)[rv$cut1$dvar == 1])
      b <- rownames(pData(rv$cut1)[rv$cut1$dvar == 2])
      c <- rownames(pData(rv$cut2)[rv$cut1$dvar == 1])
      d <- rownames(pData(rv$cut2)[rv$cut1$dvar == 2])
      cl.data <- cbind(rowMedians(hm.data[,a]),rowMedians(hm.data[,b]),rowMedians(hm.data[,c]),rowMedians(hm.data[,d]))
      plot(fviz_nbclust(cl.data, hcut, method = "wss"))
    }, height = 300)
    output$clusts <- renderPlot({
      if(input$nclusters != "null"){
        cl.fit <- hclust(dist(cl.data))
        cl.test <- cutree(cl.fit, k=input$nclusters)
        selcol <- brewer.pal(input$nclusters,"Set1")
        heatmap.2(cl.data[cl.fit$order,], 
                  trace = 'none', 
                  Colv = NULL,
                  dendrogram = 'row',
                  srtCol = 90, 
                  adjCol = c(0.8,0.6), 
                  cexCol = 1, 
                  labRow = NA,
                  margins = c(7,1),
                  key = F,
                  RowSideColors = selcol[cl.test[cl.fit$order]])}
    }, height = 600)
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