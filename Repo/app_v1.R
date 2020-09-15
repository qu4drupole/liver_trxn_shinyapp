library(shiny)
library(GEOquery)
library(limma)
library(gplots)

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
                                 "GSE97429" = "GSE97429.RData"))),
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
                                 "GSE97429" = "GSE97429.RData"))),
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
      tableOutput("test"),
      plotOutput("heatmap")
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
                        choices = c("null",levels(rv$d1$Treatment)))
      updateSelectInput(session, "time1a",
                        choices = c("null",levels(rv$d1$Time)))
      updateSelectInput(session, "trmt1b",
                        choices = c("null",levels(rv$d1$Treatment)))
      updateSelectInput(session, "time1b",
                        choices = c("null",levels(rv$d1$Time)))
    }
  })
  
  observeEvent(input$set1,{
    if(input$trmt1a != "null" & input$time1a != "null" & input$trmt1b != "null" & input$time1b != "null"){
      c1 <- which(rv$d1$Treatment == input$trmt1a & rv$d1$Time == input$time1a)
      c2 <- which(rv$d1$Treatment == input$trmt1b & rv$d1$Time == input$time1b)
      rv$cut1 <- rv$d1[, c(c1,c2)]
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
      rv$tT <- subset(rv$tT, select=c("ID","Gene.symbol","GenBank.Accession","logFC"))
      #output$test <- renderTable({head(rv$tT)})
    }
    })
  
  observeEvent(input$data2, {
    if(input$data2 != "null"){
      load(paste("Data/",input$data2,sep=""))
      rv$d2 <- gset
      output$Title2 <- renderText({paste("Title:",experimentData(rv$d2)@title)})
      updateSelectInput(session, "trmt2a",
                        choices = c("null",levels(rv$d2$Treatment)))
      updateSelectInput(session, "time2a",
                        choices = c("null",levels(rv$d2$Time)))
      updateSelectInput(session, "trmt2b",
                        choices = c("null",levels(rv$d2$Treatment)))
      updateSelectInput(session, "time2b",
                        choices = c("null",levels(rv$d2$Time)))
    }
  })
  
  observeEvent(input$set2,{
    if(input$trmt2a != "null" & input$time2a != "null" & input$trmt2b != "null" & input$time2b != "null"){
      c1 <- which(rv$d2$Treatment == input$trmt2a & rv$d2$Time == input$time2a)
      c2 <- which(rv$d2$Treatment == input$trmt2b & rv$d2$Time == input$time2b)
      rv$cut2 <- rv$d2[, c(c1,c2)]
      rv$cut2$dvar <- factor(match(paste(rv$cut2$Treatment,rv$cut2$Time),unique(paste(rv$cut2$Treatment,rv$cut2$Time))))
      #this assumes every data set will have a "geo_accession" slot
      fData(rv$cut2)$logFC <- as.numeric(apply(exprs(rv$cut2),1,function(x) mean(x[rv$cut2$geo_accession[rv$cut2$dvar==2]]) - mean(x[rv$cut2$geo_accession[rv$cut2$dvar==1]])))
      #output$test <- renderTable({head(subset(fData(rv$cut2), select = c("ID","logFC")))})
    }
  })
  
  observeEvent(input$GO,{
    fData(rv$cut2)$Gene.symbol[fData(rv$cut2)$Gene.symbol == ""] <- NA
    rv$tT$logFC.d2 <- fData(rv$cut2)$logFC[match(rv$tT$Gene.symbol,fData(rv$cut2)$Gene.symbol)]
    tT.trim <- subset(rv$tT, select = c("Gene.symbol","logFC","logFC.d2"))
    output$heatmap <- renderPlot({heatmap.2(as.matrix(tT.trim[,-1]),
                                            trace = 'none', 
                                            dendrogram = 'row', 
                                            labCol = c("Dataset 1","Dataset 2"),
                                            srtCol = 0, 
                                            adjCol = c(0.5,0.6), 
                                            cexCol = 3, 
                                            labRow = NA)}, height = 800)
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