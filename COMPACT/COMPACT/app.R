library(shiny)
library(limma)
library(Biobase)
library(dplyr)
library(tidyverse)
library(ggplot2)


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
      selectInput("data","Choose data:",
                  c("null" = "null",
                    "GSE20426" = "GSE20426.RData",
                    "GSE100616" = "GSE100616.RData",
                    "GSE33785" = "GSE33785.RData")),
      p(textOutput("Title")),
      fluidRow(
        column(5,selectInput("Var1","Condition 1:",c("null"="null"))),
        column(5,selectInput("Var2","Condition 2:",c("null"="null")))),
      selectInput("baseline", "Choose baseline:",
                  c("null" = "null")),
      sliderInput("threshold", "fold-change threshold",
                  min = 1.2, max = 5, value = 2),
      numericInput("minSet", "minimum pattern set size:", value = 15),
      actionButton("GO", "Construct COMPACT matrix"),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Working...",id="loadmessage"))
    ),
    mainPanel(
      tableOutput("test"),
      plotOutput("COMPACT")
    )
  )
)

server <- function(input, output, session) {
  eDcollapse <- function(g){
    #browser()
    #count unique trmt*time combinations
    cols = 0
    for(i in levels(g$Time)){
      cols = cols + length(unique(g$Trmt[g$Time==i]))
    }
    res.avg <- data.frame(matrix(nrow = nrow(fData(g)), ncol = cols))
    ci = 1
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
  
  rv <- reactiveValues(
    dat = 1,
    groups = 1,
    contrast = 1,
    discrete = 1,
    cm = 1
  )
  
  observeEvent(input$data, {
    if(input$data != "null"){
      load(paste("Data/",input$data,sep=""))
      rv$dat <- gset
      output$Title <- renderText({paste("Title:",experimentData(rv$dat)@title)})
      updateSelectInput(session, "Var1",
                        choices = c("null", levels(rv$dat$Trmt)))
      updateSelectInput(session, "Var2",
                        choices = c("null", levels(rv$dat$Trmt)))
      updateSelectInput(session, "baseline",
                        choices = c("null", levels(rv$dat$Time)))
    }
  })
  
  observeEvent(input$GO, {
    rv$groups <- eDcollapse(rv$dat)
    rv$contrast <- eDcontrast(rv$groups, as.numeric(input$baseline))
    rv$discrete <- eDdiscretize(rv$contrast, input$minSet, input$threshold)
    rv$cm <- makeCOMPACT(rv$discrete, input$Var1, input$Var2)

    #making the heatmap
    hm.data <- rv$cm %>% rownames_to_column('t2') %>% gather(rownames(rv$cm), value, -t2)
    #NEED TO DOUBLE CHECK THAT NAMING!
    a <- gsub(" ","_",input$Var1)
    b <- gsub(" ","_",input$Var2)
    names(hm.data) <- c(a,b,"value")
    output$COMPACT <- renderPlot({
      ggplot(hm.data, aes_string(a, b)) +
        geom_tile(aes(fill = value)) +
        geom_text(aes(label = value)) +
        scale_fill_gradient(low = "white", high = "red")
    })
  })
}

shinyApp(ui = ui, server = server)