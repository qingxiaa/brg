#' @title Shiny pwrEWAS
#'
#' @description pwrBRIDGE pwrBRIDGE performs power analysis and bridging sample size calculation via empirical simulation for batch-confounded longitudinal study. 
#' pwrBRIDGE generates data from multivariate normal distribution considering the subject and sample dependency with batch effects incorporating based 
#' on location and scale (L/S) model. Then the batch effects are corrected using any of BRIDGE, ComBat and longitudinal ComBat. The corrected data are 
#' used for downstream hypothesis testing and power/bridging sample size assessment.
#' 
#' @return pwrBRIDGE_shiny user-interface
#' 
#' @export

#'@import dplyr shiny ggplot2 gridExtra BiocManager shinyBS sva longCombat




pwrBRIDGE_shiny <- function(){

methodslist <- c("BRIDGE-1","BRIDGE-2","ComBat","LongComBat")

###############################################################Shiny Start Here #####################################################################################
ui <- shinyUI(
  fluidPage(
    tags$head(tags$style(HTML(".shiny-notification {
            height: 150px;
            width: 400px;
            position:fixed;
            font-size: 150%;
            top: calc(50% - 35px);;
            left: calc(50% - 100px);;}"))),
    tags$style(type='text/css', '#log {text-align: left;}'),
    
    titlePanel("pwrBRIDGE: Power Analysis and Bridging Sample Size Calculations"),
    HTML("pwrBRIDGE pwrBRIDGE performs power analysis and bridging sample size calculation via empirical simulation for batch-confounded longitudinal study. 
              pwrBRIDGE generates data from multivariate normal distribution considering the subject and sample dependency with batch effects incorporating based 
              on location and scale (L/S) model. Then the batch effects are corrected using any of BRIDGE, ComBat and longitudinal ComBat. The corrected data are 
              used for downstream hypothesis testing and power/bridging sample size assessment. Detailed description of in-/outputs, instructions and an example, 
              as well as interpretations of the example results are provided in the following paper: "),
    tags$a(href="https://www.dropbox.com/s/ru9z82z5c7pqc39/Paper%202_1005.docx?dl=0", "pwrBRIDGE"),
    
    
    HTML("</br></br>Authors: Qing Xia, Jeffrey A. Thompson, Devin C. Koestler </br>"),
    HTML("Department of Biostatistics and Data Science, University of Kansas School of Medicine"),
    
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("Methods","Choose methods:", choices = methodslist,selected = methodslist[3]),
        numericInput("mbe", "Multiplicative Batch Effects:", value = 1.5, step = NA), 
        numericInput("abe", "Additive Batch Effects:", value = 3, step = NA),
        numericInput("es", "Effect Size:", value = 0.5,step = NA), 
        numericInput("sd1", "Standard deviation for timepoint 1", value = 0.5,step = NA),
        numericInput("sd2", "Standard deviation for timepoint 2", value = 0.5,step = NA),
        numericInput("rho1", "Within-sample Correlation between batches", value = 0.9,min= 0, max=1,step = NA),
        numericInput("rho2", "Within-subject Correlation between timepoints", value = 0.5,min= 0, max=1,step = NA),
        numericInput("J", "Total Genes", value = 1000, min= 200, max=5000, step = NA),
        numericInput("G", "Total Differential Genes", value = 200, min= 100,max=5000,step = NA),
        numericInput("N", "Total Participants", value = 50, min=20, step = NA),
        numericInput("t1", "FDR", value = 0.05, max=1, step = NA),
        numericInput("pwr", "Desired Power", value = 0.8, min=1, step = NA),
        fluidRow(HTML('<h5><b>Bridging Sample Size Range</b></h5>'),
                 
                 column(4, numericInput("m1","From", value = 5)),
                 
                 column(4, numericInput("m2","To", value = 15)),
                 
                 column(4, numericInput("by","by", value = 5))
        ),
        numericInput("iter", "Iterations", value = 5,max=100,step = NA),
        numericInput("seed", "set seed:", value = 2021,step = NA),
        
        submitButton("Start!"),
        bsTooltip("mbe", "The multiplicative batch effects here is a ratio parameter describing the ratio of multiplicative batch effect (scaling the variance of gene expressioin) between batch 1 and batch 2, e.g. 2 means multiplicative batch effects of batch 2 is twice as that of batch 1",
                  "right", options = list(container = "body")),
        bsTooltip("abe", "The additive batch effects here is a parameter describing the difference of additive batch effect (shift the variance of gene expressioin) between batch 1 and batch 2, e.g. 2 means additive batch effects of batch 2 is 2 more than that of batch 1 in log 2 fold change",
                  "right", options = list(container = "body")),
        bsTooltip("es", "The effect size here is a difference in log 2 fold change between two time points, e.g. 2 means gene expression level for a single gene of time point 2 is 2 more than that of time point 1 in log 2 fold change",
                  "right", options = list(container = "body")),
        bsTooltip("sd1", "The standard deviation of gene expression level for timepoint 1",
                  "right", options = list(container = "body")),
        bsTooltip("sd2", "The standard deviation of gene expression level for timepoint 2",
                  "right", options = list(container = "body")),
        bsTooltip("rho1", "The correlation coefficient between bridging data (repeated measurements of same sample), that is within-sample correlation between batches",
                  "right", options = list(container = "body")),
        bsTooltip("rho2", "The correlation coefficient between two timepoints (longitudinal data of same subject), that is within-subject correlation between timepoints",
                  "right", options = list(container = "body")),
        bsTooltip("J", "Total number of genes. We recommend a value between 500~1000. A large number results in longer running time.",
                  "right", options = list(container = "body")),
        bsTooltip("G", "The number of differential genes. We recommend a value > 100, e.g. 200 out of 1000 genes are differentially expressed.",
                  "right", options = list(container = "body")),
        bsTooltip("N", "The number of subjects, e.g. 100 subjects are followed at two timepoints",
                  "right", options = list(container = "body")),
        bsTooltip("t1", "The cutoff type I error rate used for calculating type I error, e.g. FDR = 0.05, we calculate the proortion of null genes being detected as differentially expressed, as type I error.",
                  "right", options = list(container = "body")),
        bsTooltip("pwr", "The desired power as a horizontal line in power plot",
                  "right", options = list(container = "body")),
      ),
      
      mainPanel(
        
        tableOutput("table"),
        plotOutput("plot"),
        textOutput("time")
        
      )
    )
  )
)


# Define server
server <- function(session,input, output) {
  Inference_comparison <- function(para, setSeeds, method = 1, iter){
    set.seed(setSeeds)
    nprogress <- input$iter*(input$m2-input$m1)/input$by
    Stats_inference <- function(df, imputation, geneind, alpha){
      
      Obs_Y1t1 <- with(df, which(batch == "Batch1"& time == "Time1"))
      Obs_Y2t1 <- with(df, which(batch %in% c("Batch2","Imputed") & time == "Time1"))
      Obs_Y2t2 <- with(df, which(batch == "Batch2"& time == "Time2"))
      
      Y1t1 <- df[Obs_Y1t1, -c(1:3)]
      Y2t1 <- df[Obs_Y2t1, -c(1:3)]
      Y2t2 <- df[Obs_Y2t2, -c(1:3)]
      
      G <- length(geneind)
      N <- nrow(Y1t1)
      J <- ncol(Y1t1)
      pvalue <-  NULL      
      if (imputation == FALSE){ # method 1 or ComBat y1t1 vs y2t2
        pvalue <- sapply(1:J,function(x) t.test(Y1t1[,x],Y2t2[,x],paired=TRUE,alternative = "two.sided")$p.value)      
      }
      else if(imputation == TRUE){ # method 2 y2t1 vs y2t2
        pvalue <- sapply(1:J,function(x) t.test(Y2t1[,x],Y2t2[,x],paired=TRUE,alternative = "two.sided")$p.value)
      }
      # rejections based on the adjusted p values
      rej <- p.adjust(pvalue, method = "hochberg", n = J) < alpha
      
      err1 <- sum(rej[-geneind])/(J-G)
      pw <- sum(rej[geneind])/G
      return(list("TypeI" = err1,
                  "Power" = pw))
    }
    err1.df <- pw.df <- NULL
    for (i in 1:iter){
      set.seed(setSeeds+i)
      # method == 1--BRIDGE-1; method ==  2--BRIDGE-2; method == 3--ComBat; method == 4--longComBat
      # Simulate data
      simulation.out <- do.call(SimulateData, para)
      # Data N*J
      data <- t(simulation.out[[1]][,-c(1:3)])
      
      # bridging samples index and truly differential genes index
      bridging.sample.index <- simulation.out[[2]]
      gene.index <- simulation.out[[3]]
      
      N = para$N
      M = para$M
      id<-c(c(1:N),bridging.sample.index,c(1:N))
      batch<-c(rep("Batch1", N), rep("Batch2",N+M))
      time<-c(rep("Time1",N+M), rep("Time2",N))
      
      factor.df<-data.frame(id,batch,time)
      un.df <- data.frame(factor.df,t(data))
      
      if (method == 1){
        ## Adjuste Data
        # BRIDGE-1 adjusted data with batch, time, id information
        method1.out <- BRIDGE_pwr(infile=data, id, batch, time, data=factor.df,imputation = FALSE, optimization = FALSE)
        method1.df <- data.frame(factor.df, t(method1.out[[2]]))
        err1.df[i] <- Stats_inference(method1.df, imputation = FALSE, gene.index, input$t1)[[1]]
        pw.df[i] <- Stats_inference(method1.df, imputation = FALSE, gene.index, input$t1)[[2]]
      }
      # BRIDGE-2 adjusted data with batch, time, id information
      if (method == 2){
        method2.out <- BRIDGE_pwr(infile=data, id, batch, time,data=factor.df, imputation = TRUE, optimization = FALSE)
        factor.df2 <- method2.out[[1]]
        method2.df <- data.frame(factor.df2, t(method2.out[[2]]))
        err1.df[i] <- Stats_inference(method2.df, imputation = TRUE, gene.index, input$t1)[[1]]
        pw.df[i] <- Stats_inference(method2.df, imputation = TRUE, gene.index, input$t1)[[2]]
      }
      # Combat adjusted data with batch, time, id information
      if (method == 3){
        mdl <- model.matrix(~time, data = factor.df)
        combat.adj <- ComBat(dat = data,
                             batch=factor.df$batch,
                             mod=mdl,
                             par.prior=TRUE,
                             prior.plot=FALSE)
        combat.df <- data.frame(factor.df, t(combat.adj))
        err1.df[i] <- Stats_inference(combat.df, imputation = FALSE, gene.index, input$t1)[[1]]
        pw.df[i] <- Stats_inference(combat.df, imputation = FALSE, gene.index, input$t1)[[2]]
        
      }
      if (method == 4){
        # longitudinal Combat
        
        lcombat.out <- longCombat(idvar='id',
                                  timevar='time',
                                  batchvar='batch',
                                  features=colnames(un.df)[-c(1:3)],
                                  formula="time",
                                  ranef='(1|id)',
                                  data=un.df)
        
        lcombat.df <- lcombat.out$data_combat
        err1.df[i] <- Stats_inference(lcombat.df, imputation = FALSE, gene.index, input$t1)[[1]]
        pw.df[i] <- Stats_inference(lcombat.df, imputation = FALSE, gene.index, input$t1)[[2]]
        
      }
      incProgress(1/nprogress, detail = paste("Doing M = ",para$M," iteration", i,". The run time is estimated at ", es_runtime()/60, "mins"))
      
    }
    
    return(list("TypeI"= mean(err1.df), "Power" = mean(pw.df)) )
  }
  # estimate run time
  es_runtime <- reactive({
    t(as.integer(methodslist %in% input$Methods)*input$iter*input$J/1000*((input$m2-input$m1)/input$by+1))%*%c(10.045,13.019,0.487,21.533)
  })
  
  pwrBRIDGE <- function(params, setSeeds, method = 1,iter=20){
    withProgress(message = paste(methodslist[method]), value = 0, {
      #nprogress = nrow(params)
      HPC_out <- NULL
      for(n in 1:nrow(params)){
        para <- as.list(params[n,])
        HPC_out <- rbind(data.frame(HPC_out),data.frame(Inference_comparison(para, setSeeds, method=method,iter = iter)))
      } 
      return(HPC_out)
    })
  }
  
  create.table <- reactive({
    startTime <- Sys.time()
    params <- expand.grid(J = input$J,
                          N = input$N,
                          G = input$G,
                          M = seq(input$m1,input$m2,by=input$by),
                          TimeEff0 = input$es,
                          mu01 = 2,
                          sigma01 = input$sd1,
                          sigma02 = input$sd2,
                          delta2 = input$mbe,
                          rho1 = input$rho1,
                          rho2 = input$rho2,
                          phi2 = input$abe)
    
    
    
    b1 <- b2 <- cb <- lcb  <-  NULL
    
    
    if ("BRIDGE-1" %in% input$Methods){
      method = 1
      b1 <- pwrBRIDGE(params, input$seed, method = 1,iter = input$iter)
      b1$M <- params$M
      b1$methods <- "BRIDGE-1"
    }
    if ("BRIDGE-2" %in% input$Methods){
      method = 2
      b2 <- pwrBRIDGE(params, input$seed, method = 2,iter = input$iter)
      b2$M <- params$M
      b2$methods <- "BRIDGE-2"
    }
    if ("ComBat" %in% input$Methods){
      method = 3
      cb <- pwrBRIDGE(params, input$seed, method = 3,iter = input$iter)
      cb$M <- params$M
      cb$methods <- "ComBat"
    }
    if ("LongComBat" %in% input$Methods){
      method = 4
      lcb <- pwrBRIDGE(params, input$seed, method = 4,iter = input$iter)
      lcb$M <- params$M
      lcb$methods <- "LongComBat"
    }
    
    out <- rbind(b1,b2,cb,lcb)
    endTime <- difftime(Sys.time(), startTime, units = "mins") 
    
    
    out <- list("out"=out, "runtime"=endTime)
  })
  create.plot <- reactive({
    
    p1<-ggplot(create.table()$out, aes(x = M, y=TypeI, group = methods, colour = methods)) +
      geom_line(size=1) +
      geom_point()+
      # geom_errorbar(aes(ymin=TypeI-sd, ymax=TypeI+sd), width=.5,position=position_dodge(0.05))+
      theme_classic() +
      geom_text(data = create.table()$out, aes(label=round(TypeI,3), y= TypeI+0.03),position = position_dodge(3),size = 4)+
      labs(title="Type I error rate as bridging sample size increases")
    p2<-ggplot(create.table()$out, aes(x = M, y=Power, group = methods, colour = methods)) +
      geom_line(size=1) +  geom_hline(yintercept=input$pwr,linetype="dashed", color = "black",size=1) +
      
      geom_point()+
      # geom_errorbar(aes(ymin=TypeI-sd, ymax=TypeI+sd), width=.5,position=position_dodge(0.05))+
      theme_classic() +
      geom_text(data = create.table()$out, aes(label=round(Power,3), y= Power+0.01),position = position_dodge(3),size = 4)+
      labs(title="Power as bridging sample size increases")
    
    p=grid.arrange(p1,p2,ncol=1)
    return(p)
  })
  
  
  output$table = renderTable({
    create.table()$out 
  })
  
  
  
  output$plot = renderPlot({
    create.plot()  
  })
  output$time <- renderText(print (paste("Runing time is ",round(create.table()$runtime,3),' mins')))
  session$onSessionEnded(function() {
    stopApp()
  })
}

shinyApp(ui = ui, server = server)
}