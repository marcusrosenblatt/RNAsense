#' @import ggplot2
#' @import parallel
#' @importFrom NBPSeq estimate.norm.factors prepare.nbp estimate.disp exact.nb.test
#' @import qvalue
#' @import SummarizedExperiment
#' @importFrom stats fisher.test pchisq time var
#' @importFrom utils write.table data
#' @importFrom methods is

#' @title Detect switching genes
#' @description For each gene, time-resolved RNA-seq measurements are analyzed for occurence of switches (up or down)
#'
#' @param dataset Object of class SummarizedExperiment, output of \link{SummarizedExperiment}, as assays use a numeric matrix with your RNAseq count data, rows correspond to different genes, columns correspond to different experiments, as rowData provide a \link{DataFrame} with columns name (geneID) and genename (the gene names), as colData provide a \link{DataFrame} with columns condition, time and replicate
#' @param experimentStepDetection Character, Name of condition for which switch detection is performed
#' @param pValueSwitch Numeric, pValue for switch detection
#' @param cores Numeric, Number of cores for parallelization, default 1 for no parallelization
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#'
#' @return Data.frame containing gene names and results of switch detection, information about switch time point and direction
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=20),]
#' resultSwitch <- getSwitch(dataset = mydata,
#' experimentStepDetection = "WT",
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' @export
getSwitch <- function(dataset = mydata, experimentStepDetection = "WT", pValueSwitch = 0.05, cores = 1, mytimes=times){
    stopifnot(is(dataset, "SummarizedExperiment"))
    data <- dataset[,colData(dataset)$condition==experimentStepDetection]
    do.call(rbind, mclapply(seq(1,ceiling(nrow(data)/500)), function(index){
        mydatasub <- data[seq((index-1)*500+1,min(index*500,nrow(data))),]
        do.call(rbind, lapply(seq(1,nrow(mydatasub)), function(i){
            temp <- data.frame(value = assays(mydatasub)[[1]][i,], time = as.numeric(colData(data)$time))
            out <- cbind(do.call(rbind, lapply(sort(unique(temp$time))[-length(sort(unique(temp$time)))], function(t){
                temp1 <- subset(temp, time <= t)$value
                temp2 <- subset(temp, time > t)$value
                data.frame(name=rowData(mydatasub)$name[i],
                           genename=rowData(mydatasub)$genename[i],
                           timepoint=t,
                           value=var(temp1)*(length(temp1)-1)/(length(temp1)) + var(temp2)*(length(temp2)-1)/length(temp2))
            })), var=var(temp$value)*(length(temp$value)-1)/7)
            out <- out[which(out$value==min(out$value))[1],]
            pValue <- NA
            if(out$value==0){out <- cbind(out, switch="none"); out$timepoint = NA} else {
                pValue <- pchisq(out$var/out$value,1)
                if ((1-pValue) < pValueSwitch) {
                    temp1 <- subset(temp, time <= out$timepoint)$value
                    temp2 <- subset(temp, time > out$timepoint)$value
                    if(mean(temp1) > mean(temp2)) out <- cbind(out, switch="down")
                        else out <- cbind(out, switch="up")
                } else {out <- cbind(out, switch="none"); out$timepoint = NA}
            }
            return(cbind(out[,c("name","genename", "timepoint", "switch")], pvalueSwitch = (1-pValue), experiment=experimentStepDetection))
        }))
    }, mc.cores = cores))
}

#' @title Detect switching genes
#' @description For each gene, time-resolved RNA-seq measurements are analyzed for occurence of switches (up or down)
#'
#' @param dataset Object of class SummarizedExperiment, output of \link{SummarizedExperiment}, as assays use a numeric matrix with your RNAseq count data, rows correspond to different genes, columns correspond to different experiments, as rowData provide a \link{DataFrame} with columns name (geneID) and genename (the gene names), as colData provide a \link{DataFrame} with columns condition, time and replicate
#' @param experimentStepDetection Character, Name of condition for which switch detection is performed
#' @param pValueSwitch Numeric, A threshold for counting cells as being invaded or not. When cells move towards negative z-direction, threshold should be negative.
#' @param cores Numeric, Number of cores for parallelization, default 1 for no parallelization
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#' @param chooseFirst boolean, if TRUE (default), the earliest time point is chosen for which a switch could be detected, if FALSE, the time point with the best likelihood for the one-step model is chosen
#'
#' @return Data.frame containing gene names and results of switch detection, information about switch time point and direction
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=10),]
#' resultSwitch <- getSwitch(dataset = mydata,
#' experimentStepDetection = "WT",
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' @export
getSwitchCorrect <- function(dataset = mydata, experimentStepDetection = "WT", pValueSwitch = 0.05, cores = 1, mytimes=times, chooseFirst=TRUE){
  stopifnot(is(dataset, "SummarizedExperiment"))
  data <- dataset[,colData(dataset)$condition==experimentStepDetection]
  do.call(rbind, mclapply(seq(1,ceiling(nrow(data)/500)), function(index){
    mydatasub <- data[seq((index-1)*500+1,min(index*500,nrow(data))),]
    do.call(rbind, lapply(seq(1,nrow(mydatasub)), function(i){
      temp <- data.frame(value = assays(mydatasub)[[1]][i,], time = as.numeric(colData(data)$time))
      out <- cbind(do.call(rbind, lapply(sort(unique(temp$time))[-length(sort(unique(temp$time)))], function(t){
        temp1 <- subset(temp, time <= t)$value
        temp2 <- subset(temp, time > t)$value
        data.frame(name=rowData(mydatasub)$name[i],
                   genename=rowData(mydatasub)$genename[i],
                   timepoint=t,
                   value=var(temp1)*(length(temp1)-1) + var(temp2)*(length(temp2)-1))
      })), var=var(temp$value)*(length(temp$value)-1))
      if(chooseFirst){
        out2 <- out[which(1-pchisq(out$var/out$value,1) < pValueSwitch),]
        pValue <- NA
        if(dim(out2)[1]==0){out <- cbind(out[1,], switch="none"); out$timepoint = NA} else {
          #print(out2)
          out <- out2[which(out$timepoint==min(out$timepoint))[1],]
          pValue <- pchisq(out$var/out$value,1)
          if ((1-pValue) < pValueSwitch) {
            temp1 <- subset(temp, time <= out$timepoint)$value
            temp2 <- subset(temp, time > out$timepoint)$value
            if(mean(temp1) > mean(temp2)) out <- cbind(out, switch="down")
            else out <- cbind(out, switch="up")
          } else {out <- cbind(out, switch="none"); out$timepoint = NA}
        }
        return(cbind(out[,c("name","genename", "timepoint", "switch")], pvalueSwitch = (1-pValue), experiment=experimentStepDetection))
      } else {
        out <- out[which(out$value==min(out$value))[1],]
        pValue <- NA
        if(out$value==0){out <- cbind(out[1,], switch="none"); out$timepoint = NA} else {
          pValue <- pchisq(out$var/out$value,1)
          if ((1-pValue) < pValueSwitch) {
            temp1 <- subset(temp, time <= out$timepoint)$value
            temp2 <- subset(temp, time > out$timepoint)$value
            if(mean(temp1) > mean(temp2)) out <- cbind(out, switch="down")
            else out <- cbind(out, switch="up")
          } else {out <- cbind(out, switch="none"); out$timepoint = NA}
        }
        return(cbind(out[,c("name","genename", "timepoint", "switch")], pvalueSwitch = (1-pValue), experiment=experimentStepDetection))
      }
    }))
  }, mc.cores = cores))
}



#' @title Detect switching genes (new version)
#' @description For each gene and for each time point, RNA-seq count data is analyzed for fold changes between two experimental conditions. This functions bases on functions from the R package NBPSeq package for fold change analysis
#'
#' @param dataset Object of class SummarizedExperiment, output of \link{SummarizedExperiment}, as assays use a numeric matrix with your RNAseq count data, rows correspond to different genes, columns correspond to different experiments, as rowData provide a \link{DataFrame} with columns name (geneID) and genename (the gene names), as colData provide a \link{DataFrame} with columns condition, time and replicate
#' @param experimentStepDetection Character, Name of condition for which switch detection is performed
#' @param cores Numeric, Number of cores for parallelization, default 1 for no parallelization
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#' @param pValueSwitch Numeric, pValue for switch detection
#'
#' @return Data.frame containing gene names, log fold change and p-values calculated from NBPSeq, each gene appears as often as available time points
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=20),]
#' resultFC <- getSwitchNew(dataset = mydata,
#' experimentStepDetection = "WT",
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6),
#' pValueSwitch=0.01)
#' @export
getSwitchNew <- function(dataset = mydata, experimentStepDetection = "WT", cores = 1, mytimes = times, pValueSwitch = 0.01){
  stopifnot(is(dataset, "SummarizedExperiment"))
  # auxiliary function getD
  getD <- function(value, FC, thFoldChange=NA, pValueSwitch=pValueSwitch){
    if(is.na(value) | is.na(FC)){return("none")} else {
      if(is.na(thFoldChange) & is.na(pValueSwitch)){
        return("none")
      } else if(!is.na(thFoldChange) & is.na(pValueSwitch)){
        if(FC < -log2(thFoldChange)) {return("down")}
        else if(FC > log2(thFoldChange)) {return("up")}
        else return("none")
      } else if(is.na(thFoldChange) & !is.na(pValueSwitch)){
        if(value < pValueSwitch){
          if(FC < 0) {return("down")}
          if(FC > 0) {return("up")}
        } else return("none")
      } else {
        if((value < pValueSwitch) & (FC < -log2(thFoldChange))) {return("down")}
        else if((value < pValueSwitch) & (FC > log2(thFoldChange))) {return("up")}
        else return("none")
      }
      
    }
  }
  out <- do.call(rbind, mclapply(mytimes[-length(mytimes)], function(t){
    data <- assays(dataset[,colData(dataset)$condition==experimentStepDetection])[[1]]
    
    ## Specify treatment groups
    grp.ids = as.character((colData(dataset)$time <= t))[colData(dataset)$condition==experimentStepDetection]  # Numbers or strings are both OK
    
    ## Estimate normalization factors
    norm.factors = estimate.norm.factors(data);
    
    ## Prepare an NBP object, adjust the library sizes by thinning the counts.
    set.seed(999)
    obj = prepare.nbp(data, grp.ids, lib.sizes=colSums(data), norm.factors=norm.factors, print.level = 0);
    
    ## Fit a dispersion model (NBQ by default)
    obj = estimate.disp(obj, print.level = 0);
    
    ## Perform exact NB test
    grp1 = "TRUE";
    grp2 = "FALSE";
    
    obj = exact.nb.test(obj, grp1, grp2, print.level = 0);
    
    # ## Output results
    out <- data.frame(name=rowData(dataset)$name,
                      genename=rowData(dataset)$genename,
                      logFoldChangeSwitch = obj$log.fc,
                      pvalueSwitch = obj$p.values,
                      timepoint=t,
                      experiment=experimentStepDetection)
    cbind(out, switch=vapply(seq(1,dim(out)[1]), function(i){
      getD(out$pvalueSwitch[i],
           out$logFoldChangeSwitch[i],
           thFoldChange = 2,  ## ignored, if NA
           pValueSwitch = pValueSwitch  ## ignored, if NA
      )
    }, c("up")))
  }, mc.cores = cores))
  out <- do.call(rbind, mclapply(unique(out$name), function(myname){
    sub <- subset(out, name==myname)
    if(Reduce("&",sub$switch=="none")){
      sub[1,]
    } else {
      temp <- data.frame(value = assays(dataset[rowData(dataset)$name==myname,
                                                colData(dataset)$condition==experimentStepDetection])[[1]][1,],
                         time = as.numeric(colData(dataset)$time[which(colData(dataset)$condition==experimentStepDetection)]))
      modelvalue <- NA
      tout <- sub$timepoint[1]
      for(myt in sub$timepoint){
        temp1 <- subset(temp, time <= myt)$value
        temp2 <- subset(temp, time > myt)$value
        if(!is.na(modelvalue)){
          testvalue <- var(temp1)*(length(temp1)-1) + var(temp2)*(length(temp2)-1)
          if(modelvalue > testvalue){
            tout <- myt
            modelvalue <- testvalue
          }
        } else modelvalue <- var(temp1)*(length(temp1)-1) + var(temp2)*(length(temp2)-1)
      }
      sub[which(sub$timepoint==tout),]
    }
  }, mc.cores = cores))
  return(out)
}

#' @title Detect fold changes
#' @description For each gene and for each time point, RNA-seq count data is analyzed for fold changes between two experimental conditions. This functions bases on functions from the R package NBPSeq package for fold change analysis
#'
#' @param dataset Object of class SummarizedExperiment, output of \link{SummarizedExperiment}, as assays use a numeric matrix with your RNAseq count data, rows correspond to different genes, columns correspond to different experiments, as rowData provide a \link{DataFrame} with columns name (geneID) and genename (the gene names), as colData provide a \link{DataFrame} with columns condition, time and replicate
#' @param myanalyzeConditions Character vector, Name of experimental conditions
#' @param cores Numeric, Number of cores for parallelization, default 1 for no parallelization
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#' @param pValueFC Numeric, p-value for fold change detection
#'
#' @return Data.frame containing gene names, log fold change and p-values calculated from NBPSeq, each gene appears as often as available time points
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=20),]
#' resultFC <- getFC(dataset = mydata,
#' myanalyzeConditions = c("WT", "MZsox"),
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6),
#' pValueFC=0.01)
#' @export
getFC <- function(dataset = mydata, myanalyzeConditions = analyzeConditions, cores = 1, mytimes = times, pValueFC = 0.01){
    stopifnot(is(dataset, "SummarizedExperiment"))
    # auxiliary function getD
    getD <- function(value, FC, thFoldChange=NA, pValueFC=0.01){
        if(is.na(value) | is.na(FC)){return("none")} else {
            if(is.na(thFoldChange) & is.na(pValueFC)){
                return("none")
        } else if(!is.na(thFoldChange) & is.na(pValueFC)){
            if(FC < -log2(thFoldChange)) {return(paste0(myanalyzeConditions[1], "<", myanalyzeConditions[2]))}
                else if(FC > log2(thFoldChange)) {return(paste0(myanalyzeConditions[1], ">", myanalyzeConditions[2]))}
                    else return(paste0(myanalyzeConditions[1], "=", myanalyzeConditions[2]))
        } else if(is.na(thFoldChange) & !is.na(pValueFC)){
            if(value < pValueFC){
                if(FC < 0) {return(paste0(myanalyzeConditions[1], "<", myanalyzeConditions[2]))}
                if(FC > 0) {return(paste0(myanalyzeConditions[1], ">", myanalyzeConditions[2]))}
            } else return(paste0(myanalyzeConditions[1], "=", myanalyzeConditions[2]))
        } else {
            if((value < pValueFC) & (FC < -log2(thFoldChange))) {return(paste0(myanalyzeConditions[1], "<", myanalyzeConditions[2]))}
                else if((value < pValueFC) & (FC > log2(thFoldChange))) {return(paste0(myanalyzeConditions[1], ">", myanalyzeConditions[2]))}
                    else return(paste0(myanalyzeConditions[1], "=", myanalyzeConditions[2]))
        }

        }
    }
    out <- do.call(rbind, mclapply(mytimes, function(t){
        data <- assays(dataset[,colData(dataset)$time==t])[[1]]

        ## Specify treatment groups
        grp.ids = colData(dataset)$condition[colData(dataset)$time==t]  # Numbers or strings are both OK

        ## Estimate normalization factors
        norm.factors = estimate.norm.factors(data);

        ## Prepare an NBP object, adjust the library sizes by thinning the counts.
        set.seed(999)
        obj = prepare.nbp(data, grp.ids, lib.sizes=colSums(data), norm.factors=norm.factors, print.level = 0);

        ## Fit a dispersion model (NBQ by default)
        obj = estimate.disp(obj, print.level = 0);

        ## Perform exact NB test
        grp1 = myanalyzeConditions[2];
        grp2 = myanalyzeConditions[1];

        obj = exact.nb.test(obj, grp1, grp2, print.level = 0);

        # ## Output results
        out <- data.frame(name=rowData(dataset)$name,
                          genename=rowData(dataset)$genename,
                          logFoldChange = obj$log.fc,
                          pValue = obj$p.values,
                          time=t)
        cbind(out, FCdetect=vapply(seq(1,dim(out)[1]), function(i){
            getD(out$pValue[i],
                out$logFoldChange[i],
                thFoldChange = 2,  ## ignored, if NA
                pValueFC = pValueFC  ## ignored, if NA
            )
        }, c("WT > condition")))
    }, mc.cores = cores))
    out$name <- rep(rowData(dataset)$name, length(mytimes))
    return(out)
}

#' @title Combine results
#' @description Results of switch and fold change analysis are collected in one data.frame
#'
#' @param myresultSwitch data.frame, output of \link{getSwitch}
#' @param myresultFC data.frame, output of \link{getFC}
#' @param nrcores Numeric, Number of cores for parallelization, default 1 for no parallelization
#'
#' @return Data.frame containing information on switch and fold change detection for each gene
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=20),]
#' resultFC <- getFC(dataset = mydata,
#' myanalyzeConditions = c("WT", "MZsox"),
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' resultSwitch <- getSwitch(dataset = mydata,
#' experimentStepDetection = "WT",
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' combineResults(resultSwitch, resultFC)
#' @export
combineResults <- function(myresultSwitch = resultSwitch, myresultFC = resultFC, nrcores = 1){
    ## auxiliary function getFCupdown
    getFCupdown <- function(gene, myresultFC = resultFC){
        gene <- factor(gene, levels=levels(myresultFC$name))
        sub <- subset(myresultFC, name==gene)
        f <- subset(sub, grepl(">", FCdetect))$time
        g <- subset(sub, grepl("<", FCdetect))$time
        f <- as.character(format(f, nsmall = 1))
        g <- as.character(format(g, nsmall = 1))
        if(length(f)>0){f <- do.call(paste,as.list(paste0(f, "hpf")))} else{f <- ""}
        if(length(g)>0){g <- do.call(paste,as.list(paste0(g, "hpf")))} else{g <- ""}
        data.frame(FCdown=f,FCup=g)
    }
    temp <- do.call(rbind, mclapply(myresultSwitch$genename, function(gene){
        getFCupdown(gene, myresultFC)
    }, mc.cores = nrcores))
    return(cbind(resultSwitch, temp))
}

#' @title plot SSGS gene classes
#' @description Genes are sorted into groups with respect to switch time and time point of fold change detection. For each group, results of wild type and knockdown-condition are compared by means of fisher's exact test to show whether the knocked down gene enhances or suppresses the respective gene group.
#'
#' @param myresultCombined data.frame, output of \link{combineResults}
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#' @param myanalyzeConditions character vector, the conditions that were analyzed
#'
#' @return SSGS color plot in ggplot format
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' library(ggplot2)
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=20),]
#' resultFC <- getFC(dataset = mydata,
#' myanalyzeConditions = c("WT", "MZsox"),
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' resultSwitch <- getSwitch(dataset = mydata,
#' experimentStepDetection = "WT",
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' resultCombined <- combineResults(resultSwitch, resultFC)
#' plotSSGS(myresultCombined = resultCombined,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6),
#' myanalyzeConditions = c("WT", "MZsox"))
#' @export
plotSSGS <- function(myresultCombined = resultCombined, mytimes = times, myanalyzeConditions = analyzeConditions){
    # auxiliary function for application of fisher test
    getFT <- function(myresult=result, myswitch="up",switchtime=3, xaxis="2.5hpf", identifier="FCdown"){
        if(identifier == "FCdown"){
            a=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                             grepl(xaxis, FCdown)))[1]
            b=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                             grepl(xaxis, FCdown)))[1]
            c=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                             !grepl(xaxis, FCdown)))[1]
            d=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                             !grepl(xaxis, FCdown)))[1]
        } else {
            a=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                             grepl(xaxis, FCup)))[1]
            b=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                             grepl(xaxis, FCup)))[1]
            c=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                             !grepl(xaxis, FCup)))[1]
            d=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                             !grepl(xaxis, FCup)))[1]
        }

        pvalue_suppress <- fisher.test(rbind(c(a,b), c(c,d)), alternative="less")$p.value
        pvalue_enhance <- fisher.test(rbind(c(a,b), c(c,d)), alternative="greater")$p.value
        if(pvalue_suppress < 1e-20) cluster <- "1E-20 Suppression" else
            if(pvalue_suppress < 1e-10) cluster <- "1E-10 Suppression" else
                if(pvalue_suppress < 1e-5) cluster <- "1E-05 Suppression" else
                    if(pvalue_suppress < 1e-2) cluster <- "1E-02 Suppression" else
                        if(pvalue_enhance < 1e-20) cluster <- "1E-20 Enhancement" else
                            if(pvalue_enhance < 1e-10) cluster <- "1E-10 Enhancement" else
                                if(pvalue_enhance < 1e-5) cluster <- "1E-05 Enhancement" else
                                    if(pvalue_enhance < 1e-2) cluster <- "1E-02 Enhancement" else cluster <- "none"
        data.frame(pvalue_enhance=pvalue_enhance, pvalue_suppress=pvalue_suppress, cluster=cluster, myswitch=myswitch)
    }

    out <- do.call(rbind, lapply(mytimes, function(t){
        do.call(rbind, lapply(paste0(format(mytimes, nsmall = 1),"hpf"), function(x){
            do.call(rbind, lapply(c("FCdown", "FCup"), function(ident){
                do.call(rbind, lapply(c("up", "down"), function(myswitch){
                    cbind(getFT(myresult=subset(myresultCombined, experiment=="WT"), myswitch=myswitch, switchtime = t, xaxis = x, identifier = ident),
                          time=t, xaxis=x, identifier=ident, experiment="WT")
                }))
            }))
        }))
    }))
    out$cluster <- factor(out$cluster, levels = c("none", "1E-20 Enhancement","1E-10 Enhancement", "1E-05 Enhancement","1E-02 Enhancement",
                                                  "1E-20 Suppression","1E-10 Suppression", "1E-05 Suppression", "1E-02 Suppression"))

    mylabeller <- c("FCdown" = paste0(myanalyzeConditions[1]," > ",myanalyzeConditions[2]),
                    "FCup" = paste0(myanalyzeConditions[2]," > ",myanalyzeConditions[1]),
                    "up" = "up",
                    "down" = "down")

    P <- ggplot(out, aes(x=xaxis, y=time, fill=cluster)) + geom_tile(color="black") +
        xlab("Stage-specific gene sets WT") + ylab("Switch Time") +theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
        scale_x_discrete(breaks=paste0(format(seq(2.5,6, by=0.5),nsmall=1), "hpf"),labels=format(seq(2.5,6, by=0.5),nsmall=1)) +
        scale_y_reverse(breaks=seq(2.5,6, by=0.5)) +
        facet_grid(myswitch~identifier, labeller = as_labeller(mylabeller)) +
        scale_fill_manual(name="Fisher Test Result",
                            values=c("none" = "grey", "1E-20 Enhancement" = "darkred", "1E-10 Enhancement" = "red", "1E-05 Enhancement"= "orange",
                                                    "1E-02 Enhancement" = "yellow", "1E-20 Suppression"= "violet", "1E-10 Suppression" = "darkblue",
                                                    "1E-05 Suppression" = "blue", "1E-02 Suppression" = "lightblue" ))
    return(P)
}

#' @title Output gene tables
#' @description Output information on switching genes (up/down) in tabular format (gene identifier/gene name) are created as .txt file and written to the specified working directory.
#' Two of the generated files (geneNamelist) contain gene lists with gene name for genes that switch up and down respectively. The other two (genelist) contain exactly the same output but with gene identifiers instead of gene names depending on what you prefer for further analysis. Each column corresponds to a combination of switch time point, fold change direction and time point of fold change. All genes for which fold change was detected at the indicated time point and switch was detected at the indicated time point are listed in the corresponding column. Note that a single gene may appear multiple times. The fifth .txt file (switchList) contains information on detected switches in a different format. The output consists of table with six columns with each row corresponding to one gene. Detected switches are indicated by 1, -1 and 0 for switch up, switch down and no switch, respectively. If a switch was detected, the column timepoint indicated the corresponding time point of switch detection.
#'
#' @param myresultCombined data.frame, output of \link{combineResults}
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#' @param myanalyzeConditions character vector, the conditions that were analyzed
#' @param mywd character, working directory to which results will be written, if NULL the current working directory is used
#'
#' @return Working directory where results have been written to
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' @examples
#' library(ggplot2)
#' data(MZsox)
#' mydata <- MZsox[seq(1,nrow(MZsox), by=20),]
#' resultFC <- getFC(dataset = mydata,
#' myanalyzeConditions = c("WT", "MZsox"),
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' resultSwitch <- getSwitch(dataset = mydata,
#' experimentStepDetection = "WT",
#' cores = 1,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6))
#' resultCombined <- combineResults(resultSwitch, resultFC)
#' outputGeneTables(resultCombined,
#' mytimes = c(2.5,3,3.5,4,4.5,5,5.5,6),
#' myanalyzeConditions = c("WT", "MZsox"))
#' @export
outputGeneTables <- function(myresultCombined = resultCombined, mytimes = times, myanalyzeConditions = analyzeConditions, mywd = NULL){
    if(is.null(mywd)) mywd <- getwd()
    genetableUp <- c()
    genetableDown <- c()
    geneNametableUp <- c()
    geneNametableDown <- c()
    for (t in mytimes[-length(mytimes)]){
        for(identifier in c("FCdown", "FCup")){
            for(x in paste0(format(mytimes, nsmall = 1),"hpf")){
                    if(identifier=="FCdown"){
                        genesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCdown))$name))
                        genesUp <- c(genesUp, rep("",3000-length(genesUp)))
                        genetableUp <- cbind(genetableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(myanalyzeConditions[1]," > ",myanalyzeConditions[2]),x,genesUp))
                        geneNamesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCdown))$genename))
                        geneNamesUp <- c(geneNamesUp, rep("",3000-length(geneNamesUp)))
                        geneNametableUp <- cbind(geneNametableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(myanalyzeConditions[1]," > ",myanalyzeConditions[2]),x,geneNamesUp))

                        genesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCdown))$name))
                        genesDown <- c(genesDown, rep("",3000-length(genesDown)))
                        genetableDown <- cbind(genetableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(myanalyzeConditions[1]," > ",myanalyzeConditions[2]),x,genesDown))
                        geneNamesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCdown))$genename))
                        geneNamesDown <- c(geneNamesDown, rep("",3000-length(geneNamesDown)))
                        geneNametableDown <- cbind(geneNametableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(myanalyzeConditions[1]," > ",myanalyzeConditions[2]),x,geneNamesDown))
                    }
                    if(identifier=="FCup"){
                        genesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCup))$name))
                        genesUp <- c(genesUp, rep("",3000-length(genesUp)))
                        genetableUp <- cbind(genetableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(myanalyzeConditions[2]," > ",myanalyzeConditions[1]),x,genesUp))
                        geneNamesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCup))$genename))
                        geneNamesUp <- c(geneNamesUp, rep("",3000-length(geneNamesUp)))
                        geneNametableUp <- cbind(geneNametableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(myanalyzeConditions[2]," > ",myanalyzeConditions[1]),x,geneNamesUp))

                        genesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCup))$name))
                        genesDown <- c(genesDown, rep("",3000-length(genesDown)))
                        genetableDown <- cbind(genetableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(myanalyzeConditions[2]," > ",myanalyzeConditions[1]),x,genesDown))
                        geneNamesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCup))$genename))
                        geneNamesDown <- c(geneNamesDown, rep("",3000-length(geneNamesDown)))
                        geneNametableDown <- cbind(geneNametableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(myanalyzeConditions[2]," > ",myanalyzeConditions[1]),x,geneNamesDown))
                    }
            }
        }
    }

    nametag <- myanalyzeConditions[2]

    write.table(genetableUp, file=paste0(mywd,"/genelist_switchUp_",nametag,".txt"), sep = "\t", row.names=FALSE, col.names = FALSE)
    write.table(genetableDown, file=paste0(mywd,"/genelist_switchDown_",nametag,".txt"), sep = "\t", row.names=FALSE, col.names = FALSE)
    write.table(geneNametableUp, file=paste0(mywd,"/geneNamelist_switchUp_",nametag,".txt"), sep = "\t", row.names=FALSE, col.names = FALSE)
    write.table(geneNametableDown, file=paste0(mywd,"/geneNamelist_switchDown_",nametag,".txt"), sep = "\t", row.names=FALSE, col.names = FALSE)

    resultOut1 <- subset(myresultCombined[,c("name","genename", "pvalueSwitch", "switch", "timepoint", "experiment")], experiment==myanalyzeConditions[1])
    resultOut2 <- subset(myresultCombined[,c("name","genename", "pvalueSwitch", "switch", "timepoint", "experiment")], experiment==myanalyzeConditions[2])
    resultOut1$timepoint <- format(resultOut1$timepoint, nsmall = 1)
    resultOut2$timepoint <- format(resultOut2$timepoint, nsmall = 1)
    if(dim(resultOut2)[1] > 0) resultOut <- cbind(resultOut1, resultOut2[,-1]) else resultOut <- resultOut1
    resultOut$switch <- sub("none", "0", resultOut$switch)
    resultOut$switch <- sub("up", "1", resultOut$switch)
    resultOut$switch <- sub("down", "-1", resultOut$switch)
    write.table(resultOut, file=paste0(mywd,"/switchList_",nametag,".txt"), sep = "\t", row.names=FALSE)
    return(paste("Results written to", mywd))
}
