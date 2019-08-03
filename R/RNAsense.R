#' Detect switching genes
#'
#' @description For each gene, time-resolved RNA-seq measurements are analyzed for occurence of switches (up or down)
#'
#' @param dataset data.frame, rows correspond to different genes, first column contains gene identifiers, second column contains the gene name, columns 3 to n contain the RNAseq count data, column names should start with the condition identifier (e.g. "WT") followed by the time separated by "_"
#' @param experimentStepDetection Character, Name of condition for which switch detection is performed
#' @param pValueSwitch Numeric, A threshold for counting cells as being invaded or not. When cells move towards negative z-direction, threshold should be negative.
#' @param cores Numeric, Number of cores for parallelization, default 1 for no parallelization
#'
#' @return Data.frame containing gene names and results of switch detection, information about switch time point and direction
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}

getSwitch <- function(dataset = mydata, experimentStepDetection = "WT", pValueSwitch = 0.05, cores = 1){
  data <- mydata[,c(1,2,which(grepl(experimentStepDetection, colnames(dataset))))]
  mytimes <- as.numeric(do.call(rbind,strsplit(names(data)[-c(1,2)], "_"))[,2])
  do.call(rbind, mclapply(1:ceiling(nrow(data)/500), function(index){
    mydatasub <- subset(data, name%in%unique(data$name)[seq((index-1)*500+1,min(index*500,nrow(data)))])
    do.call(rbind, lapply(1:nrow(mydatasub), function(i){
      temp <- data.frame(value = as.numeric(mydatasub[i,-c(1,2)]), time = mytimes)
      out <- cbind(do.call(rbind, lapply(sort(unique(temp$time))[-length(sort(unique(temp$time)))], function(t){
        temp1 <- subset(temp, time <= t)$value
        temp2 <- subset(temp, time > t)$value
        data.frame(name=mydatasub[i,1], genename=mydatasub[i,2], timepoint=t, value=var(temp1)*(length(temp1)-1) + var(temp2)*(length(temp2)-1))
      })), var=var(temp$value)*(length(temp$value)-1))
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

#' Detect fold changes
#'
#' @description For each gene and for each time point, RNA-seq count data is analyzed for fold changes between two experimental conditions. This functions bases on functions from the R package NBPSeq package for fold change analysis
#'
#' @param dataset data.frame, rows correspond to different genes, first column contains gene identifiers, second column contains the gene name, columns 3 to n contain the RNAseq count data, column names should start with the condition identifier (e.g. "WT") followed by the time separated by "_"
#' @param myanalyzeConditions Character vector, Name of experimental conditions
#' @param cores Numeric, Number of cores for parallelization, default 1 for no parallelization
#' @param mytimes Numeric vector, Time points of the time-resolved RNA-seq data
#'
#' @return Data.frame containing gene names, log fold change and p-values calculated from NBPSeq, each gene appears as often as available time points
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}

getFC <- function(dataset = mydata, myanalyzeConditions = analyzeConditions, cores = 1, mytimes = times){
  out <- do.call(rbind, mclapply(mytimes, function(t){
    data <- as.matrix(dataset[,c(which(
      grepl(paste0(myanalyzeConditions[1],"_",t,"_"), names(dataset)) | grepl(paste0(myanalyzeConditions[2],"_",t,"_"), names(dataset))
    ))])

    ## Specify treatment groups
    grp.ids = do.call(rbind,strsplit(colnames(data),"_"))[,1]  # Numbers or strings are both OK

    ## Estimate normalization factors
    norm.factors = estimate.norm.factors(data);

    ## Prepare an NBP object, adjust the library sizes by thinning the counts.
    set.seed(999);
    obj = prepare.nbp(data, grp.ids, lib.size=colSums(data), norm.factors=norm.factors, print.level = 0);

    ## Fit a dispersion model (NBQ by default)
    obj = estimate.disp(obj, print.level = 0);

    ## Perform exact NB test
    grp1 = analyzeConditions[2];
    grp2 = analyzeConditions[1];

    obj = exact.nb.test(obj, grp1, grp2, print.level = 0);

    # ## Output results
    out <- data.frame(name=attr(obj$log.fc, "names"), logFoldChange = obj$log.fc, pValue = obj$p.values, time=t)
    cbind(out, FCdetect=sapply(1:dim(out)[1], function(i){
      getD(out$pValue[i],
           out$logFoldChange[i],
           thFoldChange = 2,  ## ignored, if NA
           pValueFC = 0.05  ## ignored, if NA
      )
    }))
  }, mc.cores = cores))
  out$name <- rep(dataset$genename, length(times))
  return(out)
}

#' Auxiliary function for getFC()
#'
#' @param value output pValue from exact.nb.test
#' @param FC output logFoldChange from exact.nb.test
#' @param thFoldChange Numeric, threshold at which fold change is counted to be detected
#' @param pValueFC Numeric, p-value at which fold change is counted to be detected
#'
#' @return input plus decision whether fold change has been detected or not
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}


getD <- function(value, FC, thFoldChange=NA, pValueFC=0.05){
  if(is.na(value) | is.na(FC)){return(NA)} else {
    if(is.na(thFoldChange) & is.na(pValueFC)){
      return(NA)
    } else if(!is.na(thFoldChange) & is.na(pValueFC)){
      if(FC < -log2(thFoldChange)) {return(paste0(analyzeConditions[1], "<", analyzeConditions[2]))}
      else if(FC > log2(thFoldChange)) {return(paste0(analyzeConditions[1], ">", analyzeConditions[2]))}
      else return(paste0(analyzeConditions[1], "=", analyzeConditions[2]))
    } else if(is.na(thFoldChange) & !is.na(pValueFC)){
      if(value < pValueFC){
        if(FC < 0) {return(paste0(analyzeConditions[1], "<", analyzeConditions[2]))}
        if(FC > 0) {return(paste0(analyzeConditions[1], ">", analyzeConditions[2]))}
      } else return(paste0(analyzeConditions[1], "=", analyzeConditions[2]))
    } else {
      if((value < pValueFC) & (FC < -log2(thFoldChange))) {return(paste0(analyzeConditions[1], "<", analyzeConditions[2]))}
      else if((value < pValueFC) & (FC > log2(thFoldChange))) {return(paste0(analyzeConditions[1], ">", analyzeConditions[2]))}
      else return(paste0(analyzeConditions[1], "=", analyzeConditions[2]))
    }

  }
}


combineResults <- function(myresultSwitch = resultSwitch, myresultFC = resultFC){
  temp <- do.call(rbind, mclapply(myresultSwitch$genename, function(gene){
    getFCupdown(gene, myresultFC)
  }, mc.cores = nrcores))
  return(cbind(resultSwitch, temp))
}

plotSSGS <- function(myresultCombined = resultCombined, mytimes = times, myanalyzeConditions = analyzeConditions){
  out <- do.call(rbind, lapply(mytimes, function(t){
    do.call(rbind, lapply(paste0(format(seq(2.5,6,by=0.5), nsmall = 1),"hpf"), function(x){
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
    scale_fill_manual(name="Fisher Test Result", values=c("none" = "grey",
                                                          "1E-20 Enhancement" = "darkred", "1E-10 Enhancement" = "red", "1E-05 Enhancement"= "orange", "1E-02 Enhancement" = "yellow",
                                                          "1E-20 Suppression"= "violet", "1E-10 Suppression" = "darkblue","1E-05 Suppression" = "blue", "1E-02 Suppression" = "lightblue" ))
  return(P)
}

outputGeneTables <- function(myresultCombined = resultCombined, mytimes = times){
  genetableUp <- c()
  genetableDown <- c()
  geneNametableUp <- c()
  geneNametableDown <- c()
  for (t in mytimes){
    for(identifier in c("FCdown", "FCup")){
      for(x in paste0(format(seq(2.5,6,by=0.5), nsmall = 1),"hpf")){
        for(exp in analyzeConditions){

          if(identifier=="FCdown"){
            genesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCdown))$name))
            genesUp <- c(genesUp, rep("",3000-length(genesUp)))
            genetableUp <- cbind(genetableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(analyzeConditions[1]," > ",analyzeConditions[2]),x,genesUp))
            geneNamesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCdown))$genename))
            geneNamesUp <- c(geneNamesUp, rep("",3000-length(geneNamesUp)))
            geneNametableUp <- cbind(geneNametableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(analyzeConditions[1]," > ",analyzeConditions[2]),x,geneNamesUp))

            genesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCdown))$name))
            genesDown <- c(genesDown, rep("",3000-length(genesDown)))
            genetableDown <- cbind(genetableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(analyzeConditions[1]," > ",analyzeConditions[2]),x,genesDown))
            geneNamesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCdown))$genename))
            geneNamesDown <- c(geneNamesDown, rep("",3000-length(geneNamesDown)))
            geneNametableDown <- cbind(geneNametableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(analyzeConditions[1]," > ",analyzeConditions[2]),x,geneNamesDown))
          }
          if(identifier=="FCup"){
            genesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCup))$name))
            genesUp <- c(genesUp, rep("",3000-length(genesUp)))
            genetableUp <- cbind(genetableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(analyzeConditions[2]," > ",analyzeConditions[1]),x,genesUp))
            geneNamesUp <- as.character(unique(subset(myresultCombined, switch=="up" & timepoint==t & grepl(x, FCup))$genename))
            geneNamesUp <- c(geneNamesUp, rep("",3000-length(geneNamesUp)))
            geneNametableUp <- cbind(geneNametableUp,c(paste0("Switch Up at ", format(t, nsmall=1)),paste0(analyzeConditions[2]," > ",analyzeConditions[1]),x,geneNamesUp))

            genesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCup))$name))
            genesDown <- c(genesDown, rep("",3000-length(genesDown)))
            genetableDown <- cbind(genetableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(analyzeConditions[2]," > ",analyzeConditions[1]),x,genesDown))
            geneNamesDown <- as.character(unique(subset(myresultCombined, switch=="down" & timepoint==t & grepl(x, FCup))$genename))
            geneNamesDown <- c(geneNamesDown, rep("",3000-length(geneNamesDown)))
            geneNametableDown <- cbind(geneNametableDown,c(paste0("Switch Down at ", format(t, nsmall=1)),paste0(analyzeConditions[2]," > ",analyzeConditions[1]),x,geneNamesDown))
          }
        }
      }
    }
  }

  nametag <- analyzeConditions[2]

  write.table(genetableUp, file=paste0("genelist_switchUp_",nametag,".txt"), sep = "\t", row.names=F, col.names = F)
  write.table(genetableDown, file=paste0("genelist_switchDown_",nametag,".txt"), sep = "\t", row.names=F, col.names = F)
  write.table(geneNametableUp, file=paste0("geneNamelist_switchUp_",nametag,".txt"), sep = "\t", row.names=F, col.names = F)
  write.table(geneNametableDown, file=paste0("geneNamelist_switchDown_",nametag,".txt"), sep = "\t", row.names=F, col.names = F)


  resultOut1 <- subset(myresultCombined[,c("name","genename", "pvalueSwitch", "switch", "timepoint", "experiment")], experiment==analyzeConditions[1])
  resultOut2 <- subset(myresultCombined[,c("name","genename", "pvalueSwitch", "switch", "timepoint", "experiment")], experiment==analyzeConditions[2])
  resultOut1$timepoint <- format(resultOut1$timepoint, nsmall = 1)
  resultOut2$timepoint <- format(resultOut2$timepoint, nsmall = 1)
  if(dim(resultOut2)[1] > 0) resultOut <- cbind(resultOut1, resultOut2[,-1]) else resultOut <- resultOut1
  resultOut$switch <- sub("none", "0", resultOut$switch)
  resultOut$switch <- sub("up", "1", resultOut$switch)
  resultOut$switch <- sub("down", "-1", resultOut$switch)
  write.table(resultOut, file=paste0("switchList_",nametag,".txt"), sep = "\t", row.names=F)

}


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

