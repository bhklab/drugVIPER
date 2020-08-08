start<-function()
{
  # DrugVIPER manuscript code 
  #*************Libraries####################
  library("PharmacoGx")
  library(corrplot)
  library(Rarity)
  library(Biobase)
  library(calibrate)
  library(ggrepel)
  require(stringi)
  library(dplyr)
  # ************Parameters####
  flag_Fig2<-1
  flag_FigS3<-0
  flag_Fig2<-1
  flag_VIPER=1
  flag_dataType<-1
  flag_corr2mat<-0
  flag_Regression<-0
  flag_Sensetivity<-1
  flag_QINdata<-0
  flag_save=0
  flag_aracneAllGenes<-0
  paramDrug=""
  #paramCellline="HCC1569"
  n<-30
  #************gCSI##################
  biobase_gCSI <- readRDS("./Data/gCSI_gdrive.rds")
  rnaseqData <- PharmacoGx::summarizeMolecularProfiles(biobase_gCSI,mDataType = "rnaseq",fill.missing = F)
  eset_gCSI_common <- exprs(rnaseqData)#19142   146 
  ensCommon<-rownames(eset_gCSI_common)
  genes_mapping <- featureInfo(biobase_gCSI,"rnaseq")
  GeneNames <- genes_mapping[which(rownames(genes_mapping) %in% rownames(eset_gCSI_common)),"gene_name"]
  rownames(eset_gCSI_common)<-GeneNames
  source("fixCellLinesNames.R")
  colnames(eset_gCSI_common)<-fixCellLinesNames(colnames(eset_gCSI_common),"./cell_annotation_all.csv")
  if(flag_save==1){write.csv(eset_gCSI_common,"./op/gCSI_eset.csv")}
  
  
  #************PDX  ########################################################################################
  flag_TCGARegulons<-1
  library(RCurl) 
  library(Biostrings)
  library(rentrez)
  library(annotate)
  library(org.Hs.eg.db)
  library(annotate)
  library(Xeva)
  if(flag_TCGARegulons==1)
  {
    library(aracne.networks)
    data(package = "aracne.networks")$results[, "Item"]
    #data(regulonblca)
    # z1<-load(file ="./aracneAdj/regulons/skcm-tf-generegul.rda")
    # z2<-load(file ="./aracneAdj/regulons/skcm-sig-generegul.rda")
    z3 <- load(file = "./aracneAdj/regulons/skcm_allregulon.rda")
    regulonskcm <- regul
  }
  load("./Data/pdxe.rda")
  objXeva <- pdxe
  #objXeva@drug$treatment.target[1]<-"FGFR1,FGFR2,FGFR3"  esetXeva[DR5,] TF_pdxe_all[]
  dfDrugs <- read.csv("./data/pdxeDrugs_revised.csv")
  objXeva@drug$treatment.target <- dfDrugs$GeneTarget
  uniqueTissues <- unique(objXeva@molecularProfiles$RNASeq$tissue)
 
   if(!file.exists("./Data/TF_pdxe_all.csv"))
   {
     for (tID in 1:15)
     {
       print(tID)
       if (tID == 1)
       {
         reg <- regulonbrca
       } else if (tID == 2)
       {
         reg <- regulonread
       } else if (tID == 3)
       {
         reg <- regulonlusc
       } else if (tID == 4)
       {
         reg <- regulonkirc
       } else if (tID == 15)
       {
         reg <- regulonkirp
       } else if (tID == 5)
       {
         reg <- regulonpaad
       } else if (tID == 6)
       {
         reg <- regulonov
       } else if (tID == 7)
       {
         reg <- regulonsarc
       } else if (tID == 8)
       {
         reg <- regulonskcm
       } else if (tID == 9)
       {
         reg <- regulonucec
       } else if (tID == 10)
       {
         reg <- regulonlihc
       } else if(tID == 13)
       {
         reg<-regulonhnsc
       }
      
       if (!tID %in% c(11, 12, 14))
       {
         names(reg) <- lookUp(c(names(reg)), 'org.Hs.eg', 'SYMBOL')
         for (i in 1:length(names(reg)))
         {
           names(reg[[i]]$tfmode) <-
             lookUp(c(names(reg[[i]]$tfmode)), 'org.Hs.eg', 'SYMBOL')
         }
         # regucec<-reg# regbrca reghnsc reglihc reglusc regov regpaad regread regsarc reguloncesc
         if (!tID == 15)
         {
           tissueID <-
             which(objXeva@molecularProfiles$RNASeq$tissue %in% uniqueTissues[tID])
         } else
         {
           tissueID <-
             which(objXeva@molecularProfiles$RNASeq$tissue %in% uniqueTissues[4])
         }
         if (!tID %in% c(10, 13))
         {
           eset_tissue <-
             objXeva@molecularProfiles$RNASeq@assayData$exprs[, tissueID]
         }else
         {
           eset_tissue <-
             objXeva@molecularProfiles$RNASeq@assayData$exprs[, c(1, tissueID)]
         }
         TF_activities <-
           viper(eset_tissue, reg, verbose = TRUE, minsize = 1)
       }
       if(tID==1)
       {
         TF_BRCA<-TF_activities
       }else if (tID==2)
       {
         TF_read<-TF_activities
       }else if(tID==3)
       {
         TF_lusc<-TF_activities
       }else if(tID==4)
       {
         TF_kirc<-TF_activities
         
       }else if(tID==15)
       {
         TF_kirp<-TF_activities
       }else if(tID==5)
       {
         TF_paad<-TF_activities
       }else if(tID==6)
       {
         TF_ov<-TF_activities
       }else if(tID==7)
       {
         TF_sarc<-TF_activities
       }else if(tID==8)
       {
         TF_skcm<-TF_activities
       }else if(tID==9)
       {
         TF_ucec<-TF_activities
       }else if(tID==10)
       {
         TF_lihc<-TF_activities[,-1]#X-3769
       }else if(tID==13)
       {
         TF_hnsc<-TF_activities[,-1]#X-4819
       }
     }
     common <- intersect(rownames(TF_BRCA), rownames(TF_read))
     common <- intersect(common, rownames(TF_lusc))
     common <- intersect(common, rownames(TF_kirc))
     common <- intersect(common, rownames(TF_kirp))
     common <- intersect(common, rownames(TF_paad))
     common <- intersect(common, rownames(TF_ov))
     common <- intersect(common, rownames(TF_sarc))
     common <- intersect(common, rownames(TF_skcm))
     common <- intersect(common, rownames(TF_ucec))
     common <- intersect(common, names(TF_hnsc))
     common <- intersect(common, names(TF_lihc))
     TF_pdxe_all <-
       cbind(
         TF_BRCA[common, ],
         TF_read[common, ],
         TF_lusc[common, ],
         TF_kirc[common, ],
         TF_paad[common, ],
         TF_ov[common, ],
         TF_sarc[common, ],
         TF_skcm[common, ],
         TF_ucec[common, ],
         TF_lihc[common],
         TF_hnsc[common]
       )
     colnames(TF_pdxe_all)[395] <- "X-3769"
     colnames(TF_pdxe_all)[396] <- "X-4819"
     colsPDXE<-colnames(TF_pdxe_all)
     write.csv(TF_pdxe_all, "TF_pdxe_all.csv")
   }else
   {
     TF_pdxe_all <- read.csv("TF_pdxe_all.csv",check.names=FALSE)
     TF_pdxe_all<-as.matrix(TF_pdxe_all)
     rownames(TF_pdxe_all)<-TF_pdxe_all[,1]
     TF_pdxe_all<-TF_pdxe_all[,-1]
   }

  eset_pdxe<-objXeva@molecularProfiles$RNASeq@assayData$exprs[,colnames(TF_pdxe_all)]
  esetXeva<-eset_pdxe
  dfResXeva<-data.frame()
  
  for (i in 1: length(objXeva@drug$drug.id))
  {
    tryCatch({
    
    drugXeva<-objXeva@drug$drug.id[i]
    targetXeva<-unlist(strsplit(as.character(objXeva@drug$treatment.target[i]),","))
    sensetivityXeva<-Xeva::summarizeMolecularProfiles(object=objXeva, drug=drugXeva, mDataType="RNASeq")
    for(j in 1:length(targetXeva))
    {
      
      predictionXevaVector<-esetXeva[targetXeva[j],as.vector(sensetivityXeva$biobase.id)]
      predictionTFXevaVector<-as.double(TF_pdxe_all[targetXeva[j],as.vector(sensetivityXeva$biobase.id)])

      k<-1
      for (s in 1:k)
      {
        if(s==1)
        {
          sensetivityXevaVector<-sensetivityXeva$AUC
        }
        if(s==2)
        {
          sensetivityXevaVector<-sensetivityXeva$slope
        }
        if(s==3)
        {
          sensetivityXevaVector<-sensetivityXeva$best.average.response
        }
        wCI<-wCI::paired.concordance.index(predictions = predictionXevaVector,observations = sensetivityXevaVector,delta.obs = 0.12,delta.pred = 0,outx = FALSE)
        print(paste0(round(wCI$cindex,4)," ",round(wCI$p.value,4)))
        wCI_TF<-wCI::paired.concordance.index(predictions = predictionTFXevaVector,observations = sensetivityXevaVector,delta.obs = 0.12,delta.pred = 0,outx = FALSE)
        print(paste0(round(wCI_TF$cindex,4)," ",round(wCI_TF$p.value,4)))
        
        if(wCI_TF$cindex>0.5 | wCI$cindex>0.5)
        {
          dfResXeva<-rbind(dfResXeva,data.frame("drug"=drugXeva,"target"=targetXeva[j],"targetTF"=targetXeva[j],"ci"=round(wCI$cindex,4),"p"=round(wCI$p.value,4),"ciTF"=round(wCI_TF$cindex,4),"pTF"=round(wCI_TF$p.value,4)))
        }
      }
    
    }},error=function(e) {}, finally = {})
    
  }
  fdr<-p.adjust(dfResXeva$p)
  fdrTF<-p.adjust(dfResXeva$pTF)
  deltaCI<-dfResXeva$ciTF-dfResXeva$ci
  dfResXeva<-data.frame(dfResXeva,fdr,fdrTF,deltaCI)
  write.csv(dfResXeva,"dfResXeva13May.csv")
  boxplot(dfResXeva$ciTF, dfResXeva$ci,
          main = "CI comparision",
          #at = c(1,2,4,5),
          names = dfResXeva$drug,
          las = 2,
          col = c("orange","red"),
          border = "brown",
          horizontal = TRUE,
          notch = TRUE
  )
  
  dfResXevaDeltaCI<-aggregate(dfResXeva[, c(4,6)], list(dfResXeva$drug), median)
  dfResXevaDeltaCI<-data.frame(dfResXevaDeltaCI,"Delta"=dfResXevaDeltaCI$ciTF-dfResXevaDeltaCI$ci)
  dfResXevaDeltaCI <-dfResXevaDeltaCI[order(-dfResXevaDeltaCI$Delta),]
  
  
  alphaciTF<-c()
  for(i in dfResXevaDeltaCI$ciTF)
  {alphaciTF<-c(alphaciTF,(i-min(dfResXevaDeltaCI$ciTF))/(max(dfResXevaDeltaCI$ciTF)-min(dfResXevaDeltaCI$ciTF)))}
  alphaci<-c()
  for(i in dfResXevaDeltaCI$ci)
  {alphaci<-c(alphaci,(i-min(dfResXevaDeltaCI$ci))/(max(dfResXevaDeltaCI$ci)-min(dfResXevaDeltaCI$ci)))}
  dfResXevaDeltaCI<-data.frame(dfResXevaDeltaCI,"alphaci"=alphaci)
  dfResXevaDeltaCI<-data.frame(dfResXevaDeltaCI,"alphaciTF"=alphaciTF)

  significance<-c()
  for (i in dfResXevaDeltaCI$Group.1)
  {
    dfdrugSet<-subset(dfResXeva,dfResXeva$drug %in% i)
      obj1<-min(dfdrugSet$pTF)
      obj2<-min(dfdrugSet$p)
      significance[i]<-min(c(obj1,obj2))
  }
  
  library(RColorBrewer)
  coul = c("orange","purple")#brewer.pal(1, "Pastel2") 
  lbls<-dfResXevaDeltaCI$Group.1
  par(mar=c(9, 4.1,6.1, 1.1))
  plotdata<-dfResXevaDeltaCI$Delta
  mids<-barplot(plotdata,col=ifelse(dfResXevaDeltaCI$Delta>0,rgb(1, 0.3, 0, alpha=dfResXevaDeltaCI$alphaciTF),rgb(0.2, 0, 0.2, alpha=dfResXevaDeltaCI$alphaci)) ,beside = TRUE,ylab = "Delta CI",names.arg = lbls,las=2, xpd = FALSE,ylim=c(-0.2, 0.2),font = 2, yaxt = "n")
  axis(side = 2)
  legend("topright",
         c("VIPER","RNA-Seq"),
         fill = coul,cex = 0.8,border = FALSE)

  library(RColorBrewer)
  coul = c("orange","purple")#brewer.pal(1, "Pastel2") 
  lbls<-dfResXevaDeltaCI$Group.1
  par(mar=c(9, 4.1,6.1, 1.1))
  plotdata<-dfResXevaDeltaCI$Delta
  mids<-barplot(plotdata,col=ifelse(dfResXevaDeltaCI$Delta>0,coul[1],coul[2]) ,beside = TRUE,ylab = "Delta CI",names.arg = lbls,las=2, xpd = FALSE,ylim=c(-0.2, 0.2))
  text(x=mids , y= ifelse(dfResXevaDeltaCI$Delta>0,plotdata+0.01,plotdata-0.01), labels=ifelse(significance<0.01,'**',ifelse(significance<0.05,'*','')), xpd=FALSE)
 legend("topright",
          c("VIPER","rnaSeq"),
          fill = coul,cex = 0.8,border = FALSE)

  dfResXeva <-dfResXeva[order(-dfResXeva$ciTF),]
  CI<-c(dfResXeva$ciTF,dfResXeva$ci)
  DataType<-c(rep("VIPER",length(dfResXeva$ci)),rep("RNASeq",length(dfResXeva$ci)))
  Drug<-rep(dfResXeva$drug,2)
  p<-c(dfResXeva$pTF,dfResXeva$p)
  dfResXevaPlot<-data.frame("Drug"=Drug,"CI"=CI,"Type"=DataType,"p"=p)
  #######
  ggplot(dfResXevaPlot, aes(x=Drug, y=CI, fill=Type)) +  
    #geom_point(position = pd, size = 2) +
    geom_boxplot()+  
    scale_fill_manual(values=c("purple", "orange")) +
    geom_jitter(pch = 21,aes(shape = Type), position = position_jitterdodge())+
    geom_vline(xintercept = c(1.5:length(Drug)), color = "gray", size=0.1)+
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #######
  ggplot(dfResXevaPlot,aes(fill=dfResXevaPlot$Type,y=dfResXevaPlot$CI,x=dfResXevaPlot$Drug))+geom_bar(position="stack",stat="identity",alpha=0.5)+
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_manual(values=c("purple", "orange")) +
    geom_bar(position = "fill")+
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #######
  targets<-as.character(unique(dfResXeva$target))
  drugs<-unique(dfResXeva$drug)
  matrix_ci <- matrix(, nrow = length(drugs), ncol = length(targets))
  rownames(matrix_ci)<-drugs
  colnames(matrix_ci)<-targets
  for(i in 1:length(drugs)){
    for(j in 1:length(targets)){
      
      val<-dfResXeva[which(dfResXeva$drug==drugs[i] & dfResXeva$target==targets[j]& dfResXeva$p<0.05),"ci"]
      if(length(val)>0){
      matrix_ci[i, j] <-val }    }
  }
  library(circlize)
  chordDiagram(matrix_ci)
  #######
  targets<-as.character(unique(dfResXeva$target))
  drugs<-unique(dfResXeva$drug)
  matrix_ci <- matrix(, nrow = length(drugs), ncol = length(targets))
  rownames(matrix_ci)<-drugs
  colnames(matrix_ci)<-targets
  for(i in 1:length(drugs)){
    for(j in 1:length(targets)){
      
      val<-dfResXeva[which(dfResXeva$drug==drugs[i] & dfResXeva$target==targets[j]& dfResXeva$pTF<0.05),"ciTF"]
      if(length(val)>0){
        matrix_ci[i, j] <-val }    }
  }
  library(circlize)
  chordDiagram(matrix_ci)
  #######
  df<-as.matrix(df$drug,df$target)
  library(RColorBrewer)
  coul = brewer.pal(5, "Pastel2") 
  lbls<-paste0(df$Biomarker,"_",df$Drug,sep = "")
  par(mar=c(9, 4.1,6.1, 1.1))
  plotdata<-rbind(df$mci.RPPA,df$mci.VIPER,df$mci.RNASEQ,df$mci.CNV,df$mci.Mutation)
  fdrdata<-rbind(df$FDR.RPPA,df$FDR.VIPER,df$FDR.RNASEQ,df$FDR.CNV,df$FDR.Mutation)
  mids<-barplot(plotdata,col=coul ,beside = TRUE,ylab = "mCI",names.arg = paste(df$Biomarker,"-",df$Drug,sep = ""),las=2, xpd = FALSE,ylim=c(0.5, 1.0))
  text(x=mids , y= plotdata+0.01, labels=ifelse(fdrdata<0.01,'**',ifelse(fdrdata<0.05,'*','')), xpd=FALSE)
  legend("topright",
         c("CCLE-RPPA","CCLE-VIPER","rnaseq","CNV","Mutation"),
         fill = coul,cex = 0.5)
  
  ggplot(data =dfResXevaPlot,aes(x= Drug, y = CI, shape=group2, colour=group2,group = group2)) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl),colour = "black", width = 0.5, position = pd) + 
    geom_point(position = pd, size = 4) +  
    scale_colour_manual(values=c("red","blue","red","blue")) +   
    scale_shape_manual(values=c(19,19,17,17)) +
    scale_fill_hue(name="Treatment & State",  #Legend label, use darker colors
                   labels=c("Control", "Exclosure"),l=40)
  load("./Data/pdxe.rda")
  objXeva<-pdxe
  objXeva<-readRDS("./Data/PDX_TNBC_Nor_Xeva_Obj.rda")
  #function
  esetXeva<-objXeva@molecularProfiles$RNASeq@assayData$exprs#TN
  if(flag_save==1) {write.csv(esetXeva,"./op/esetXevaPDXE.csv")}#for VIPER
  if(length(rownames(esetXeva)>30000))
  {GeneNames <- genes_mapping[which(rownames(genes_mapping) %in% rownames(esetXeva)),"gene_name"]
  rownames(esetXeva)<-GeneNames}
  drugXeva<-"ERIBULIN"
  targetXeva<-"TUBB1"
  sensetivityXeva<-Xeva::summarizeMolecularProfiles(object=objXeva, drug=drugXeva, mDataType="RNASeq")
  predictionXevaVector<-esetXeva[targetXeva,sensetivityXeva$biobase.id]
  k<-1
  for (i in 1:k)
  {
    if(i==1)
    {
      sensetivityXevaVector<-sensetivityXeva$AUC
    }
    if(i==2)
    {
      sensetivityXevaVector<-sensetivityXeva$slope
    }
    if(i==3)
    {
      sensetivityXevaVector<-sensetivityXeva$best.average.response
    }
    wCI<-wCI::paired.concordance.index(predictions = predictionXevaVector,observations = sensetivityXevaVector,delta.obs = 0.12,delta.pred = 0,outx = FALSE)
    print(paste0(round(wCI$cindex,4)," ",round(wCI$p.value,4)))
  }
  id<-"xevaUHN_BRCAregulon"#"xevaUHN_BRCAregulon"  "xevaUHN"
  id_aracne<-"BRCA_TCGA"#"BRCA_TCGA" "xevaUHN"
  file_mrs_path<-paste("./op/",id,"","_mrs.csv",sep="",collapse = "")
  mrs<-as.matrix(read.csv(file_mrs_path, row.names=1, check.names=FALSE))
  file_aracne<-paste("./aracneAdj/",id_aracne,"_aracne",".txt",sep="",collapse = "")
  my_data<-read.table(file_aracne,header=TRUE,sep = '\t')
  target_mrs<-my_data[which(my_data$Target %in% targetXeva),]#"Regulator"
  max<--100
  max_mrs<-""
  for(i in 1:length(target_mrs$Regulator))
  {
    predictionXevaVector<-mrs[sensetivityXeva$biobase.id,target_mrs$Regulator[i]]
    wCI_mrs<-wCI::paired.concordance.index(predictions = predictionXevaVector,observations = sensetivityXevaVector,delta.obs = 0.12,delta.pred = 0,outx = FALSE)
    print(paste0(round(wCI_mrs$cindex,4)," ",round(wCI_mrs$p.value,4)," ",round(target_mrs$MI[i],4))) 
    if(max< wCI_mrs$cindex)
    {
      max<-wCI_mrs$cindex
      max_mrs<-target_mrs$Regulator[i]
    }
  }
  print(paste0(max_mrs," ",round(max,4)))
  print(paste0(round(wCI$cindex,4)," ",round(wCI$p.value,4)))
  
  # ************Load omics and sensitivity Pset (RNASeq, AAC) ####
  id<-"cclenoleukemia"#"ccle_single_mutation"
  load("./Data/CCLE_kallisto.RData")
  AAC <- t(PharmacoGx::summarizeSensitivityProfiles(pSet = CCLE,sensitivity.measure = "auc_recomputed",fill.missing = F)) #ic50_recomputed
  load("./Data/CTRPv2.RData")
  AAC_CTRPv2 <- t(PharmacoGx::summarizeSensitivityProfiles(pSet = CTRPv2,sensitivity.measure = "auc_recomputed",fill.missing = F)) #ic50_recomputed
  CRTPV2Classes<- read.csv("./Data/CRTPV2DrugClassesGDSC1000.csv")
  #if(flag_dataType==1) {
  samplesIDs <- rownames(CCLE@cell)[!is.na(CCLE@cell$tissueid) & CCLE@cell$tissueid != "haematopoietic_and_lymphoid_tissue"]
  dataType<-"CNV"
  cnvData <- PharmacoGx::summarizeMolecularProfiles(pSet = CCLE,mDataType = "cnv",fill.missing = F)#,summary.stat = "median"
  eset_CNV <- t(exprs(cnvData))# 742 24960
  geneSymbols<-colnames(eset_CNV)#fData(cnvData)[,"Symbol"]
  idxRemovedLeukemiaCelllines<-which(CCLE@cell$tissueid %in% "haematopoietic_and_lymphoid_tissue" & CCLE@cell$cellid %in% rownames(eset_CNV))#length(intersect(CCLE@cell$cellid,rownames(exprMatSymbols)))
  removedLeukemiaCelllines<-CCLE@cell$cellid[idxRemovedLeukemiaCelllines]#145 CL
  temp_celllines<-rownames(eset_CNV) [! rownames(eset_CNV) %in% removedLeukemiaCelllines]
  eset_CNV<-eset_CNV[temp_celllines,]#597 24960
  colnames(eset_CNV)<-geneSymbols
  #} else if(flag_dataType==2){
  dataType<-"mutation"
  mutationData <- PharmacoGx::summarizeMolecularProfiles(pSet = CCLE,mDataType = "mutation",fill.missing = F,summary.stat = "and")#or
  eset_Mutation <- t(exprs(mutationData))#1044 1667
  geneSymbols<-fData(mutationData)[,"Symbol"]
  idxRemovedLeukemiaCelllines<-which(CCLE@cell$tissueid %in% "haematopoietic_and_lymphoid_tissue" & CCLE@cell$cellid %in% rownames(eset_Mutation))#length(intersect(CCLE@cell$cellid,rownames(exprMatSymbols)))
  removedLeukemiaCelllines<-CCLE@cell$cellid[idxRemovedLeukemiaCelllines]#186 CL
  temp_celllines<-rownames(eset_Mutation) [! rownames(eset_Mutation) %in% removedLeukemiaCelllines]
  eset_Mutation<-eset_Mutation[temp_celllines,]#858 1667
  colnames(eset_Mutation)<-geneSymbols
  #} else if(flag_dataType==0){
  dataType<-"rnaseq"
  rnaseqData <- PharmacoGx::summarizeMolecularProfiles(pSet = CCLE,mDataType = "rnaseq",fill.missing = F)#,summary.stat = "median"#rnaseqData <- summarizeMolecularProfiles(pSet = CCLE,mDataType = "rnaseq",fill.missing = F)
  eset_rnaseq <- t(exprs(rnaseqData))#933 49922
  geneSymbols<-fData(rnaseqData)[,"Symbol"]
  idxRemovedLeukemiaCelllines<-which(CCLE@cell$tissueid %in% "haematopoietic_and_lymphoid_tissue" & CCLE@cell$cellid %in% rownames(eset_rnaseq))#length(intersect(CCLE@cell$cellid,rownames(exprMatSymbols)))
  removedLeukemiaCelllines<-CCLE@cell$cellid[idxRemovedLeukemiaCelllines]#171 CL
  temp_celllines<-rownames(eset_rnaseq) [! rownames(eset_rnaseq) %in% removedLeukemiaCelllines]
  eset_rnaseq<-eset_rnaseq[temp_celllines,]#
  colnames(eset_rnaseq)<-geneSymbols
  
  oCCLErnaseq<-eset_rnaseq #762 49922
  oCCLEcnv<-eset_CNV #597 24960
  oCCLEmutation<-eset_Mutation #858 1667
  
  eset<-eset_rnaseq
  if (flag_save==1)
  {
    write.csv(t(eset), file = file.path(getwd(),paste0("./aracneAdj/ccle_single_",dataType,"_expr4Aracne.csv")), quote =
                FALSE)
  }
  # Read ccleRPPA  ----------------------------------------------------------------
  ccleRPPA<- as.matrix(read.csv("./Data/ccleRPPA.csv", row.names=1,check.names=FALSE))#899 214 #row.names=1 #check.names=FLASE to remove x front from column name 
  ccleRPPATarget<- read.csv("./Data/ccleRPPAtarget.csv")
  ccleRPPATarget$Antibody_Name<-normalizeName(ccleRPPATarget$Antibody_Name)
  rownames(ccleRPPA)<-gsub("_.*","",rownames(ccleRPPA))
  source("fixCellLinesNames.R")
  rownames(ccleRPPA)<-fixCellLinesNames(rownames(ccleRPPA),"./cell_annotation_all.csv")
  colnames(ccleRPPA)<-normalizeName(colnames(ccleRPPA)) #setdiff(colnames(ccleRPPA),ccleRPPATarget$Antibody_Name)
  # Read MCLP ---------------------------------------------------------------
  mclpData<- as.matrix(read.csv("./Data/MCLP.csv", row.names=1,check.names=FALSE))#651 452
  mclpTarget<- read.csv("./Data/MCLP-Targets.csv")# read as dataframe
  mclpTarget$Antibody_Name<-normalizeName(mclpTarget$Antibody_Name)
  source("fixCellLinesNames.R")
  rownames(mclpData)<-fixCellLinesNames(rownames(mclpData),"./cell_annotation_all.csv")
  colnames(mclpData)<-normalizeName(colnames(mclpData))#setdiff(colnames(mclpData),mclpTarget$Antibody_Name)
  # Common cell lines--------------------
  wrong_CL_CCLERRPPA<-c("T173","ISHIKAWAHERAKLIO02ER","HS688AT","PECAPJ41CLONED2","KO52","MDAMB435S","SKNBE2","PECAPJ34CLONEC12","KPNSI9S","MDAPCA2B")
  replaced_CL_CCLERRPPA<-c("T1-73","Ishikawa (Heraklio) 02 ER-","Hs 688(A).T","PE/CA-PJ41 (clone D2)","KO52","MDA-MB-435","SK-N-BE(2)","PE/CA-PJ34 (clone C12)","KP-N-SI9s","MDA PCa 2b")
  rownames(ccleRPPA)[which(rownames(ccleRPPA)%in%wrong_CL_CCLERRPPA)]<-replaced_CL_CCLERRPPA
  idxRemovedLeukemiaCelllines<-which(CCLE@cell$tissueid %in% "haematopoietic_and_lymphoid_tissue" & CCLE@cell$cellid %in% rownames(ccleRPPA))#length(intersect(CCLE@cell$cellid,rownames(exprMatSymbols)))
  removedLeukemiaCelllines<-CCLE@cell$cellid[idxRemovedLeukemiaCelllines]#156 CL
  temp_celllines<-rownames(ccleRPPA) [! rownames(ccleRPPA) %in% removedLeukemiaCelllines]
  ccleRPPA<-ccleRPPA[temp_celllines,]
  oCCLERPPA<-ccleRPPA#741 214
  z<-setdiff(rownames(ccleRPPA),rownames(eset))
  wrong_CL_MCLP<-c()#("T173","ISHIKAWAHERAKLIO02ER","HS688AT","PECAPJ41CLONED2","KO52","MDAMB435S","SKNBE2","PECAPJ34CLONEC12","KPNSI9S","MDAPCA2B")
  replaced_CL_MCLP<-c()#c("T1-73","Ishikawa (Heraklio) 02 ER-","Hs 688(A).T","PE/CA-PJ41 (clone D2)","KO52","MDA-MB-435","SK-N-BE(2)","PE/CA-PJ34 (clone C12)","KP-N-SI9s","MDA PCa 2b")
  rownames(mclpData)[which(rownames(mclpData)%in%wrong_CL_MCLP)]<-replaced_CL_MCLP
  idxRemovedLeukemiaCelllines<-which(CCLE@cell$tissueid %in% "haematopoietic_and_lymphoid_tissue" & CCLE@cell$cellid %in% rownames(mclpData))#length(intersect(CCLE@cell$cellid,rownames(exprMatSymbols)))
  removedLeukemiaCelllines<-CCLE@cell$cellid[idxRemovedLeukemiaCelllines]#70 CL
  temp_celllines<-rownames(mclpData) [! rownames(mclpData) %in% removedLeukemiaCelllines]
  mclpData<-mclpData[temp_celllines,]
  oMCLPData<-mclpData#580 452
  
  #************SFig genomics/proteomics correlation####################
  if (flag_FigS3==1)
  {
    load("./Data/gCSI_kallisto.RData")
    AAC_gCSI <- t(PharmacoGx::summarizeSensitivityProfiles(pSet = gCSI,sensitivity.measure = "auc_recomputed",fill.missing = F)) #ic50_recomputed
    rnaseqData_gCSI<- PharmacoGx::summarizeMolecularProfiles(pSet = gCSI,mDataType = "rnaseq",fill.missing = F,summary.stat = "median")#rnaseqData <- summarizeMolecularProfiles(pSet = CCLE,mDataType = "rnaseq",fill.missing = F)
    eset_gCSI <- t(exprs(rnaseqData_gCSI)) #672 49922
    colnames(eset_gCSI)<-fData(rnaseqData_gCSI)[,"Symbol"]
    idxRemovedLeukemiaCelllines<-which(CCLE@cell$tissueid %in% "haematopoietic_and_lymphoid_tissue" & CCLE@cell$cellid %in% rownames(eset_gCSI))#length(intersect(CCLE@cell$cellid,rownames(exprMatSymbols)))
    removedLeukemiaCelllines<-CCLE@cell$cellid[idxRemovedLeukemiaCelllines]# 81 CL
    temp_celllines<-rownames(eset_gCSI) [! rownames(eset_gCSI) %in% removedLeukemiaCelllines]
    eset_gCSI<-eset_gCSI[temp_celllines,]
    ogCSI<-eset_gCSI #591 49922
    
    cl_gCSICCLE<-intersect(rownames(oCCLErnaseq),rownames(eset_gCSI))#382    not365
    cl_RPPAMCLP<-intersect(rownames(oMCLPData),rownames(oCCLERPPA))#294 not283
    cl_commonAll<-intersect(cl_gCSICCLE,cl_RPPAMCLP)#214
    esetcommonCL_cclerppa<-oCCLERPPA[cl_commonAll,]
    esetcommonCL_mclpdata<-oMCLPData[cl_commonAll,]
    esetcommonCL_CCLE<-oCCLErnaseq[cl_commonAll,]
    esetcommonCL_gcsi<-eset_gCSI[cl_commonAll,]
    eset_CCLE<-esetcommonCL_CCLE#eset
    eset_gCSI<-esetcommonCL_gcsi#eset_gCSI
    eset_MCLP<-esetcommonCL_mclpdata#mclpData
    eset_CCLERPPA<-esetcommonCL_cclerppa #ccleRPPA
    
    drugSavePath<-paste(dirMain,"./op/proteomics/","FigS1",sep="")
    
    saveFile<-paste(drugSavePath,"/RNASeqCCLE_VS_gCSI",sep="")
    df1<-correlate2Mat(eset_CCLE,eset_gCSI)#df1<-correlate2Mat(exprMatSymbols,exprMatSymbols2)
    df1<-data.frame(df1,"RNASeq (CCLE vs gCSI)")
    
    proteins_intersect_CCLE_MCLP<-intersect(colnames(eset_CCLERPPA),colnames(eset_MCLP))
    genes_intersect_CCLE_MCLP<-normalizeName( targetlst.get(ccleRPPATarget[which(ccleRPPATarget$Antibody_Name %in% proteins_intersect_CCLE_MCLP),"Target_Genes"]))#163
    genes_intersect_MCLP_CCLE<-normalizeName( targetlst.get(mclpTarget[which(mclpTarget$Antibody_Name %in% proteins_intersect_CCLE_MCLP),"Target_Genes"]))#171
    genes_MCLP_CCLE<-unique(c(genes_intersect_CCLE_MCLP,genes_intersect_MCLP_CCLE))
    genes_MCLP_CCLE<-intersect(genes_MCLP_CCLE,colnames(eset_CCLE))
    genes_MCLP_CCLE<-c(genes_MCLP_CCLE,"TIGAR","FOSL1","HSPA1B","RPS6KB1","NKX2-1")#length=178 genes - I used gene card to map the protein target genes to their genesymbol in CCLE-RNASeq
    df3<-correlate2Mat(eset_CCLE[,genes_MCLP_CCLE],eset_gCSI[,genes_MCLP_CCLE])
    df3<-data.frame(df3,"RNASeq (CCLE vs GCSI) of (CCLERPPA and MCLP) target genes")#"RNASeq (CCLE-GCSI) of (CCLERPPA-MCLP) target genes"
    
    genes_cclerppa<-normalizeName(unique(targetlst.get(ccleRPPATarget$Target_Genes)))
    genes_cclerppa<-intersect(genes_cclerppa,colnames(eset_CCLE))
    genes_cclerppa<-c(genes_cclerppa,"TIGAR")#length=173 genes
    df5<-correlate2Mat(eset_CCLE[,genes_cclerppa],eset_gCSI[,genes_cclerppa])
    df5<-data.frame(df5,"RNASeq (CCLE vs gCSI) of CCLERPPA target genes")
    
    genes_mclp<-unique(targetlst.get(mclpTarget$Target_Genes))#length=311
    genes_mclp<-intersect(genes_mclp,colnames(eset_CCLE))
    genes_mclp<-c(genes_mclp,c("HIST3H3","PTGDS","HSPB1","ATP5F1","ABL1","DPP4","FOSL1","H2AFX","HSPA1B","MT-CO2","RPS6KB1","PKM","RPA2","MAPT","TIGAR","ERCC4"))
    df6<-correlate2Mat(eset_CCLE[,genes_mclp],eset_gCSI[,genes_mclp])
    df6<-data.frame(df6,"RNASeq (CCLE vs gCSI) of MCLP target genes")

    saveFile<-paste(drugSavePath,"/mclp_VS_cclerppa",sep="")
    df2<-correlate2Mat.alias(eset_MCLP,eset_CCLERPPA,alias.prepare(data.frame("Antibody_Name"=ccleRPPATarget$Antibody_Name,"Target_Genes"=ccleRPPATarget$Antibody_Name)))
    df2<-data.frame(df2,"RPPA (CCLE vs MCLP)")
    if(flag_save==1){write.csv(df2,paste0(saveFile,".csv"))}
   
    dfplot1<-data.frame("id"=df1$id,"e"=df1$e,"p"=df1$p,"lbl"=df1$X.RNASeq..CCLE.vs.gCSI..)
    #dfplot4<-data.frame("id"=df4$id,"e"=df4$e,"p"=df4$p,"lbl"=df4$X.RNASeq.CCLE.vs.GCSI..for.RPPA.Target..)
    dfplot2<-data.frame("id"=df2$id,"e"=df2$e,"p"=df2$p,"lbl"=df2$X.RPPA..CCLE.vs.MCLP..)
    dfplot3<-data.frame("id"=df3$id,"e"=df3$e,"p"=df3$p,"lbl"=df3$X.RNASeq..CCLE.vs.GCSI..of..CCLERPPA.and.MCLP..target.genes.)
    dfplot5<-data.frame("id"=df5$id,"e"=df5$e,"p"=df5$p,"lbl"=df5$X.RNASeq..CCLE.vs.gCSI..of.CCLERPPA.target.genes.)
    dfplot6<-data.frame("id"=df6$id,"e"=df6$e,"p"=df6$p,"lbl"=df6$X.RNASeq..CCLE.vs.gCSI..of.MCLP.target.genes.)
    dfplot<-rbind(dfplot1,dfplot5,dfplot6,dfplot3,dfplot2)#df1=44792, df2=199, df3=176, df5=173, df6=307
    
    dfplot<-rbind(dfplot5,dfplot6,dfplot3,dfplot2)#df1=44792, df2=199, df3=176, df5=173, df6=307
    
    library(ggplot2)
    myplot<-ggplot(dfplot, aes(x = lbl, y = e, fill = lbl)) + 
      xlab("Data type")+ylab("Pearson Correlation")+
      geom_violin()+
      geom_boxplot(width = 0.2) 
    myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))    
    #wilcox.test(df1$e,df5$e)
    library(ggpubr)
    # Visualize: Specify the comparisons you want
    lst<-unique(dfplot$lbl)
    my_comparisons <- list( c(lst[1], lst[2]), c(lst[1], lst[3]), c(lst[1], lst[4]))
    p<-ggerrorplot(dfplot, x = "lbl", y = "e",
                desc_stat = "mean_sd", color = "lbl",
                add = "violin", add.params = list(color = "darkgray"))+ 
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ylab("Pearson Correlation")+xlab("Data type")+
      stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = 1.6)                  # Add global p-value
    p
    p + rremove("legend")+ rremove("x.text")+rremove("xlab")
  }
  
  # find common celllines between 
  eset<-oCCLErnaseq[intersect(rownames(oCCLErnaseq),rownames(oCCLERPPA)),]#702 49922
  ccleRPPA<-oCCLERPPA[intersect(rownames(oCCLErnaseq),rownames(oCCLERPPA)),]#702 214
  
  data_genes_proteomics_genomics<-intersect(targetlst.get(ccleRPPATarget$Target_Genes),colnames(eset))
  data_TFs<-as.matrix(read.csv("./Data/TFs.csv"))[,1]
  data_genes_study<-unique(c(data_genes_proteomics_genomics,data_TFs))
  #************ viperActivity ################################  
  library(viper)
  filedummyName<-""
  file_mrs_path<-paste("./op/",id,filedummyName,"_mrs.csv",sep="",collapse = "")
  file_mrs_gene_path<-paste("./op/",id,filedummyName,"_mrs_genes.csv",sep="",collapse = "")
  file_Session<-paste("./op/sessions/",id,filedummyName,".RData",sep="",collapse = "")
  file_mrs_regulonTarget_path<-paste("./op/",id,filedummyName,"_regulonTarget.csv",sep="",collapse = "")
  
  mrs<-as.matrix(read.csv(file_mrs_path, row.names=1, check.names=FALSE))
  mrsTargets<-as.matrix(read.csv(file_mrs_regulonTarget_path, row.names=1, check.names=FALSE)) 
  mrsTargetslst<-c()
  for(i in 1:length(mrsTargets))
  {
    mrsTargetslst<-c(mrsTargetslst,unlist(strsplit(mrsTargets[i], split=" ")))
  }
  mrsTargetslst<-unique(mrsTargetslst)
  if(flag_aracneAllGenes==1)
  {
    file_mrs_path_all<-"./op/cclenoleukemiaallgenes_mrs.csv"
    file_mrs_regulonTarget_path_all<-"./op/cclenoleukemiaallgenes_regulonTarget.csv"
  }
 
  mrs_all<-as.matrix(read.csv(file_mrs_path_all, row.names=1, check.names=FALSE))
  mrsTargets_all<-as.matrix(read.csv(file_mrs_regulonTarget_path_all,  check.names=FALSE))
  mrsTargetslst_all<-c()
  for(i in 1:length(mrsTargets_all))
  {
    mrsTargetslst_all<-c(mrsTargetslst_all,unlist(strsplit(mrsTargets_all[i], split=" ")))
  }
  mrsTargetslst_all<-unique(mrsTargetslst_all)
  
  #************SFig genomics/proteomics correlation########################################################################################
  if (flag_Fig2==1)
  { 
    # Correlation between. 1) Gene expression vs Proteomics 2) VIPER and Proteomics (CCLE-RPPA) 3) VIPER and Gene expression
    drugSavePath<-paste(dirMain,"./op/proteomics/","Fig2",sep="")
    ########### RNASeq-exprMat Vs CCLE-RPPA ########### 
    saveFile<-paste(drugSavePath,"/RNASeq Vs CCLE-RPPA",sep="")
    df1<-correlate2Mat.alias(eset,ccleRPPA,alias.prepare(ccleRPPATarget))
    df1<-data.frame(df1,"RNASeq Vs CCLE-RPPA")
    if(flag_save==1){write.csv(df1,paste0(saveFile,".csv"))}
    #z<-length(intersect(unique(targetlst.get(ccleRPPATarget$Target_Genes)),colnames(exprMatSymbols)))#226 data points based on 172 unique genes but 
    ###########  VIPER(all genes) Vs RNASeq ########### 
    saveFile<-paste(drugSavePath,"/VIPER_all Vs RNASeq",sep="")
    df5<-correlate2Mat.alias(eset,mrs_all,alias.prepare(data.frame("Antibody_Name"=colnames(mrs_all),"Target_Genes"=colnames(mrs_all))))
    df5<-data.frame(df5,"VIPER(all genes) Vs RNASeq")
    if(flag_save==1){write.csv(df5,paste0(saveFile,".csv"))}
    #z<-length(intersect(mrsTargetslst,colnames(exprMatSymbols)))#39117
    ###########  VIPER Vs RNASeq ########### 
    saveFile<-paste(drugSavePath,"/VIPER Vs RNASeq",sep="")
    df2<-correlate2Mat.alias(eset,mrs,alias.prepare(data.frame("Antibody_Name"=colnames(mrs),"Target_Genes"=colnames(mrs))))
    df2<-data.frame(df2,"VIPER(TF) Vs RNASeq")
    if(flag_save==1){write.csv(df2,paste0(saveFile,".csv"))}
    #z<-length(intersect(mrsTargetslst,colnames(exprMatSymbols)))#39117
    ########### VIPER vs CCLE-RPPA ########### 
    saveFile<-paste(drugSavePath,"/VIPER Vs CCLE-RPPA",sep="")
    #df3<-correlate2Mat.alias(mrs,ccleRPPA,alias.prepare(data.frame("Antibody_Name"=ccleRPPATarget$Antibody_Name,"Target_Genes"=ccleRPPATarget$Target_Genes)))
    commonTargets<-intersect(unique(mrsTargetslst),unique(targetlst.get(ccleRPPATarget$Target_Genes)))
    aracneTarget<-read.table(file_aracne,header=TRUE,sep = '\t')
    df3<-correlate2Mat.RPPA_VIPER(mrs,ccleRPPA,aracneTarget,ccleRPPATarget,commonTargets)
    df3<-data.frame(df3,"VIPER Vs CCLE-RPPA target genes")
    if(flag_save==1){write.csv(df3,paste0(saveFile,".csv"))}
    z<-length(intersect(unique(mrsTargetslst),unique(targetlst.get(ccleRPPATarget$Target_Genes))))#
    ########### VIPER vs CCLE-RPPA ########### 
    saveFile<-paste(drugSavePath,"/VIPER Vs CCLE-RPPA",sep="")
    df4<-correlate2Mat.alias(mrs,ccleRPPA,alias.prepare(data.frame("Antibody_Name"=ccleRPPATarget$Antibody_Name,"Target_Genes"=ccleRPPATarget$Target_Genes)))
    df4<-data.frame(df4,"VIPER(TF) Vs CCLE-RPPA")
    if(flag_save==1){write.csv(df4,paste0(saveFile,".csv"))}
    
    dfplot1<-data.frame("id"=df1$id,"e"=df1$e,"p"=df1$p,"lbl"=df1$X.RNASeq.Vs.CCLE.RPPA.)
    dfplot5<-data.frame("id"=df5$id,"e"=df5$e,"p"=df5$p,"lbl"=df5$X.VIPER.all.genes..Vs.RNASeq.)
    dfplot2<-data.frame("id"=df2$id,"e"=df2$e,"p"=df2$p,"lbl"=df2$X.VIPER.TF..Vs.RNASeq.)
    dfplot3<-data.frame("id"=df3$id,"e"=df3$e,"p"=df3$p,"lbl"=df3$X.VIPER.Vs.CCLE.RPPA.target.genes.)
    dfplot4<-data.frame("id"=df4$id,"e"=df4$e,"p"=df4$p,"lbl"=df4$X.VIPER.TF..Vs.CCLE.RPPA.)
    
    dfplot<-rbind(dfplot5,dfplot2,dfplot1,dfplot4,dfplot3)
    
    #violin-plot
    library(ggplot2)
    myplot<-ggplot(dfplot, aes(x = lbl, y = e, fill = lbl)) + 
      xlab("Data type")+ylab("Pearson Correlation")+
      geom_violin()+
      geom_boxplot(width = 0.2) 
    
    myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(text = element_text(size=12))
    
    #wilcox.test
    library(ggpubr)
    lst<-unique(dfplot$lbl)
    my_comparisons <- list( c(lst[1], lst[2]), c(lst[1], lst[3]), c(lst[1], lst[4]) , c(lst[1], lst[5]))
    #--
    ggballoonplot(dfplot, x = "lbl", y = "e",
                  color = "lbl", palette = "jco")+ theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ylab("Pearson Correlation")+xlab("Data type")+
      stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = 1.8)    # Add global p-value
    #--
    p<-ggerrorplot(dfplot, x = "lbl", y = "e",
                desc_stat = "mean_sd", color = "lbl",
                add = "violin", add.params = list(color = "darkgray"))+ 
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ylab("Pearson Correlation")+xlab("Data type")+
      stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = 1.8)                  # Add global p-value
    p
    p + rremove("legend")+ rremove("x.text")+rremove("xlab")
  }
  mrsFig2<-mrs
  mrsTargetsFig2<-mrsTargets
  mrsTargetslstFig2<-mrsTargetslst
  # find common celllines between 
  cl_common<-intersect(rownames(oCCLErnaseq),rownames(oCCLEcnv))
  cl_common<-intersect(cl_common,rownames(oCCLEmutation))
  cl_common<-intersect(cl_common,rownames(oCCLERPPA))
  eset<-oCCLErnaseq[cl_common,]
  ccleRPPA<-oCCLERPPA[cl_common,]
  eset_CNV<-oCCLEcnv[cl_common,]
  eset_Mutation<-oCCLEmutation[cl_common,]
  file_mrs_path<-paste("./op/",id,filedummyName,"_omics_mrs.csv",sep="",collapse = "")
  file_mrs_gene_path<-paste("./op/",id,filedummyName,"_omics_mrs_genes.csv",sep="",collapse = "")
  mrs<-as.matrix(read.csv(file_mrs_path, row.names=1, check.names=FALSE))
  mrsTargets<-as.matrix(read.csv(file_mrs_gene_path, row.names=1, check.names=FALSE)) 
  
  mrsTargetslst<-c()
  for(i in 1:length(mrsTargets))
  {
    mrsTargetslst<-c(mrsTargetslst,unlist(strsplit(mrsTargets[i], split=" ")))
  }
  mrsTargetslst<-unique(mrsTargetslst)
  
  
  #************Fig forest & 4 barplot of validation biomarkers##############################################
  Biomark4valid<- read.csv("./Data/Ben-biomarkers-validation.csv", row.names=1)
  
  biomarkers.uniq<-unique(Biomark4valid$gene)                                #NQO1    BRAF    ALK     EGFR    HGF     MET     MDM2    ERBB2   FLT3  BCR_ABL
  biomarkers.mrTF<-intersect(Biomark4valid$gene,mrsTargetslst)              #"NQO1"  "BRAF"  "ALK"   "EGFR"  "HGF"   "MET"   "MDM2"  "ERBB2" "FLT3"
  biomarkers.rppa<-intersect(Biomark4valid$gene,ccleRPPATarget$Target_Genes)#         "BRAF"         "EGFR"          "MET"   "MDM2"  "ERBB2"
  biomarkers.mclp<-intersect(Biomark4valid$gene,mclpTarget$Target_Genes)              #"BRAF"        "EGFR"          "MET"   "MDM2"  "ERBB2" "FLT3" 
  
  exprMatSymbols<-eset_Mutation#eset
  uniqueCompounds<-unique(Biomark4valid$compound)
  flag_save<-1
  drugSavePath<-paste(dirMain,"./op/proteomics/",Sys.Date(),sep="")
  if (!file.exists(drugSavePath)){dir.create(drugSavePath)}
  for (drugindx in 1: length(uniqueCompounds)) 
  {
    ##### Hassan#### drugindx<-11  #comment
    print(drugindx)
    paramDrug<-as.character(uniqueCompounds[drugindx])
    validationSet<-unique(Biomark4valid[which(Biomark4valid[,"compound"]==paramDrug),"gene"])
    
    drugSavePath<-paste(dirMain,"./op/proteomics/",Sys.Date(),"/", paramDrug,sep="")
    if (!file.exists(drugSavePath)){dir.create(drugSavePath)}
    
    if(flag_Sensetivity==1)
    {
      ########### 1A- AAC vs exprs ########### 
      saveFile<-paste(drugSavePath,"/1a-AAC_VS_exprs","_RNASeq",sep="")
      df<-correlation.mci.AAC(AAC,eset,paramDrug)#df<-data.frame("e"=df$e,"p"=df$p,"id"=df$id,"lbl"=df$lbl,"p_fdr"=df$p_fdr)
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      plotVolc.ggrepel(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
     
 #############################################################################
      
      ########### 1B- AAC vs CNV ########### 
      saveFile<-paste(drugSavePath,"/1b-AAC_VS_exprs","_CNV",sep="")
      df<-correlation.mci.AAC(AAC,eset_CNV,paramDrug)#df<-data.frame("e"=df$e,"p"=df$p,"id"=df$id,"lbl"=df$lbl,"p_fdr"=df$p_fdr)
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
      #plot.df(df,0.5,validationSet,"mCI Gene-expression and AAC",n)
      
      ########### 1C- AAC vs Mutation ########### 
      saveFile<-paste(drugSavePath,"/1c-AAC_VS_exprs","_Mutation",sep="")
      df<-correlation.mci.AAC(AAC,eset_Mutation,paramDrug)#df<-data.frame("e"=df$e,"p"=df$p,"id"=df$id,"lbl"=df$lbl,"p_fdr"=df$p_fdr)
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
      #plot.df(df,0.5,validationSet,"mCI Gene-expression and AAC",n)
      
      #}
      #}  
      
      ########### 2- AAC vs ccleRPPA ########### 
      saveFile<-paste(drugSavePath,"/2-AAC_VS_ccleRPPA",sep="")
      df<-correlation.mci.AAC(AAC,ccleRPPA,paramDrug)
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0.5,validationSet,"mCI CCLE-RPPA,AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"mCI CCLE-RPPA,AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
      
      df<-correlation.mci.AAC.targets(AAC,ccleRPPA,paramDrug,alias.prepare(ccleRPPATarget))
      if(flag_save==1){write.csv(df,paste0(saveFile,"_targets",".csv"))}
      ##plotdf(df,0.5,validationSet,"mCI_Targets CCLE-RPPA,AAC",n,paste0(saveFile,"_target_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"mCI_Targets CCLE-RPPA,AAC",n,paste0(saveFile,"_target_volc.jpg",sep=""))
      
      #dfforest<-subset(df, id %in% validationSet)
      #plotforest(subset(df, id %in% validationSet),"Title")#/+++++++++++++++++++++++++++++11/3/19
      #commonCCLERPPAgenes<-intersect(fData(rnaseqData)[,"Symbol"],ccleRPPATarget$Antibody_Name)
      ########### 3- AAC vs MCLP ########### 
      saveFile<-paste(drugSavePath,"/3-AAC_VS_MCLP",sep="")
      df<-correlation.mci.AAC(AAC,mclpData,paramDrug)
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0.5,validationSet,"mCI MCLP-RPPA,AAC",n,paste0(saveFile,"_plot.tiff",sep="")) 
      plotVolc.ggrepel(df,0.5,validationSet,"mCI MCLP-RPPA,AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
      
      df<-correlation.mci.AAC.targets(AAC,mclpData,paramDrug,alias.prepare(mclpTarget))
      if(flag_save==1){write.csv(df,paste0(saveFile,"_targets",".csv"))}
      plotVolc.ggrepel(df,0.5,validationSet,"mCI_Targets MCLP-RPPA,AAC",n,paste0(saveFile,"_volc_plot.jpg",sep=""))
      ##plotdf(df,0.5,validationSet,"mCI_Targets MCLP-RPPA,AAC",n,paste0(saveFile,"_target_plot.tiff",sep=""))  
    }
    ########### 4- AAC vs mrs ########### 
    flag_VIPER<-1
    if(flag_VIPER==1)
    {
      ##hassan only for mrs #########    validationSet<-validationSet[2]    #####important############
      #drugindx<-6 drugindx<-11
      # validationSet<-validationSet[2]
      for(m in 1:length(validationSet))
      {
        interestSet<-c()
        for(i in 1:dim(mrsTargets)[1])
        {
          validIntersect<-intersect(unlist(strsplit(as.character(mrsTargets[i,1]), split=" ")),validationSet[m])
          if (length(validIntersect>0))
          {
            interestSet<-c(interestSet,rownames(mrsTargets)[i])
            print(rownames(mrsTargets)[i])
          }
        }
        
        saveFile<-paste(drugSavePath,"/4-AAC_VS_mrs",sep="")
        df<-correlation.mci.AAC(AAC,mrs,paramDrug)
        if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
        ##plotdf(df,0.5,validationSet,"mCI mrs_TF,AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
        plotVolc.ggrepel(df,0.5,validationSet[m],"mCI mrs_TF,AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
        
        for (i in 1: length(df$lbl)) 
        {
          intersectSet<-intersect(interestSet,df$lbl[i])
          if(length(intersectSet)==0)
          {
            df$lbl[i]<-""
          } else{  
          }
        }
        if(flag_save==1){write.csv(df,paste0(saveFile,"_interesting_targets_",validationSet[m],".csv"))}
      }
    }
    if(flag_corr2mat==1)
    {
      ########### 5-rnaSeq-exprMat vs MCLP########### 
      saveFile<-paste(drugSavePath,"/5-expr_VS_mclp",sep="")
      df<-correlate2Mat.alias(exprMatSymbols,mclpData,alias.prepare(mclpTarget))#  df<-correlate2Mat(exprMatSymbols,mclpData)
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      plotVolc.ggrepel(df,0.5,validationSet,"Correlation MCLP,Gene-expression",n,paste0(saveFile,"_volc.jpg",sep=""))
      ##plotdf(df,0,validationSet,"Correlation MCLP,Gene-expression",n,paste0(saveFile,"_plot.tiff",sep=""))
      
      #df<-correlate2MatmCI.alias(exprMatSymbols,mclpData,colnames(mclpData))
      #plotdf(df,0,validationSet,"mCI MCLP,Gene-expression")
      #plot.df(df,0,validationSet,"mCI MCLP,Gene-expression")
      ########### 6- rnaSeq-exprMat vs RPPA-CCLE ########### 
      saveFile<-paste(drugSavePath,"/6-expr_VS_cclerppa",sep="")
      df<-correlate2Mat.alias(exprMatSymbols,ccleRPPA,alias.prepare(ccleRPPATarget))
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0,validationSet,"Correlation RPPA,Gene-expression",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"Correlation RPPA,Gene-expression",n,paste0(saveFile,"_volc.jpg",sep=""))
      
      #df<-correlate2MatmCI.alias(exprMatSymbols,ccleRPPA,alias.prepare(ccleRPPATarget))
      #plotdf(df,0,validationSet,"mCI RPPA,Gene-expression")
      #plot.df(df,0,validationSet,"mCI RPPA,Gene-expression")  
      ########### 7- rnaSeq-exprMat vs mrs ########### 
      saveFile<-paste(drugSavePath,"/7-expr_VS_mrs",sep="")
      df<-correlate2Mat.alias(exprMatSymbols,mrs,alias.prepare(data.frame("Antibody_Name"=colnames(mrs),"Target_Genes"=colnames(mrs))))
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0,validationSet,"Correlation mrs,Gene-expression",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"Correlation mrs,Gene-expression",n,paste0(saveFile,"_volc.jpg",sep=""))
      
      #df<-correlate2MatmCI.alias(exprMatSymbols,mrs,colnames(mrs))
      #plotdf(df,0,validationSet,"mCI mrs,Gene-expression")
      #plot.df(df,0,validationSet,"mCI mrs,Gene-expression")  
      ########### 8- MCLP vs RPPA-CCLE ########### 
      saveFile<-paste(drugSavePath,"/8-mclp_VS_cclerppa",sep="")
      df<-correlate2Mat.alias(mclpData,ccleRPPA,alias.prepare(data.frame("Antibody_Name"=ccleRPPATarget$Antibody_Name,"Target_Genes"=ccleRPPATarget$Antibody_Name)))
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0,validationSet,"Correlation MCLP,RPPA",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"Correlation MCLP,RPPA",n,paste0(saveFile,"_volc.jpg",sep=""))
      
      # Create R ggplot Violin plot
      # Importing the ggplot2 library
      # library(ggplot2)
      # # Create a Violin plot
      # ggplot(diamonds, aes(x = cut, y = price)) + 
      #   geom_violin()
      # 
      # ggplot(df, aes(x ="RPPA Vs MCLP" , y = e)) + 
      #   geom_violin()
      
      ###########9- mrs vs RPPA-CCLE ########### 
      saveFile<-paste(drugSavePath,"/9-mrs_VS_cclerppa",sep="")
      #df<-correlate2Mat.alias(mrs,ccleRPPA,alias.prepare(data.frame("Antibody_Name"=ccleRPPATarget$Antibody_Name,"Target_Genes"=ccleRPPATarget$Antibody_Name)))
      df<-correlate2Mat.alias(mrs,ccleRPPA,alias.prepare(data.frame("Antibody_Name"=ccleRPPATarget$Antibody_Name,"Target_Genes"=ccleRPPATarget$Target_Genes)))
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0,validationSet,"Correlation mrs,RPPA",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"Correlation mrs,RPPA",n,paste0(saveFile,"_volc.jpg",sep=""))
      
      ###########10- mrs vs MCLP ########### 
      saveFile<-paste(drugSavePath,"/10-mrs_VS_mclp",sep="")
      #df<-correlate2Mat.alias(mrs,mclpData,alias.prepare(data.frame("Antibody_Name"=mclpTarget$Antibody_Name,"Target_Genes"=mclpTarget$Antibody_Name)))
      df<-correlate2Mat(mrs,mclpData,alias.prepare(data.frame("Antibody_Name"=mclpTarget$Antibody_Name,"Target_Genes"=mclpTarget$Target_Genes)))
      if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
      ##plotdf(df,0,validationSet,"Correlation mr,mclp",n,paste0(saveFile,"_plot.tiff",sep=""))
      plotVolc.ggrepel(df,0.5,validationSet,"Correlation mrs,MCLP",n,paste0(saveFile,"_volc.jpg",sep=""))
    }
  }
  #************Fig The Elastic net model results based on combining omics for known biomarkers####
  Fig6<-1 
  drugindxxx<-c(2,4,5,6,7,9,12,13)
  if(Fig6==1)
  {
    
    viperFixedBiomarkers<-c("NR4A2","TSC22D1" ,"TSC2DD1","MTERF2","ZBED2","AHR","ZXDA","ZNF800","GRHL1","TMF1","CPXCR1","ZBED2","ATOH8")
    # Read validation biomarkers ----------------------------------------------
    Biomark4valid<- read.csv("./Data/Ben-biomarkers-validation.csv", row.names=1)
    biomarkers<-Biomark4valid[-c(6,10,12,13),]
    biomarkers<-cbind(biomarkers,as.data.frame(viperFixedBiomarkers))
    #uniqueCompounds<-unique(Biomark4valid$compound)
    flag_save<-1
    flagAvgCV<-0
    #if (!file.exists(drugSavePath)){dir.create(drugSavePath)}
    for (drugindx in 1: length(biomarkers)) 
    {
      ##### Hassan#### drugindx<-11  #comment
      #paramDrug<-as.character(uniqueCompounds[drugindx])
      paramDrug<-as.character(biomarkers$compound[drugindx])
      validationSet<-as.character(biomarkers$gene[drugindx])
      viperbiomarker<-as.character(biomarkers$viperFixedBiomarkers[drugindx])
      ##################      
      A<-AAC
      B<-eset
      commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
      x_exp=B[commonSample,validationSet,drop=F]
      y_exp = A[commonSample,paramDrug,drop=F]
      colnames(x_exp)<-paste0(colnames(x_exp),"_exp")
     ##################      
      A<-AAC
      B<-eset_CNV
      commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
      x_cnv=B[commonSample,validationSet,drop=F]
      y_cnv = A[commonSample,paramDrug,drop=F]
      colnames(x_cnv)<-paste0(colnames(x_cnv),"_cnv")
      ###############
      A<-AAC
      B<-eset_Mutation
      commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
      x_mutation=B[commonSample,validationSet,drop=F]
      y_mutation = A[commonSample,paramDrug,drop=F]
      colnames(x_mutation)<-paste0(colnames(x_mutation),"_mutation")
      x_mutation<-as.factor(x_mutation)
      ##############
      A<-AAC
      B<-ccleRPPA
      commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
      alias<-alias.prepare(ccleRPPATarget)
      for (i in 1:dim(B)[2])
      {
        tryCatch({target <-as.character(alias[colnames(B)[i],1])},error=function(e) {target<-colnames(B)[i]}, finally = {})
        target<-unlist(strsplit(target, split=" "))
        for(j in 1:length(target))#more than one target
        {
          if(target %in% validationSet)
          {
            x_cclerppa<-cbind(B[commonSample,colnames(B)[i]])
            colnames(x_cclerppa)<-append(colnames(x_cclerppa),paste0(target,"_(",colnames(B)[i],"_RPPA)"))
            
          }
        }
      }
      y_cclerppa = A[commonSample,paramDrug,drop=F]
     ###############
      A<-AAC
      B<-mrs
      commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
      x_mrs=B[commonSample,viperbiomarker,drop=F]
      y_mrs = A[commonSample,paramDrug,drop=F]
      
      y<-y_exp
      
      rci<-c()
      p<-c()
      l<-c()
      u<-c()
      corr<-c()
      
      
      x<-x_exp#
      x = cbind(x, 0)
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-x_cnv#
      x = cbind(x, 0)
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-x_mutation#
      x = cbind(x, 0)
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-x_cclerppa#
      x = cbind(x, 0)
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-x_mrs#
      x = cbind(x, 0)
      df<-elasticnet_RCI(x,y)#0.9022299
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv)# 0.7989454
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_mutation)# 0.8225363
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cclerppa)#0.8025581
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_mrs)#0.773066
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_mutation)#0.7578536
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_cclerppa)#0.8586757
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_mrs)#0.7743463 
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_mutation,x_cclerppa)#0.8471523 
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_mutation,x_mrs)#0.5
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cclerppa,x_mrs)#0.850145
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_mutation)#0.8704798
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_cclerppa)#0.8465934
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_mrs)#0.7994639
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_mutation,x_cclerppa)#0.8753646
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_mutation,x_mrs)# 0.8427383
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cclerppa,x_mrs)#0.8048257
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_mutation,x_cclerppa)#0.773066
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_mutation,x_mrs)#0.7578536
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_cclerppa,x_mrs)#0.8638557
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_mutation,x_cclerppa,x_mrs)#0.8602242
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_mutation,x_cclerppa)# 0.8748432
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_mutation,x_mrs)#0.8704798
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_cclerppa,x_mrs)#0.852375
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_mutation,x_cclerppa,x_mrs)#0.8753646
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_cnv,x_mutation,x_cclerppa,x_mrs)#0.8345003
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      x<-cbind(x_exp,x_cnv,x_mutation,x_cclerppa,x_mrs)#0.8748432
      df<-elasticnet_RCI(x,y)
      if(flagAvgCV==1)
      {
        rci<-append(rci,mean(df$rCI, na.rm=TRUE))#0.773066
        p<-append(p,mean(df$p, na.rm=TRUE))
        l<-append(l,mean(df$l, na.rm=TRUE))
        u<-append(u,mean(df$u, na.rm=TRUE))
        corr<-append(corr,mean(df$corr, na.rm=TRUE))
      }
      else
      {
        rci<-append(rci,df$rCI)#0.773066
        p<-append(p,df$p)
        l<-append(l,df$l)
        u<-append(u,df$u)
        corr<-append(corr,df$corr)
      }
      
      rnames<-c("RNASeq","CNV","Mutation","RPPA","VIPER",
                "RNASeq,CNV",
                "RNASeq,Mutation",
                "RNASeq,RPPA",
                "RNASeq,VIPER",
                "CNV,Mutation",
                "CNV,RPPA",
                "CNV,VIPER",
                "Mutation,RPPA",
                "Mutation,VIPER",
                "RPPA,VIPER",
                "RNASeq,CNV,Mutation",
                "RNASeq,CNV,RPPA",
                "RNASeq,CNV,VIPER",
                "RNASeq,Mutation,RPPA",
                "RNASeq,Mutation,VIPER",
                "RNASeq,RPPA,VIPER",
                "CNV,Mutation,RPPA",
                "CNV,Mutation,VIPER",
                "CNV,RPPA,VIPER",
                "Mutation,RPPA,VIPER",
                "RNASeq,CNV,Mutation,RPPA",
                "RNASeq,CNV,Mutation,VIPER",
                "RNASeq,CNV,RPPA,VIPER",
                "RNASeq,Mutation,RPPA,VIPER",
                "CNV,Mutation,RPPA,VIPER",
                "RNASeq,CNV,Mutation,RPPA,VIPER")
      guide<-rbind(c(0,0,1,0,0),
               c(0,0,0,1,0),
               c(0,0,0,0,1),
               c(1,0,0,0,0),
               c(0,1,0,0,0),
               c(0,0,1,1,0),
               c(0,0,1,0,1),
               c(1,0,1,0,0),
               c(0,1,1,0,0),
               c(0,0,0,1,1),
               c(1,0,0,1,0),
               c(0,1,0,1,0),
               c(1,0,0,0,1),
               c(0,1,0,0,1),
               c(1,1,0,0,0),
               c(0,0,1,1,1),
               c(1,0,1,1,0),
               c(0,1,1,1,0),
               c(1,0,1,0,1),
               c(0,1,1,0,1),
               c(1,1,1,0,0),
               c(1,0,0,1,1),
               c(0,1,0,1,1),
               c(1,1,0,1,0),
               c(1,1,0,0,1),
               c(1,0,1,1,1),
               c(0,1,1,1,1),
               c(1,1,1,1,0),
               c(1,1,1,0,1),
               c(1,1,0,1,1),
               c(1,1,1,1,1))
      result_elasticnet<-cbind(rci,p,l,u,corr)
      rownames(result_elasticnet)<-rnames
      
      
      
      sortedRCI<-sort(rci,index.return=TRUE,decreasing = TRUE)
      guideh<-t(guide)
      rownames(guideh)<-c("RPPA","VIPER","RNASeq","CNV","Mutation")
      corrplot(guideh[,sortedRCI$ix], method="circle", cl.pos = "n")
      #plot(sortedRCI$x,xlab = rnames[sortedRCI$ix],ylab = "rCI")
      df<-data.frame("y"=sortedRCI$x,"x"=rnames[sortedRCI$ix],"upper"=u,"lower"=l)
      ggplot(df, aes(y=df$y, x=factor(df$x, levels=unique(as.character(df$x)) ))) +ylab("rCI")+xlab("")+ggtitle(paste0(" Elastic net (10 CV) - ",paramDrug, "(", validationSet,")"))+
        geom_point(shape=1) +theme_bw() +theme(axis.text.x = element_text(angle = 90,size=12))+   # Use hollow circles,face = "bold"
        geom_errorbar(aes(ymin = df$lower[sortedRCI$ix], ymax = df$upper[sortedRCI$ix]),width = 0.2,position = "dodge")
    }
  }
  #************Fig CRTPV2 drug classes#############
  Fig5<-1 
  if(Fig5==1)
  {
    Fig5path<-paste(dirMain,"./op/proteomics/","Fig5_",Sys.Date(),sep="")
    if (!file.exists(Fig5path)){dir.create(Fig5path)}
    dfCRTPV2Classes<- read.csv("./Data/CRTPV2DrugClassesGDSC1000.csv")#z<-intersect(targetlst,colnames(exprMatSymbols))#zz<-setdiff(targetlst,colnames(exprMatSymbols))
    names(dfCRTPV2Classes)<-c('Drugs','Targets','Class','PossibleGDSC1000Class')
    #uniqueCRTPV2Classes<-unique(dfCRTPV2Classes$PossibleGDSC1000Class)
    uniqueCRTPV2Classes<-names(sort(table(dfCRTPV2Classes$PossibleGDSC1000Class),decreasing = TRUE))
    #View(sort(table(dfCRTPV2Classes$PossibleGDSC1000Class),decreasing = TRUE))
    
    targetlst<-c()
    for (i in 1:length(dfCRTPV2Classes$Targets))
    {
      targetlst<-c(targetlst,unlist(strsplit(as.character(dfCRTPV2Classes$Targets[i]), split=";")))
    }
    targetlst<-unique(targetlst)
    biomarkers.mrTF<-intersect(targetlst,mrsTargetslst)              #"NQO1"  "BRAF"  "ALK"   "EGFR"  "HGF"   "MET"   "MDM2"  "ERBB2" "FLT3"
    biomarkers.rppa<-intersect(targetlst,ccleRPPATarget$Target_Genes)#         "BRAF"         "EGFR"          "MET"   "MDM2"  "ERBB2"
    biomarkers.mclp<-intersect(targetlst,mclpTarget$Target_Genes)              #"BRAF"        "EGFR"          "MET"   "MDM2"  "ERBB2" "FLT3" 
    flag_save<-1
    exprsFlag<-1
    RPPAFlag<-1
    MCLPFlag<-0
    flag_VIPER<-1
    CNVFlag<-1
    MutationFlag<-1
    cutoffdrugclass<-9#10
    #error tanespimycin c<-3 d<-21  /// c<-11 d<-4 //c<-22 //c<23
    #for (c in 1: length(uniqueCRTPV2Classes)) 
    for (c in 10: 15) 
    {
      classname<-stri_replace_all_fixed(uniqueCRTPV2Classes[c], "(", " ")
      classname<-stri_replace_all_fixed(classname, ")", " ")
      classname<-stri_replace_all_fixed(classname, "(", " ")
      classname<-stri_replace_all_fixed(classname, ":", " ")
      classname<-stri_replace_all_fixed(classname, "/", " ")
      #classname<-stri_replace_all_fixed(classname, "@\", " ")
      
      #if(!classname=="other" & !classname %in% c("apoptosis regulation","DNA replication","metabolism","RTK signaling","chromain  histone acetylation"))
      #{
      classItems<-dfCRTPV2Classes[which(dfCRTPV2Classes$PossibleGDSC1000Class %in% uniqueCRTPV2Classes[c]),]
      #if(length(classItems$Drugs)>cutoffdrugclass)#15
      if(classname %in% c("cell cycle","EGFR signaling","ERK MAPK signaling","Genome integrity","mitosis","p53 pathway"))
      {
        drugclassPath<-paste(Fig5path,"/",classname,sep="")
        if (!file.exists(drugclassPath)){dir.create(drugclassPath)}
        for (d in 1: length(classItems$Drugs)) 
        {
          print(paste0(c," ",classname," ",d))
          ##### Hassan#### drugindx<-11  #comment
          paramDrug<-as.character(classItems$Drugs[d])
          validationSet<-unlist(strsplit(as.character(classItems$Targets[d]), split=";"))      # get drug targets 
          
          if(flag_VIPER==1)
          {
            ##hassan only for mrs #########    validationSet<-validationSet[2]    #####important############
            
            interestSet<-c()
            for(i in 1:dim(mrsTargets)[1])
            {
              validIntersect<-intersect(unlist(strsplit(as.character(mrsTargets[i,1]), split=" ")),validationSet)
              if (length(validIntersect>0))
              {
                interestSet<-c(interestSet,rownames(mrsTargets)[i])
                #print(rownames(mrsTargets)[i])
              }
            }
          }
          drugSavePath<-paste(drugclassPath,"/", as.character(d),sep="")
          if (!file.exists(drugSavePath)){dir.create(drugSavePath)}
          
          #yy<-subset(Biomark4valid, compound=paramDrug)
          if(exprsFlag==1)
          {
            ########### 1- A}, error=function(e){})AC vs exprs ########### 
            tryCatch({
              saveFile<-paste(drugSavePath,"/1-AAC_VS_exprs",sep="")
              if(!file.exists(file=paste0(saveFile,".csv")))
              {
                df<-correlation.mci.AAC(AAC_CTRPv2,eset,paramDrug)#df<-data.frame("e"=df$e,"p"=df$p,"id"=df$id,"lbl"=df$lbl,"p_fdr"=df$p_fdr)
                if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))
                  print(saveFile)}
                ##plotdf(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
                plotVolc.ggrepel(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
                #plot.df(df,0.5,validationSet,"mCI Gene-expression and AAC",n)
              }
            }, error=function(e){print("error")})
          }
          if(CNVFlag==1)
          {
            ########### 1- A}, error=function(e){})AC vs exprs ########### 
            tryCatch({
              saveFile<-paste(drugSavePath,"/1-AAC_VS_CNV",sep="")
              if(!file.exists(file=paste0(saveFile,".csv")))
              {
                df<-correlation.mci.AAC(AAC_CTRPv2,eset_CNV,paramDrug)#df<-data.frame("e"=df$e,"p"=df$p,"id"=df$id,"lbl"=df$lbl,"p_fdr"=df$p_fdr)
                if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))
                  print(saveFile)}
                ##plotdf(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
                plotVolc.ggrepel(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
                #plot.df(df,0.5,validationSet,"mCI Gene-expression and AAC",n)
              }
            }, error=function(e){print("error")})
          }
          if(MutationFlag==1)
          {
            ########### 1- A}, error=function(e){})AC vs exprs ########### 
            tryCatch({
              saveFile<-paste(drugSavePath,"/1-AAC_VS_Mut",sep="")
              if(!file.exists(file=paste0(saveFile,".csv")))
              {
                df<-correlation.mci.AAC(AAC_CTRPv2,eset_Mutation,paramDrug)#df<-data.frame("e"=df$e,"p"=df$p,"id"=df$id,"lbl"=df$lbl,"p_fdr"=df$p_fdr)
                if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))
                  print(saveFile)}
                ##plotdf(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
                plotVolc.ggrepel(df,0.5,validationSet,"mCI Gene-expression and AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
                #plot.df(df,0.5,validationSet,"mCI Gene-expression and AAC",n)
              }
            }, error=function(e){print("error")})
          }
          if(RPPAFlag==1)
          { 
            ########### 2- AAC vs ccleRPPA ########### 
            saveFile<-paste(drugSavePath,"/2-AAC_VS_ccleRPPA",sep="")
            df<-correlation.mci.AAC(AAC_CTRPv2,ccleRPPA,paramDrug)
            if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
            ##plotdf(df,0.5,validationSet,"mCI CCLE-RPPA,AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
            plotVolc.ggrepel(df,0.5,validationSet,"mCI CCLE-RPPA,AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
            
            df<-correlation.mci.AAC.targets(AAC_CTRPv2,ccleRPPA,paramDrug,alias.prepare(ccleRPPATarget))
            if(flag_save==1){write.csv(df,paste0(saveFile,"_targets",".csv"))}
            ##plotdf(df,0.5,validationSet,"mCI_Targets CCLE-RPPA,AAC",n,paste0(saveFile,"_target_plot.tiff",sep=""))
            plotVolc.ggrepel(df,0.5,validationSet,"mCI_Targets CCLE-RPPA,AAC",n,paste0(saveFile,"_target_volc.jpg",sep=""))
            
            #dfforest<-subset(df, id %in% validationSet)
            #plotforest(subset(df, id %in% validationSet),"Title")#/+++++++++++++++++++++++++++++11/3/19
            #commonCCLERPPAgenes<-intersect(fData(rnaseqData)[,"Symbol"],ccleRPPATarget$Antibody_Name)
          }
          if(MCLPFlag==1)
          { 
            ########### 3- AAC vs MCLP ########### 
            saveFile<-paste(drugSavePath,"/3-AAC_VS_MCLP",sep="")
            df<-correlation.mci.AAC(AAC_CTRPv2,mclpData,paramDrug)
            if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
            ##plotdf(df,0.5,validationSet,"mCI MCLP-RPPA,AAC",n,paste0(saveFile,"_plot.tiff",sep="")) 
            plotVolc.ggrepel(df,0.5,validationSet,"mCI MCLP-RPPA,AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
            
            df<-correlation.mci.AAC.targets(AAC_CTRPv2,mclpData,paramDrug,alias.prepare(mclpTarget))
            if(flag_save==1){write.csv(df,paste0(saveFile,"_targets",".csv"))}
            plotVolc.ggrepel(df,0.5,validationSet,"mCI_Targets MCLP-RPPA,AAC",n,paste0(saveFile,"_volc_plot.jpg",sep=""))
            ##plotdf(df,0.5,validationSet,"mCI_Targets MCLP-RPPA,AAC",n,paste0(saveFile,"_target_plot.tiff",sep=""))  
          }
          if(flag_VIPER==1)
          { 
            ########### 4- AAC vs mrs ########### 
            saveFile<-paste(drugSavePath,"/4-AAC_VS_mrs",sep="")
            df<-correlation.mci.AAC(AAC_CTRPv2,mrs,paramDrug)
            if(flag_save==1){write.csv(df,paste0(saveFile,".csv"))}
            ##plotdf(df,0.5,validationSet,"mCI mrs_TF,AAC",n,paste0(saveFile,"_plot.tiff",sep=""))
            plotVolc.ggrepel(df,0.5,validationSet,"mCI mrs_TF,AAC",n,paste0(saveFile,"_volc.jpg",sep=""))
            
            for (i in 1: length(df$lbl)) 
            {
              intersectSet<-intersect(interestSet,df$lbl[i])
              if(length(intersectSet)==0)
              {
                df$lbl[i]<-""
              } else{  
              }
            }
            if(flag_save==1){write.csv(df,paste0(saveFile,"_interesting_targets",".csv"))}
            ##Hassan##df<-correlation.mci.AAC.targets(AAC,mrs,paramDrug,mrsTargets)
            ##Hassan##if(flag_save==1){write.csv(df,paste0(saveFile,"_targets",".csv"))}
            ##Hassan##plotdf(df,0.5,validationSet,"mCI_Targets mrs,AAC",n,paste0(saveFile,"_target_plot.tiff",sep=""))
            #}
          }
        }
      }
    }
  }
}

elasticnet_RCI<-function(x,y)
{
  naINDX<-which(is.na(x), arr.ind=TRUE)[,1]
  if(length(naINDX)>0)
  {x<-x[-naINDX,]
  y<-y[-naINDX]}
  
  #x<-x[!is.na(x[,1]),drop=F]
  #y<-y[!is.na(x[,1]),drop=F]
  set.seed(0123)
  library(caret)
  library(glmnet)
  data_ctrl <- trainControl(method = "cv", number = 5)
  idx<-createFolds(x[,1],k=10)
  y_all<-numeric()
  y_predicted_all<-numeric()
  for (k in 1:length(idx)) {
    fold=idx[[k]]
    testdata=as.matrix(x[fold,])
    traindata=as.matrix(x[-fold,])
    colnames(testdata)<-colnames(x)
    colnames(traindata)<-colnames(x)
    cvfit = cv.glmnet(x=traindata, y=y[-fold],nfolds=5)
    #plot(cvfit)
    cvfit$lambda.min
    coef(cvfit, s = cvfit$lambda.min)
    y_predicted=predict(cvfit, newx = testdata, s = "lambda.min")
    
    #y_predicted=predict(model_caret$finalModel, newdata = data.frame(testdata), type = "response")
    y_all<-append(y_all,y[fold])
    y_predicted_all<-append(y_predicted_all,y_predicted)
  }
  corr_person<-cor.test(y_all,y_predicted_all)
    ci<-mCI::paired.concordance.index(predictions = y_predicted_all,observations =y_all,delta.obs = 0.12,delta.pred = 0,outx = FALSE)#
  df<-data.frame("rCI"=ci$cindex,"p"=ci$p.value,"u"=ci$upper,"l"=ci$lower,"corr"=corr_person$estimate)
  return(df)
}

GET.OUTPUT.DF<-function()
{
  # read csv files for drug
  dirMain<-'D:/+BHK+/R scripts/'
  setwd(dirMain)
  Biomark4valid<- read.csv("./Data/Ben-biomarkers-validation.csv", row.names=1)
  uniqueCompounds<-unique(Biomark4valid$compound)
  df<-data.frame()
  for (drugindx in 1: length(uniqueCompounds)) 
  {    
    paramDrug<-as.character(uniqueCompounds[drugindx])
    drugSavePath<-paste(dirMain,"./op/proteomics/","2019-06-12","/", paramDrug,sep="")#Sys.Date()
    validationSet<-unique(Biomark4valid[which(Biomark4valid[,"compound"]==paramDrug),"gene"])
    for(m in 1:length(validationSet))
    {
      if(!validationSet[m]=="BCR_ABL")
      {
        print(paramDrug)
        dfcnv<-read.csv(paste(drugSavePath,"/1b-AAC_VS_exprs","_CNV.csv",sep=""))
        objcnv<-dfcnv[which(dfcnv$id %in% validationSet[m]),]
        if(dim(objcnv)[1]==0){objcnv<-data.frame("e"=NA,"p"=NA,"fdr"=NA,l=NA,u=NA)}
        if(!is.na(objcnv$e) & objcnv$e<0.5){
          print(paste("F CNV ",objcnv$l," ",objcnv$e," ",objcnv$u))
          u=(1-objcnv$l)
          l=(1-objcnv$u)
          objcnv$l=l
          objcnv$u=u
          objcnv$e=1-objcnv$e
          print(paste("T CNV ",objcnv$l," ",objcnv$e," ",objcnv$u))
        }
        
        dfrna<-read.csv(paste(drugSavePath,"/1a-AAC_VS_exprs","_RNASeq.csv",sep=""))
        objrna<-dfrna[which(dfrna$id %in% validationSet[m]),]
        if(dim(objrna)[1]==0){objrna<-data.frame("e"=NA,"p"=NA,"fdr"=NA,l=NA,u=NA)}
        if(!is.na(objrna$e) & objrna$e<0.5){
          print(paste("F rna ",objrna$l," ",objrna$e," ",objrna$u))
          u=(1-objrna$l)
          l=(1-objrna$u)
          objrna$l=l
          objrna$u=u
          objrna$e=1-objrna$e
          print(paste("T rna ",objrna$l," ",objrna$e," ",objrna$u))
        }
        
        dfmut<-read.csv(paste(drugSavePath,"/1c-AAC_VS_exprs","_Mutation.csv",sep=""))
        objmut<-dfmut[which(dfmut$id %in% validationSet[m]),]
        if(dim(objmut)[1]==0){objmut<-data.frame("e"=NA,"p"=NA,"fdr"=NA,l=NA,u=NA)}
        if(!is.na(objmut$e) & objmut$e<0.5){
          print(paste("F Mut ",objmut$l," ",objmut$e," ",objmut$u))
          
          if(objmut$l==0 & objmut$u==0){l=u=1-objmut$e}
          else{      
            u=(1-objmut$l)
            l=(1-objmut$u)
          }
          objmut$l=l
          objmut$u=u
          objmut$e=1-objmut$e
          print(paste("T Mut ",objmut$l," ",objmut$e," ",objmut$u))}
        
        
        dfcclerppa<-read.csv(paste(drugSavePath,"/2-AAC_VS_ccleRPPA_targets.csv",sep=""))
        objclerppa<-dfcclerppa[which(dfcclerppa$id %in% validationSet[m]),]
        objclerppa<-objclerppa[which.max(objclerppa$e),]
        if(dim(objclerppa)[1]==0){objclerppa<-data.frame("e"=NA,"p"=NA,"fdr"=NA,l=NA,u=NA)}
        if(!is.na(objclerppa$e) & objclerppa$e<0.5){
          print(paste("F rppa ",objclerppa$l," ",objclerppa$e," ",objclerppa$u))
          u=(1-objclerppa$l)
          l=(1-objclerppa$u)
          objclerppa$l=l
          objclerppa$u=u
          objclerppa$e=1-objclerppa$e
          print(paste("T rppa ",objclerppa$l," ",objclerppa$e," ",objclerppa$u))
        }
        
        if(length(validationSet)>1 & !paramDrug %in% c("TAE684","AZD0530") )
        {dfviper<-read.csv(paste0(drugSavePath,"/4-AAC_VS_mrs","_interesting_targets_",validationSet[m],".csv"))
        }else{dfviper<-read.csv(paste0(drugSavePath,"/4-AAC_VS_mrs","_interesting_targets.csv"))}
        
        objviper<-dfviper[which(!is.na(dfviper$lbl) & dfviper$fdr<0.05),]
        if(dim(objviper)[1]==0){objviper<-dfviper[which(!is.na(dfviper$lbl)),]}
        objviper<-objviper[which.max(objviper$e),]
        if(dim(objviper)[1]==0){objviper<-data.frame("e"=NA,"p"=NA,"fdr"=NA,l=NA,u=NA)}
        if(!is.na(objviper$e) & objviper$e<0.5){
          print(paste("F viper ",objclerppa$l," ",objclerppa$e," ",objclerppa$u))
          u=(1-objviper$l)
          l=(1-objviper$u)
          objviper$l=l
          objviper$u=u
          objviper$e=1-objviper$e
          print(paste("T viper ",objclerppa$l," ",objclerppa$e," ",objclerppa$u))
        }
        
        df<-rbind.data.frame(df,data.frame("Drug"=paramDrug,	"Biomarker"=validationSet[m],	
                                           "mci.RPPA"=objclerppa$e,"p.RPPA"=objclerppa$p,"FDR.RPPA"=objclerppa$fdr,"lower.RPPA"=objclerppa$l,"upper.RPPA"=objclerppa$u,	
                                           "mci.VIPER"=objviper$e,"p.VIPER"=objviper$p,"FDR.VIPER"=objviper$fdr,"lower.VIPER"=objviper$l,"upper.VIPER"=objviper$u,
                                           "mci.RNASEQ"=objrna$e,"p.RNASEQ"=objrna$p,	"FDR.RNASEQ"=objrna$fdr,"lower.RNASEQ"=objrna$l,"upper.RNASEQ"=objrna$u,
                                           "mci.CNV"=objcnv$e,"p.CNV"=objcnv$p,"FDR.CNV"=objcnv$fdr,"lower.CNV"=objcnv$l,"upper.CNV"=objcnv$u,
                                           "mci.Mutation"=objmut$e,"p.Mutation"=objmut$p,"FDR.Mutation"=objmut$fdr,"lower.Mutation"=objmut$l,"upper.Mutation"=objmut$u))
        #
      }
    }
  }
  
  return(df)
}
plotforestAnt<-function(dfFOREST,title)
{
  label <- c("Mutation","CNV","RNASEQ","VIPER","RPPA")#"MCLP",#drugsbiomarkers #paste0("X", 1:6)#Drug-validationSet
  
  #data=dfFOREST, aes(x=dfFOREST$label, y=dfFOREST$mean, ymin=dfFOREST$lower, ymax=dfFOREST$upper
  c_indices <- structure(
    list(
      mean  = c(NA,dfFOREST$mean),
      lower = c(NA,dfFOREST$lower),
      upper = c(NA,dfFOREST$upper)#c(NA,dfFOREST$upper)
    ),
    .Names = c("C-index    ", "lower", "upper"),
    row.names = c("omics",label), 
    class = "data.frame"
  )
  
  c_tabletext <- cbind(
    c("omics",as.character(dfFOREST$label)) ,
    c("mCI",as.character(formatC(dfFOREST$mean, digits = 2))) ,
    c("FDR",as.character(formatC(dfFOREST$fdr, digits = 2))) #common samples for each dataset
    # c("C-index", formatC(GRAYci, format = "e", digits = 2), formatC(UHNci, format = "e", digits = 2), formatC(combined_ci$estimate, format = "e", digits = 2)),
    #  c("P-value", formatC(GRAYpvalue, format = "e", digits = 2), formatC(UHNpvalue, format = "e", digits = 2), formatC(combined_ci_p, format = "e", digits = 2))
  )
  
  drugSavePath<-paste(dirMain,"./op/proteomics/",Sys.Date(),"_Fig3",sep="")
  fileName<- paste0(drugSavePath,"/",title,".pdf")
  pdf(fileName, width=7, height=3, onefile=FALSE)
  
  forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(F,F,F,F,T,F), xlab="Concordance Index", 
             title=title, zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:4, col = "#000044")),
             txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                            ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                            xlab=gpar(fontfamily = "", cex=0.8, fontface=2),
                            legend=gpar(fontfamily = "", cex = 1, fontface=1)),
             col=fpColors(box=RColorBrewer::brewer.pal(n=4, name="Set2"),
                          line=RColorBrewer::brewer.pal(n=4, name="Set2"),
                          summary="blue"),
             xticks= c( .2, .4, .5, .6, .7, .8)
  )
  
  dev.off()
}
plotforestop<-function(df,title)#drugsbiomarkers,mCI,lower,upper
{
  library(ggplot2)
  fp <- ggplot(data=df, aes(x=df$label, y=df$mean, ymin=df$lower, ymax=df$upper)) +
    geom_pointrange(aes(col=df$label)) + 
    geom_hline(yintercept=0.5, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab(title) + ylab("mCI (95% Confidence)") +#theme_bw()#+ #+  # use a white background
    theme(plot.title=element_text(size=16,face="bold"),
          #axis.text.y=element_blank(),
          #axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme_bw()
  print(fp)
}
plotforest<-function(dataf,title)#drugsbiomarkers,mCI,lower,upper
{
  label <- c("MCLP-RPPA","CCLE-RPPA","VIPER","RNASeq","CNV","Mutation")
  mean  <- c(0.445,0.468,0.451,0.446,0.522,0.789)
  lower <- c(0.391,0.431,0.416,0.410,0,0)
  upper <- c(0.499,0.506,0.486,0.482,0,0)
  
  library(forestplot)
  label <- dataf$id#drugsbiomarkers #paste0("X", 1:6)#Drug-validationSet
  mean  <- dataf$e 
  lower <- dataf$l
  upper <- dataf$u 
  
  df <- data.frame(label, mean, lower, upper)
  
  # reverses the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label))
  
  library(ggplot2)
  fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    geom_pointrange() + 
    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("PD-0325901 (BRAF)") + ylab("mCI (95% Confidence)") +
    theme_bw()  # use a white background
  print(fp)
  
}
read_op_files<-function()
{
  fls <- list.files(path="./op/proteomics/2019-06-13/",pattern = ".csv", recursive = TRUE)
  for (i in 1:length(fls)){
    filepath<-paste0("./op/proteomics/2019-06-13/",fls[i],sep="")
    df<- read.csv(filepath)
    fdr<-p.adjust(df$p,method="fdr")
    df<-data.frame("e"=df$e,"p"=df$p,"fdr"=fdr,"l"=df$l,"u"=df$u,"id"=df$id,"lbl"=df$lbl)
    write.csv(df,filepath)
  }
}

plot_op_labmeeting<-function()
{
  library(RColorBrewer)
  coul = brewer.pal(5, "Pastel2") 
  lbls<-paste0(df$Biomarker,"_",df$Drug,sep = "")
  
  par(mar=c(9, 4.1,6.1, 1.1))
  
  plotdata<-rbind(df$mci.RPPA,df$mci.VIPER,df$mci.RNASEQ,df$mci.CNV,df$mci.Mutation)
  fdrdata<-rbind(df$FDR.RPPA,df$FDR.VIPER,df$FDR.RNASEQ,df$FDR.CNV,df$FDR.Mutation)
  mids<-barplot(plotdata,col=coul ,beside = TRUE,ylab = "mCI",names.arg = paste(df$Biomarker,"-",df$Drug,sep = ""),las=2, xpd = FALSE,ylim=c(0.5, 1.0))
  text(x=mids , y= plotdata+0.01, labels=ifelse(fdrdata<0.01,'**',ifelse(fdrdata<0.05,'*','')), xpd=FALSE)
  legend("topright",
         c("CCLE-RPPA","CCLE-VIPER","rnaseq","CNV","Mutation"),
         fill = coul,cex = 0.5)
  #---------------------------------------------------------------
  lenDrugs<-length(paste(df$Biomarker,"-",df$Drug,sep = ""))
  datatype<-c(rep("CCLE-RPPA",lenDrugs),rep("CCLE-VIPER",lenDrugs),rep("rnaseq",lenDrugs),rep("CNV",lenDrugs),rep("Mutation",lenDrugs))
  drugs<-rep(paste(df$Biomarker,"-",df$Drug,sep = ""),5)
  wCI<-c(plotdata[1,],plotdata[2,],plotdata[3,],plotdata[4,],plotdata[5,])
  data<-data.frame(drugs,datatype,wCI)
  ggplot(data,aes(fill=data$datatype,y=data$wCI,x=data$drugs))+geom_bar(position="stack",stat="identity",alpha=0.5)+
    theme(axis.text.x = element_text(angle = 90)) 
  #---------------------------------------------------------------
  to_plot<-t(as.data.frame(plotdata[,-9]))
  rownames(to_plot)<-ylbls[-9]
  colnames(to_plot)<-c("RPPA","VIPER","RNASeq","CNV","Mutation")
  melted<-melt(to_plot, id="x")
  colnames(melted)<-c("Drugs","Type","mCI")
  bp<-print(ggplot(melted,aes(x=Drugs,y=mCI,fill=Type)) + 
          geom_bar(stat="identity",position = "identity", alpha=0.9))+theme_bw()+
    theme(axis.text.x = element_text(angle = 90))
  bp+ coord_cartesian( ylim = c(0.5, 1))
  #text(x=bp , y= to_plot+0.01, labels=ifelse(fdrdata<0.01,'**',ifelse(fdrdata<0.05,'*','')), xpd=FALSE)

  
  #---------------------------------------------------------------

  plotdata<-plotdata[,-9]
  barplot(plotdata[2,], col="blue",alpha=0.5,width=0.2,space=.1)#,
  barplot(plotdata[5,], col="cyan", add=TRUE,width=0.1,space=.2)#width=0.3,
  barplot(plotdata[4,], col="yellow", add=TRUE)#width=0.3,
  barplot(plotdata[3,], col="green", add=TRUE)
   barplot(plotdata[1,], col="red", add=TRUE,)#, width=0.1
#  barplot(plotdata[2,], width=0.3, space=c(0.9, 1.4, 1.4, 1.4, 1.4), col="blue", add=TRUE)
  
   barplot(plotdata[1,],border="red",density=0)
   par(new=TRUE)
   barplot(plotdata[2,],border="green",density=0)
   par(new=TRUE)
   barplot(plotdata[3,],border="blue",density=0)
   par(new=TRUE)
   barplot(plotdata[4,],border="yellow",density=0)
   par(new=TRUE)
   barplot(plotdata[5,],border="cyan",density=0)
   
   ###########################################used to Ben (transparent)
  ylbls<- paste(df$Biomarker,"-",df$Drug,sep = "")
 # ylbls<-ylbls[-9]
   barplot(plotdata[2,],col="white",border="green",las=2, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[3,],col="white",border="blue",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[4,],col="white",border="black",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[5,],col="white",border="orange",add=T,las=2,lty =4, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[1,],col="white",border="red",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0),ylab = "wCI",names.arg = ylbls)
   legend("topright",
          c("CCLE-RPPA","CCLE-VIPER","rnaseq","CNV","Mutation"),
          fill = c("red","green","blue","black","orange"),cex = 1,text.width = 12,bty = "n",border=NA)
   ###########################################
   dat <- read.table(text = "A   B   C
                      1 4 1.5
                     2 3 1 
                     3 2 1
                     1 2 3", header = TRUE)
   
   barplot(dat$A, col=rgb(1, 0, 0, 0.3)) 
   barplot(dat$B, col=rgb(0, 1, 0, 0.3), add=TRUE) 
   barplot(dat$C, col=rgb(0, 0, 1, 0.3), add=TRUE)
   
   # adding a legend
   legend('topright', bty = 'n', title = 'Legend',
          legend = c('A', 'B', 'C'), fill = c('red', 'green','blue'))
   ###########################################
   ylbls<- paste(df$Biomarker,"-",df$Drug,sep = "")
   # ylbls<-ylbls[-9]
   barplot(plotdata[2,],col="white",border="green",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[3,],col="white",border="blue",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[4,],col="white",border="black",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[5,],col="white",border="orange",add=T,las=2,lty =4, xpd = FALSE,ylim=c(0.5, 1.0))
   barplot(plotdata[1,],col="white",border="red",add=T,las=2, xpd = FALSE,ylim=c(0.5, 1.0),ylab = "wCI",names.arg = ylbls)
   legend("topright",
          c("CCLE-RPPA","CCLE-VIPER","rnaseq","CNV","Mutation"),
          fill = c("red","green","blue","black","orange"),cex = 1,text.width = 12,bty = "n",border=NA)
   barplot(plotdata[2,], col=rgb(0, 1, 0, 0.8),ylim=c(0.5, 1.0),las=2, xpd = FALSE) 
   barplot(plotdata[3,], col=rgb(1, 0, 0, 0.8), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE) 
   barplot(plotdata[4,], col=rgb(0, 0, 1, 0.8), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE)
   barplot(plotdata[5,], col=rgb(0, 1, 1, 0.8), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE)
   barplot(plotdata[1,], col=rgb(1, 1, 0, 0.8), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE,ylab = "wCI",names.arg = ylbls)
   # adding a legend
   legend('topleft', bty = 'n', title = 'Data type',
          legend = c('RPPA', 'VIPER', 'RNASeq','CNV','Mutation'), fill = c('yellow', 'green','blue','red','cyan'))
   legend('topright', bty = 'n', title = 'Legend',
          legend = c('A', 'B', 'C'), fill = c('red', 'green','blue'))
   transparency<-1
   barplot(plotdata[2,], col=rgb(0, 1, 0, transparency),ylim=c(0.5, 1.0),las=2, xpd = FALSE) 
   barplot(plotdata[3,], col=rgb(0, 0, 1, transparency), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE)
   barplot(plotdata[4,], col=rgb(1, 0, 0, transparency), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE)
   barplot(plotdata[1,], col=rgb(1, 1, 0, transparency), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE) 
   barplot(plotdata[5,], col=rgb(0, 1, 1, transparency), add=TRUE,ylim=c(0.5, 1.0),las=2, xpd = FALSE,ylab = "wCI",names.arg = ylbls)
   legend('topleft', bty = 'n', title = 'Data type',
          legend = c('RPPA', 'VIPER', 'RNASeq','CNV','Mutation'), fill = c('yellow', 'green','blue','red','cyan'))
   
   
   library(reshape2)
   data$id <- 1:nrow(data)
   DF1 <- melt(data, id.var="id")
   p <- ggplot(DF1, aes(x = Rank, y = value, fill = variable)) +
     geom_bar(stat = "identity")
   d
   p <- ggplotly(p)
   
   
   #####################ggplot
   library(ggplot2)
   library(reshape)
   x = c("Band 1", "Band 2", "Band 3")
   y1 = c("1","2","3")
   y2 = c("2","3","4")
   to_plot <- data.frame(x=x,y1=y1,y2=y2)
   melted<-melt(to_plot, id="x")
   ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat="identity", alpha=.3)
   print(ggplot(melted,aes(x=x,y=value,fill=variable)) + 
           geom_bar(stat="identity",position = "identity", alpha=.3))
}
plot_op<-function()
{
  df_old<- read.csv("D:/+BHK+/R scripts/op/proteomics/proteomics_op_forest-noleukemia.csv")#D:/+BHK+/R scripts/op/proteomics_op_results.csv
  df_old<- df_old[-c(6,10,12),]
  colnames(df_old)[1]<-"Drug"
  
  df<-GET.OUTPUT.DF()
  df<- df[-c(9),]

  library(forestplot)
  for(i in 1:dim(df)[1])
  {
    label <- c("Mutation","CNV","RNASEQ","VIPER","RPPA")#"MCLP",#drugsbiomarkers #paste0("X", 1:6)#Drug-validationSet
    mean  <- c(df[i,"mci.Mutation"],df[i,"mci.CNV"],df[i,"mci.RNASEQ"],df[i,"mci.VIPER"],df[i,"mci.RPPA"])#df[i,"mci.MCLP"], #c(1.29,0.76,2.43,1.68,1.22,1.7)#mCI
    lower <- c(df[i,"lower.Mutation"],df[i,"lower.CNV"],df[i,"lower.RNASEQ"],df[i,"lower.VIPER"],df[i,"lower.RPPA"])#df[i,"lower.MCLP"], #c(0.84,0.50,1.58,1.1,0.8,1.11)#p
    upper <- c(df[i,"upper.Mutation"],df[i,"upper.CNV"],df[i,"upper.RNASEQ"],df[i,"upper.VIPER"],df[i,"upper.RPPA"])#df[i,"upper.MCLP"], #c(1.95,1.16,3.67,2.54,1.85,2.56)#p
    fdr<-c(df[i,"FDR.Mutation"],df[i,"FDR.CNV"],df[i,"FDR.RNASEQ"],df[i,"FDR.VIPER"],df[i,"FDR.RPPA"])
    title<-paste0(df[i,"Drug"]," (",df[i,"Biomarker"],")",sep="")
    dfFOREST<- data.frame(label, mean, lower, upper,fdr)
    plotforestAnt(dfFOREST,title)
  }
  #create color palette:
  library(RColorBrewer)
  coul = brewer.pal(5, "Pastel2") # grey.colors(6)
  lbls<-paste0(df$Biomarker,"_",df$Drug,sep = "")
  par(mar=c(9, 4.1,6.1, 1.1))
  plotdata<-rbind(df$mci.RPPA,df$mci.VIPER,df$mci.RNASEQ,df$mci.CNV,df$mci.Mutation)#df$mci.MCLP,
  fdrdata<-rbind(df$FDR.RPPA,df$FDR.VIPER,df$FDR.RNASEQ,df$FDR.CNV,df$FDR.Mutation)#df$FDR.MCLP,
  plotdata_tmp<-plotdata-0.5
  mids<-barplot(plotdata_tmp,col=coul ,beside = TRUE,ylab = "mCI",names.arg = paste(df$Biomarker,"-",df$Drug,sep = ""),las=2,axes=FALSE, xpd = FALSE,ylim = c(min(plotdata_tmp,na.rm = T)-0.05,max(plotdata_tmp,na.rm = T)+0.05))#,offset = 0.5
  axis(side=2,at=seq(0,1,0.1)-0.5,labels=(seq(0,1,0.1)))
  text(x=mids , y= ifelse(plotdata_tmp>0,plotdata_tmp+0.01,plotdata_tmp-0.01), labels=ifelse(fdrdata<0.01,'**',ifelse(fdrdata<0.05,'*','')), xpd=FALSE)
  legend("topright",legend =  c("CCLE-RPPA","CCLE-VIPER","rnaseq","CNV","Mutation")
         ,bty = "n", fill = c(coul),cex = 0.3)
  legend("topleft",legend =  c("
                               ** FDR < 0.01"," * FDR < 0.05")
         ,bty = "n",cex = 0.5)
 
  
  coul = grey.colors(6)
  par(mar=c(9, 4.1, 2.1, 1.1))
  plotdata<-rbind(df$p.MCLP,df$p.RPPA,df$p.VIPER,df$p.RNASEQ,df$p.CNV,df$p.Mutation)
  val<-replacezeroPval(plotdata) #0.00005
  plotdata[which(plotdata==0)]<-val-0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
  
  plotdata<--log10(plotdata)
  mids<-barplot(plotdata,col=coul ,beside = TRUE,ylab = "-log10(P-Value)",names.arg = paste(df$Biomarker,"-",df$Drug,sep = ""),las=2, xpd = TRUE)
  legend("topright",
         c("MCLP-RPPA","CCLE-RPPA","CCLE-VIPER","rnaseq","CNV","Mutation"),
         fill = coul,cex = 0.3)
  #text(x=mids , y= plotdata+0.01, labels="*", xpd=TRUE)
  text(x=mids , y= plotdata+0.3, labels=ifelse(fdrdata<0.005,'*',ifelse(fdrdata<0.05,'+','')), xpd=TRUE)
}

replacezeroPval<-function(plotdata)
{
  sortedVals<-sort(colMins(plotdata, na.rm = TRUE))
  val<- sortedVals[which(sortedVals>0)[1]]
  return(val)
  
}

mrsTarget.prepare<-function(mrs.unnormalized,regulonEX)
{
  mrsTargets<-c()
  #indx<-0
  for (i in mrs.unnormalized)
  {
    mrsTargets<-c(mrsTargets,paste0((rownames(as.matrix(regulonEX[[i]]$tfmode))),collapse=" "))
    #mrsTargets<-c(mrsTargets,list(rownames(as.matrix(regulonEX[[i]]$tfmode))))
    #mrsTargets<-c(mrsTargets,paste(rownames(as.matrix(regulonEX[[i]]$tfmode)),sep=" "))
    #indx<-indx+1
    #names(mrsTargets[indx])<-i
  }#names(mrsTargets)<-mrs.unnormalized #(list)
  mrsTargets<-as.matrix(mrsTargets)
  return(mrsTargets)
}

alias.prepare<-function(dfTarget)
{
  alias<-as.matrix(dfTarget$Target_Genes)
  rownames(alias)<-dfTarget$Antibody_Name
  return(alias)
}
targetlst.get<-function(targets)
{
  # purpose for venn diagrams to show number of commom target_genes of CCLE-RPPA for instance and RNASeq_Genes
  #targets<-ccleRPPATarget$Target_Genes
  targetlst<-c()
  for (i in 1:length(targets))
  {
    targetlst<-c(targetlst,unlist(strsplit(as.character(targets[i]), split=" ")))
  }
  
  #z<-intersect(unique(targetlst),colnames(exprMatSymbols))
  return(unique(targetlst))
}

correlation.mci.AAC.targets<-function(A,B,paramDrug,alias)#used
{
  commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
  #UNcommonGP<-setdiff(colnames(mclpData),fData(rnaseqData)[,"Symbol"])#length=147
  id<-c()
  lbl<-c()
  e<-numeric()
  p<-numeric()
  u<-numeric()
  l<-numeric()
  df<-data.frame()
  for (i in 1:dim(B)[2])
  {
    tryCatch({target <-as.character(alias[colnames(B)[i],1])},error=function(e) {target<-colnames(B)[i]}, finally = {})
    target<-unlist(strsplit(target, split=" "))
    for(j in 1:length(target))#more than one target
    {
      id<-append(id,target[j])
      ci<-mCI::paired.concordance.index(predictions = as.numeric(B[commonSample,colnames(B)[i]]),observations = A[commonSample,paramDrug,drop=F],delta.obs = 0.12,delta.pred = 0,outx = FALSE)#,outx = FALSE
      e<-append(e,ci$cindex)
      p<-append(p,ci$p.value)
      u<-append(u,ci$upper)
      l<-append(l,ci$lower)
      lbl<-append(lbl,paste(target[j]," ", colnames(B)[i]))
    }
  }
  df<-data.frame("id"=id,e,p,lbl,u,l)
  return(df)
}
correlation.mci.AAC<-function(A,B,paramDrug)#used
{
  commonSample <- intersect(rownames(A[,paramDrug,drop=F]),rownames(B))
  df<-data.frame()
  e<-numeric()
  p<-numeric()
  u<-numeric()
  l<-numeric()
  id<-c()
  lbl<-c()
  for (protein in colnames(B))
  {
    tryCatch({
      ci<-mCI::paired.concordance.index(predictions = as.numeric(B[commonSample,protein]),observations = A[commonSample,paramDrug,drop=F],delta.obs = 0.12,delta.pred = 0,outx = FALSE)#
      id<-append(id,protein)
      e<-append(e,ci$cindex)
      p<-append(p,ci$p.value)
      u<-append(u,ci$upper)
      l<-append(l,ci$lower)
      lbl<-id
    }, error=function(e){})
  }
  df<-data.frame(e,p,u,l,id,lbl)
  return(df)
}
correlate2Mat.alias<-function(A,B,alias)#used
{
  commonSample<-intersect(rownames(A),rownames(B))
  #UNcommonGP<-setdiff(colnames(mclpData),fData(rnaseqData)[,"Symbol"])#length=147
  #dummy<-c()
  lbl<-c()
  id<-c()
  protein<-c()
  t<-numeric()
  e<-numeric()
  p<-numeric()
  df<-data.frame()
  for (i in 1:dim(B)[2])
  {
    tryCatch({target <-alias[colnames(B)[i],1]},error=function(e) {target<-colnames(B)[i]}, finally = {})
    target<-unlist(strsplit(target, split=" "))
    for(j in 1:length(target))#more than one target
    {
      if(target[j] %in% colnames(A))# if the protein target within rnaeq_genes
      {
        tryCatch({
          #dummy<-c(target[j])
          cp<-cor.test(A[commonSample,target[j]],as.numeric(B[commonSample,colnames(B)[i]]))
          e<-append(e,cp$estimate)
          t<-append(t,cp$statistic)
          p<-append(p,cp$p.value)
          protein<-append(protein,colnames(B)[i])
          id<-append(id,paste(target[j]))
          lbl<-append(lbl,paste(target[j]," ", colnames(B)[i]))},error=function(e){},finally={})
        
      }
    }
  }
  df<-data.frame("id"=id,e,p,t,protein,lbl)
  return(df)
}

correlate2Mat.RPPA_VIPER<-function(mrs,ccleRPPA,aracneTarget,ccleRPPATarget,commonTargets)#used
{
  commonSample<-intersect(rownames(mrs),rownames(ccleRPPA))
  #UNcommonGP<-setdiff(colnames(mclpData),fData(rnaseqData)[,"Symbol"])#length=147
  #dummy<-c()
  lbl<-c()
  id<-c()
  protein<-c()
  t<-numeric()
  e<-numeric()
  p<-numeric()
  df<-data.frame()
  #targets<-data.frame("mrs"=rownames(mrsTargets),"regulators"=mrsTargets)
  dfRPPA<-data.frame()
  for(i in 1:length(ccleRPPATarget$Target_Genes))
  {
    dfRPPA<-rbind(dfRPPA,data.frame(ccleRPPATarget$Antibody_Name[i],unlist(strsplit(as.character(ccleRPPATarget$Target_Genes[i]), split=" "))))
  }
  names(dfRPPA)<-c("Antibody","Targets")
  sum<-0
  for (i in 1:length(commonTargets))
  {
    
    A<-dfRPPA[which(dfRPPA$Targets %in% commonTargets[i]),"Antibody"]
    B<-aracneTarget[which(aracneTarget$Target %in% commonTargets[i]),"Regulator"]
    A<-normalizeName(A)
    B<-normalizeName(B)
    sum<-sum+(length(A)*length(B))
    for(j in 1:length(A))
    {
      for(k in 1:length(B))
      {
        #if(B[k]!="C11orf95")
        tryCatch(
          {
            cp<-cor.test(ccleRPPA[commonSample,A[j]],mrs[commonSample,as.character(B[k])])
            e<-append(e,cp$estimate)
            t<-append(t,cp$statistic)
            p<-append(p,cp$p.value)
            protein<-append(protein,paste(A[j]," ", B[k]))
            id<-append(id,paste(commonTargets[i]))
            lbl<-append(lbl,paste(commonTargets[i]," ", A[j]," ", B[k]))
          }, error=function(e){print(as.character(B[k]))})
      }
    }
    
  }
  df<-data.frame("id"=id,e,p,t,protein,lbl)
  return(df)
}

correlate2Mat<-function(A,B)
{
  commonGP<-intersect(colnames(A),colnames(B))
  commonSample<-intersect(rownames(A),rownames(B))#common cell lines
  UNcommonGenes<-setdiff(colnames(A),colnames(B))#length=147
  t<-numeric()
  e<-numeric()
  p<-numeric()
  df<-data.frame()
  for (i in commonGP)
  {
    tryCatch({
      cp<-cor.test(A[commonSample,i],as.numeric(B[commonSample,i]))
      e<-append(e,cp$estimate)
      t<-append(t,cp$statistic)
      p<-append(p,cp$p.value)
    }, error=function(e){})
  }
  df<-data.frame("id"=commonGP,e,p,t)
  return(df)
}

normalizeName<-function(item)#used
{
  item<-item%>% stringr::str_replace("_Caution", "")
  item<-item %>% stringr::str_replace("_(?!p)", "")
  item<-stri_replace_all(item, "", fixed=".")
  item<-stri_replace_all_fixed(item, "_", "")
  item<-stri_replace_all_fixed(item, " ", "")
  item<-stri_replace_all_fixed(item, "-", "")
  item<-stri_replace_all_fixed(item, ")", "")
  item<-stri_replace_all_fixed(item, "(", "")
  item <- sapply(item,toupper)
  return(item)
}

plotVolc.ggrepel<-function(df,vLine,validationSet,plot_title,n,saveFile)
{
  # library(corrplot)
  # library(Rarity)
  # library(Biobase)
  # library(calibrate)
  # library(ggrepel)
  
  df$p[df$p == 0] <- 0.000000000000000005
  p_fdr<-p.adjust(df$p,method="fdr")
  df<-cbind(df,"p_fdr"=p_fdr)
  
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  # Read data from the web
  #df = read.table(url, header=TRUE)
  data_interest<-subset(df, p_fdr<.05 & p_fdr %in% head(sort(df$p_fdr),n))
  data_validation<-subset(df, id %in% validationSet)
  
  df = mutate(df, sig=ifelse(df$p_fdr<0.05, "FDR<0.05", "Not Sig"))
  p = ggplot(df, aes(e, -log10(p))) +
    geom_point(aes(col=sig)) +
    scale_color_manual(values=c("red", "black")) + 
    ggtitle(plot_title) +
    geom_text_repel(data=data_interest, aes(label=id))+
    geom_text_repel(data=data_validation, aes(label=id),col="blue")
  ggsave(saveFile)
}
plotdf<-function(df,vLine,validationSet,title,n,saveFile)
{
  df$p[df$p == 0] <- 0.000000000000000005
  p_fdr<-p.adjust(df$p,method="fdr")
  df<-cbind(df,"p_fdr"=p_fdr)
  
  tiff(saveFile, width = 10, height = 5, units = 'in', res = 1024, compression = 'none')
  validIntersect<-intersect(df[,"id"],validationSet)
  with(df , plot(e, -log10(p), pch=20, cex=0.8,  col="gray", main=title,xlab = "Correlation"))#, xlim=c(-2.5,2))
  abline(v=vLine,col="blue")
  with(subset(df, id %in% validIntersect ), points(e, -log10(p), pch=20, col="orange"))#orange validIntersect(all) 
  with(subset(df, p_fdr<.05 ), points(e, -log10(p), pch=20, col="cyan"))
  with(subset(df, p_fdr<.05 & p_fdr %in% head(sort(df$p_fdr),n)),textxy(e, -log10(p), labs=id, cex=.5))
  #with(subset(df, p_fdr<.05 ),textxy(e, -log10(p), labs=lbl, cex=.5))#cyan pfdr<0.05 
  with(subset(df, p_fdr<.05 & id %in% validIntersect), points(e, -log10(p), pch=20, col="red"))#red validIntersect & pfdr<0.05 
  with(subset(df, p_fdr<.05 & id %in% validIntersect),textxy(e, -log10(p), labs=id, cex=.5,col="red"))
  dev.off()
  #with(subset(df, p_fdr %in% head(sort(p_fdr)[which(p_fdr>0  & p_fdr<.05 )],n)), textxy(e, -log10(p), labs=lbl, cex=.5))
  ####### with(subset(df, p_fdr %in% head(sort(p_fdr)[which(p_fdr>0  & p_fdr<.05 )],n)), textxy(e, -log10(p), labs=lbl, cex=.5))
  #with(subset(df, p_fdr>0 & p_fdr<.05 & df$p_fdr %in% head(df[order(df$p_fdr),],n)), textxy(e, -log10(p), labs=lbl, cex=.5))
  #with(subset(df, p_fdr>0 & p_fdr<.05 ), textxy(e, -log10(p), labs=lbl, cex=.5))
  #with(subset(df,   id %in% validIntersect & id!= "KRAS"), )#, offset=3#& id!= "KRAS"
  
  #legend(1,95,legend=c("FDR<0.05 & Validated","FDR<0.05","Non validated","Non significant"),
  #       col=c("red","cyan","orange","gray"),text.font = 4,bg='lightblue')
}
plot.df<-function(df,vLine,validationSet,title,n)
{
  p_fdr<-p.adjust(df$p,method="fdr")
  df<-cbind(df,"p_fdr"=p_fdr)
  validIntersect<-intersect(df[,"id"],validationSet)
  
  with(df , plot(e, -p, pch=20, cex=0.1,  col="gray", main=title,xlab = "Correlation"))#, xlim=c(-2.5,2))
  abline(v=vLine,col="blue")
  with(subset(df, id %in% validIntersect ), points(e, -p, pch=20, col="orange"))
  with(subset(df, p_fdr<.05 ), points(e, -p, pch=20, col="cyan"))
  with(subset(df, p_fdr<.05 & id %in% validIntersect), points(e, -p, pch=20, col="red"))
  #with(subset(df, p_fdr<.05 & p_fdr %in% head(sort(p_fdr),20)), points(e, p, pch=20, col="green"))
  with(subset(df, p_fdr<.05 & p_fdr %in% head(sort(p_fdr) & p_fdr>0,n)), textxy(e, -p, labs=lbl, cex=.5))
  with(subset(df,   id %in% validIntersect & id!= "KRAS" ), textxy(e, -p, labs=lbl, cex=.5,col="red"))#, offset=3#& id!= "KRAS"
  
  #legend(1,95,legend=c("FDR<0.05 & Validated","FDR<0.05","Non validated","Non significant"),
  #       col=c("red","cyan","orange","gray"),text.font = 4,bg='lightblue')
}