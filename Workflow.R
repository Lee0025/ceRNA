rm(list = ls())
options(stringsAsFactors = F)
# Attentions: this pipeline were the streamlined steps
# data preparsion
## data assess and gene ID transform
library(rjson) 
library(dplyr) 
library(limma) 
library(stringr) 
###---1.dispose json files----------------------
jsonFile <- fromJSON(file="./metadata.cart.2020-05-13.json")
filesNameToBarcode <- data.frame(filesName = c(),TCGA_Barcode = c())
for(i in 1:length(jsonFile)){
  TCGA_Barcode <- jsonFile[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <- jsonFile[[i]][["file_name"]]
  filesNameToBarcode <- rbind(filesNameToBarcode,data.frame(filesName = file_name,TCGA_Barcode = TCGA_Barcode))
}
rownames(filesNameToBarcode) <- filesNameToBarcode[,1]

###-----2.Reading miRNA data-----------
filepath <- dir(path ="./CESC-miRNAData",pattern=".mirnas.quantification.txt$",full.names = TRUE,recursive = T)
counts <- data.frame()
rpm <- data.frame()
for(wd in filepath){
  read.table(wd,header = F,sep = "\t")
  oneSampExp <- read.table(wd,header = T)
  tempPath <- unlist(strsplit(wd,"/"))
  filename <- tempPath[length(tempPath)]
  print(paste0("load data:",filename,""))
  oneSampcounts <- oneSampExp[,1:2]
  oneSampRPM <- oneSampExp[,c(1,3)]
  colnames(oneSampcounts) <- c("miRNA_ID",filesNameToBarcode[filename,"TCGA_Barcode"])
  colnames(oneSampRPM) <- c("miRNA_ID",filesNameToBarcode[filename,"TCGA_Barcode"])
  if (dim(counts)[1]== 0){
    counts <- oneSampcounts
    rpm <- oneSampRPM
  }
  else {
    counts <- merge(counts,oneSampcounts,by = "miRNA_ID")
    rpm <- merge(rpm,oneSampRPM,by = "miRNA_ID")
  }
}


#-------3.expression data dispose
DLExpData <- function(Exp){
  Exp <- as.matrix(Exp)
  rownames(Exp) <- Exp[,"miRNA_ID"]
  Exp <- Exp[,-1] %>%  avereps() %>% as.data.frame()
  TumorBarcode <-c()
  NormalBarcode <-c()
  for (id in colnames(Exp)) {
    ids <- str_split(id,'-',simplify = T)[4]
    ids <- as.numeric(substr(ids,1,1))
    
    if(ids == 0){
      TumorBarcode <- c(TumorBarcode,id)
    }else if (ids == 1) {
      NormalBarcode <-c(NormalBarcode,id)
    }
  }
  sampleNum <- data.frame(NormalNum = length(NormalBarcode),
                          TumorNum = length(TumorBarcode))
  sampleList <- data.frame(ID = c(NormalBarcode,TumorBarcode),
                           Type = c(rep("Normal",length(NormalBarcode)),
                                    rep("Tumor",length(TumorBarcode))))
  return(list(Exp = Exp[,c(NormalBarcode,TumorBarcode)],
              BarcodeSort = c(NormalBarcode,TumorBarcode),
              TumorBarcode = TumorBarcode,
              NormalBarcode = NormalBarcode,
              sampleNum = sampleNum ,sampleList = sampleList,
              note = "The column order in Exp is normal sample first,tumor sample last"))
}

miRNACountsExpdata <- DLExpData(Exp = counts)
miRNA_RPM_Expdata <- DLExpData(Exp = rpm)
miRNAdata <- list(miRNACountsExpdata = miRNACountsExpdata,
                  miRNA_RPM_Expdata = miRNA_RPM_Expdata)
save(miRNAdata,file = "miRNAdata.Rdata")

## clinical data dispose------
clinical <- read.delim("~/TCGA/TCGA_CESC/clinical.cart.2020-12-04/clinical.tsv")
clinical <- clinical[!duplicated(clinical$case_id),]
clinical <- clinical[,c("case_submitter_id","age_at_index","days_to_death","vital_status",
                        "ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_t",
                        "days_to_last_follow_up")]
colnames(clinical) <- c('ID','age','days','status','m','n','t','os')
clinical$os <- ifelse(clinical$days >0, clinical$days,clinical$os)
clinical$status <- ifelse(clinical$status=='Alive',0,1)
save(clinical,file = 'CesClinc.Rdata')

# differential expression analysis------
library(edgeR)
load("/Users/llls2012163.com/TCGA/TCGA_CESC/pcExpdata.Rdata")
exp <- pcExpdata[[1]]
x <- rownames(exp)
#### obtain integer------ 
apply(exp,2,class)
exp <- apply(exp,2,as.numeric)
exp <- round(as.matrix(exp))
rownames(exp)<- x
####edgeR statistical analysis------
library(edgeR)
group <- c(rep(1,3),rep(2,240))
y <- DGEList(counts = exp3,group = group)
###filter the expression data------
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.size=FALSE]
###calculate normal factors------
y <- calcNormFactors(y)
#estimate disperassion------
y <- estimateDisp(y)
###exact t test------
et <- exactTest(y)
et <- as.data.frame(et)
#####selections of DEGs: up(p <0.05 & logFC >2) down(p < 0.05 & logFC < -2)------
DEG <- et
DEG$Change <- ifelse(DEG3$logFC >2 & DEG3$PValue <=0.05, 'Up',ifelse(DEG3$logFC < -2 & DEG3$PValue <=0.05,'Down','Non'))
gene_id <- rownames(DEG)
DEG <- cbind(gene_id,DEG)
write.csv(DEG,'RNADEG.csv',row.names = FALSE)# the DEG lncRNAs and miRNAs were obtained similar to above steps
# the results visualizations steps were performed through ggplot2 R package which was omitted from the domain workflow

# univariate Cox regression for detecting prognostic DEGs------
## data preparsion------
library(survival)
library(survminer)
load("/Users/llls2012163.com/TCGA/TCGA_CESC/pcExpdata.Rdata")
load("~/ceRNA/CesClinc.Rdata")
exp <- pcExpdata$Exp
exp <- as.data.frame(t(exp))
exp <- cbind(ID= rownames(exp),exp)
pre_survive_data <- merge(exp,clinical[1,4,8],by.x = 'ID',by.y = "case_submitter_id")
### univariate Cox regression------
gene <- colnames(pre_survive_data)[3:ncol(pre_survive_data)]
gene <- as.vector(gene)
###originate the sheet------
outTab <- data.frame()
for(i in gene){
  expr=pre_survive_data[,i]
  cox=coxph(Surv(Day,statu)~ expr,pre_survive_data)
  coxsummary=summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,'exp(coef)'],2),
                            ###HR hazard rate
                            z=round(coxsummary$coefficients[,'z'],2),##z value
                            '95%CI'=paste(round(coxsummary$conf.int[,3],2),
                                          round(coxsummary$conf.int[,4],2),sep = '-'),##95%CI
                            pvalue=round(coxsummary$coefficients[,'Pr(>|z|)'],2)))
  ##p value
}# the steps to detecting the prognostic miRNAs and prognostic lncRNAs similar to above steps and 
# the forest plot for showing this results were omitted from the main workflow

# gene functional annotations: here we screen prognostic mRNAs as target genes------
library(clusterProfiler)
library(org.Hs.eg.db)
gene <- bitr(targetGenes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db) # id transform
GOenrich <- enrichGO(gene$ENTREZID,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',pvalueCutoff = 1,
                     pAdjustMethod = "BH",minGSSize = 4,maxGSSize = 100,ont = 'ALL') # GO annotations
save(GOenrich,file = 'Goenrich.Rdata') # results saving. The results visulization steps were omitted from pipeline to maintain simplization of codes.
# Prognostic model development------
## we first filtered the key prognostic mRNAs through a multivariate Cox regression------ 
library(survival)
library(survminer)
load("/Users/llls2012163.com/ceRNA/CesClinc.Rdata")
# mapping the target genes from pcRNA expression matrix and merge it with its clinical data
s <- t(exp)['target genes',]
s <- merge(clinical,s,by.x='ID',by.y='ID')
expr <- paste('target genes','+')
cox=coxph(Surv(Day,statu)~ expr,s)# here we abtained 4 key genes for prognostic model development
# Risk Score (RS) = (ExpOPN3 × βOPN3) + (ExpDAAM2 × βDAAM2) + (ExpHENMT1 × βHENMT1) + (ExpCAVIN3 × βCAVIN3) 
# Then an multivariate Cox regression were used to identify the individualization of our model. Code as above.


