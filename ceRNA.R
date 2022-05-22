s <- rbind(TRP_series$TRPV,TRP_series$TRPC,TRP_series$TRPM,TRP_series$TRPA1)
s <- as.data.frame(s)
s[,2:ncol(s)] <- apply(s[,2:ncol(s)],2,as.numeric)
s <- s[,-1]
s <- t(s)
clin <- ICI_logRank$data[,-2]
m <- as.data.frame(`geneClusterTcga.k=3.consensusClass`)
colnames(m) <- c('ID','Group')
m$Group <- ifelse(m$Group==1,'A',
                  ifelse(m$Group==2,'B','C'))
data <- merge(clin,m,by.x = 'ID',by.y = 'ID')
library(survminer)
library(survival)
library(ggprism)
b <- data
fit <- survfit(Surv(Day, Statu) ~ group, data = b)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_prism(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF",'red')
)


###
s <- rbind(TRP_series$TRPV,TRP_series$TRPC,TRP_series$TRPM,TRP_series$TRPA1)
s1 <- s[,-1]
s1 <- apply(s1,2,as.numeric)
rownames(s1) <- s$ID
s1 <- as.data.frame(s1)
FpkmToTPM <- function(fpkm){
  exp(log(fpkm)-log(sum(fpkm))+log(1e6))
}###Function of FPKM2TPM
s2 <- apply(s1,2,FpkmToTPM)
d=s2

w <- as.data.frame(`geneClusterTcga.k=3.consensusClass`)
colnames(w) <- c('ID','group')
w$group <- ifelse(w$group==1,'A',
                  ifelse(w$group==2,'B','C'))
w1 <- ICI_logRank$data[,-2]
data <- merge(w,w1,by.x = 'ID',by.y = 'ID')

s <- as.data.frame(z)
s < as.matrix(z)
