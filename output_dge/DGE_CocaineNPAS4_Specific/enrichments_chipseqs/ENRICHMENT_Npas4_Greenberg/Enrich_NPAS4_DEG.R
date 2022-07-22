library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)

# Brain Exp
rm(list=ls())
wgcna = list.files(pattern = 'MouseID')
tab=read.table(wgcna,sep="\t",header=T)
tab <- tab[,1:2]
colnames(tab)=c("Gene","DEFINITION")
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
load("Greenberg_GeneSets.RData")
#GeneSets <- GeneSets[1:7]
#for(i in 1:length(GeneSets))
#{
#	GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
#}
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 21914-Genes$Freq-nrow(GeneSets[[i]]) #Protein Coding bkg = 19901 Ortho bkg = 14373 Brain Expressed bkg = 15585
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
save(FisherMat,file="GeneSet_NPAS4_Enrichment.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='BH') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames


FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0

p <- as.data.frame(FisherAdj)
cor <- as.data.frame(FisherOR)


p $Module <- rownames(p )
cor$Module <- rownames(cor)

p <- melt(p)
cor <- melt(cor)
p$cor <- cor$value
p$log <- -log10(p$value)

p$log[p$log < 1.3] <- NA
p$OR<- ifelse(is.na(p$log), p$log, p$cor)

# Colors
low0 <- rgb(246, 248, 251, maxColorValue = 255)
low1 <- rgb(201, 230, 234, maxColorValue = 255)
low2 <- rgb(89, 158, 193, maxColorValue = 255)
mid <- rgb(67, 121, 180, maxColorValue = 255)
high1 <- rgb(65, 90, 158, maxColorValue = 255)
high2 <- rgb(52, 52, 106, maxColorValue = 255)
high3 <- rgb(15, 15, 25, maxColorValue = 255)

pdf("Enrichment_Npas4_Cell_Bubble.pdf",width=5,height=6)
ggscatter(p, 
            x = "variable",
            y = "Module",
            size="OR",
            color="log",
            alpha = 0.8,
            xlab = "",
            ylab = "") +
            theme_minimal() +
            rotate_x_text(angle = 45)+
        scale_size(range = c(0, 10))+ 
        gradient_color(c("lightblue","darkblue"))
dev.off()

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
#order <- c("min30","hr01","hr03","hr06")
#FisherAdj <- FisherAdj[match(order,rownames(FisherAdj)),]
#FisherOR <- FisherOR[match(order,rownames(FisherOR)),]
pdf("Enrichment_Npas4_Cell.pdf",width=5,height=6)
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
par(mar = c(8, 8, 5, 5));
labeledHeatmap(Matrix = df, xLabels = colnames(df), yLabels = rownames(df), colorLabels =FALSE,colors=colorRampPalette(c(low0, low1, low2, mid, high1))(50),textMatrix=LabelMat, setStdMargins = FALSE, cex.text = 0.5)
dev.off()

