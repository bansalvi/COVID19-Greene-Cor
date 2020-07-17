#Vikas Bansal
#Make a correlation of log2fc for cell lines and human samples from tenover and landthaler paper

set.seed(786)
setwd("/data/vikas/COV_ace/")


library(ggplot2)
library(cowplot)

library(dplyr)
library(readxl)

library(stringi)
library(purrr)
library(reshape)
library(corrplot)

#pathto.outPlots <- "/data/vikas/COV_ace/OutputPlots"

HumanCellLines <- read_excel("PublicData/DEgenesTenOverLandthaler/1-s2.0-S009286742030489X-mmc1.xlsx", sheet = 2)
HumanCellLinesv2 <- (as.data.frame(HumanCellLines))




NHEB <- read_excel("PublicData/DEgenesTenOverLandthaler/1-s2.0-S009286742030489X-mmc2.xlsx", sheet = 2)
NHEBv2 <- (as.data.frame(NHEB))


NHEB_HCL <- merge(NHEBv2[,1:2],HumanCellLinesv2[,1:10], by="GeneName")


COVID19 <- read_excel("PublicData/DEgenesTenOverLandthaler/1-s2.0-S009286742030489X-mmc4.xlsx", sheet = 2)
COVID19v2 <- (as.data.frame(COVID19))

colnames(COVID19v2)[1] <- "GeneName"
COVID19v3 <- COVID19v2[(COVID19v2$status == "OK"),]
COVID19v3$`SARS-CoV-2(COVID19)_L2FC` <- as.numeric(COVID19v3$`SARS-CoV-2(COVID19)_L2FC`)

NHEB_HCL_COVID <- merge(COVID19v3[,c(1,3)],NHEB_HCL, by="GeneName")

NHEB_HCL_COVIDtmp <- NHEB_HCL_COVID[(NHEB_HCL_COVID[,1]%in%(COVID19v3[(COVID19v3$`SARS-CoV-2(COVID19)_L2FC` < 0),1])),]




Calu3_S2_24h <- read.delim("PublicData/DEgenesTenOverLandthaler/DE_Calu3_S2.24h_mock.24h.txt")
colnames(Calu3_S2_24h)[1:2] <- c("GeneName","SARS-CoV-2(Calu-3 24hrs)_L2FC")
Calu3_S2_12h <- read.delim("PublicData/DEgenesTenOverLandthaler/DE_Calu3ser2_S2.12h_mock.12h.txt")
colnames(Calu3_S2_12h)[1:2] <- c("GeneName","SARS-CoV-2(Calu-3 12hrs)_L2FC")
Caco2_S2_24h <- read.delim("PublicData/DEgenesTenOverLandthaler/DE_Caco2_S2.24h_mock.24h.txt")
colnames(Caco2_S2_24h)[1:2] <- c("GeneName","SARS-CoV-2(Caco-2 24hrs)_L2FC")

LandthalerCellsS2 <- merge(Caco2_S2_24h[,1:2],(merge(Calu3_S2_24h[,1:2],Calu3_S2_12h[,1:2], by="GeneName")),by="GeneName")



Calu3_S1_24h <- read.delim("PublicData/DEgenesTenOverLandthaler/DE_Calu3_S1.24h_mock.24h.txt")
colnames(Calu3_S1_24h)[1:2] <- c("GeneName","SARS-CoV-1(Calu-3 24hrs)_L2FC")
Calu3_S1_12h <- read.delim("PublicData/DEgenesTenOverLandthaler/DE_Calu3ser2_S1.12h_mock.12h.txt")
colnames(Calu3_S1_12h)[1:2] <- c("GeneName","SARS-CoV-1(Calu-3 12hrs)_L2FC")
Caco2_S1_24h <- read.delim("PublicData/DEgenesTenOverLandthaler/DE_Caco2_S1.24h_mock.24h.txt")
colnames(Caco2_S1_24h)[1:2] <- c("GeneName","SARS-CoV-1(Caco-2 24hrs)_L2FC")

LandthalerCellsS1 <- merge(Caco2_S1_24h[,1:2],(merge(Calu3_S1_24h[,1:2],Calu3_S1_12h[,1:2], by="GeneName")),by="GeneName")


LandthalerCells <- merge(LandthalerCellsS1,LandthalerCellsS2,by="GeneName")

NHEB_HCL_COVID_Landthaler <- merge(LandthalerCells,NHEB_HCL_COVID,by="GeneName")


M<-cor(NHEB_HCL_COVID_Landthaler[,-1], method = "pearson")


p.mat <- cor.mtest(NHEB_HCL_COVID_Landthaler[,-1])

col=brewer.pal(n=10, name="RdYlBu")

pdf("PublicData/DEgenesTenOverLandthaler/LandthalertenOverCorFC.pdf", width=20, height=18)
corrplot(M, method="color", col=col,
          order="hclust",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", number.cex=1.5,cl.cex=2,
         # hide correlation coefficient on the principal diagonal
         diag=FALSE
)
dev.off()
