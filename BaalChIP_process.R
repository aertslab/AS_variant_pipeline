.libPaths("/staging/leuven/stg_00002/lcb/zkalender/R/3.6")

library(BaalChIP)
library(doParallel)
library(coda)
library(scales)
library(plotrix)

args <- commandArgs(trailingOnly = TRUE)
sample_name<-args[1]

work_dir<-paste(paste("/staging/leuven/stg_00002/lcb/zkalender/melanoma_WGS/",sample_name,sep=""),"/alleleseq_bowtie2/BaalChIP",sep="")

setwd(work_dir)


m_min6_with_RAF<-read.delim(file="m_min6_with_RAF_for_BaalChIP",header=F)
colnames(m_min6_with_RAF)<-c("chr","start","end","ID","nTotal","ATAC.REF.1","ATAC.ALT.1","chr_WGS","start_WGS","end_WGS","WGS_DP","WGS_VAF")


 counts<-m_min6_with_RAF[,c("ID","ATAC.REF.1","ATAC.ALT.1")]
counts$ATAC.score<-rep(1,dim(counts)[1])
counts<-counts[,c(1,4,2,3)]

m_min6_with_RAF$RMbias<-0.5
 m_min6_with_RAF$RAF<-1-m_min6_with_RAF$WGS_VAF
biastable<-m_min6_with_RAF[,c("ID","RMbias","RAF")]

Iter=5000
conf_level = 0.95
cores=16
 source("/staging/leuven/stg_00002/lcb/zkalender/R/3.6/BaalChIP/R/applyBayes.R")


if (nrow(counts) > 0) {
           Bayes_report <- runBayes(counts = counts, bias = biastable, Iter = Iter, conf_level = conf_level,
               cores = cores)
       } else {
           message("no variants left for ", ID)
           Bayes_report <- data.frame()
       }


Bayes_report$Corrected.AR<-rowMeans(Bayes_report[,c("Bayes_lower","Bayes_upper")])
Bayes_report$isASB<-ifelse(Bayes_report$Bayes_sig_A==1 | Bayes_report$Bayes_sig_B==1, "TRUE" , "FALSE")
Bayes_report$RAF<-m_min6_with_RAF$RAF

Bayes_report$AR<-m_min6_with_RAF$ATAC.REF.1/m_min6_with_RAF$nTotal

save(file="Bayes_report.RData", Bayes_report)
write.table(file="Bayes_report.txt",Bayes_report, sep="\t",quote=F,row.names=F)

reg<-lm(Bayes_report$RAF~Bayes_report$Corrected.AR)
cor_res<-cor.test(Bayes_report$RAF,Bayes_report$Corrected.AR,method = "spearman")
png(file=paste(sample_name,".Corrected.AR_versus_RAF.png",sep="") ,width=1200,height=1200, res=150,type="cairo")
plot(Bayes_report$RAF~Bayes_report$Corrected.AR,pch=20,col=alpha("black",0.5) ,cex=2,xlim=c(0,1),ylim=c(0,1),xlab="Corrected AR", ylab="RAF",bty="n", main=sample_name)
ablineclip(reg, col="blue",lwd=3, x1=0, x2=1)
mtext(paste("p-value:",round(cor_res$p.value,digits=3),"cor:",round(cor_res$estimate,digits=3),sep=" "),3,3)
dev.off()

png(file="allelic_ratio_hist.png",width=1200, height=800, res=150, type="cairo")
hist(Bayes_report$AR, breaks = 17, col=rgb(0.5,0.5,0.5,0.5),xlab="AR", main=sample_name)
dev.off()


png(file="corrected_allelic_ratio_hist.png",width=1200, height=800, res=150, type="cairo")
hist(Bayes_report$Corrected.AR, breaks = 17, col=rgb(0.5,0.5,0.5,0.5),xlab="Corrected AR", main=sample_name)
hist(subset(Bayes_report,isASB=="TRUE")$Corrected.AR, breaks = 17, col=rgb(1,0,0,0.5),add=T)
dev.off()
