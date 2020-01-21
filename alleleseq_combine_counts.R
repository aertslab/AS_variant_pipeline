


hap1_unique_counts<-read.delim(file="HET_SNVs_hap1_counts",header=F)
hap2_unique_counts<-read.delim(file="HET_SNVs_hap2_counts",header=F)
hap1_common_counts<-read.delim(file="HET_SNVs_hap1_common_counts", header=F)
hap2_common_counts<-read.delim(file="HET_SNVs_hap2_common_counts", header=F)
hap1_and_hap2_common_counts<-merge(hap1_common_counts,hap2_common_counts,by.x="V2",by.y="V2",all=T)


hap1_and_hap2_common_counts$V1.x<-as.character(hap1_and_hap2_common_counts$V1.x)
hap1_and_hap2_common_counts$V1.y<-as.character(hap1_and_hap2_common_counts$V1.y)


to_subtract<-subset(hap1_and_hap2_common_counts,V1.x!=V1.y)
# we might get a multiple allelic compositions for a single position (liftover artifact) - we filter out positions like this
hap2_common_counts_subtracted<-hap2_common_counts[! hap2_common_counts$V2 %in% to_subtract$V2,]

tmp_merg<-merge(hap1_unique_counts,hap2_unique_counts,by.x="V1",by.y="V1",all=T)
colnames(tmp_merg)<-c("SNV_ID","HAP1","HAP2")
merged_counts<-merge(tmp_merg,hap2_common_counts_subtracted, by.x="SNV_ID",by.y="V2",all=T)
merged_counts<-merged_counts[!duplicated(merged_counts$SNV_ID),]
rownames(merged_counts)<-merged_counts$SNV_ID
merged_counts<-merged_counts[,-1]
merged_counts$HAP1<-as.character(merged_counts$HAP1)
merged_counts$HAP2<-as.character(merged_counts$HAP2)
merged_counts$V1<-as.character(merged_counts$V1)
merged_counts[is.na(merged_counts)] <- "0_0_0_0"


out<-strsplit(as.character(merged_counts$HAP1),"_")
a<-do.call(rbind,out)
hap1_unique_counts<-apply(a,2,as.numeric)

out<-strsplit(as.character(merged_counts$HAP2),"_")
a<-do.call(rbind,out)
hap2_unique_counts<-apply(a,2,as.numeric)

out<-strsplit(as.character(merged_counts$V1),"_")
a<-do.call(rbind,out)
common_counts<-apply(a,2,as.numeric)

combined_counts<-hap1_unique_counts+hap2_unique_counts+common_counts
rownames(combined_counts)<-rownames(merged_counts)
colnames(combined_counts)<-c("nA","nC","nG","nT")

heterozygote_SNVs_in_peaks<-read.delim(file="heterozygote_SNVs_in_peaks.BL_filtered.bed",header=F)
colnames(heterozygote_SNVs_in_peaks)<-c("chr","start","end","ID","REF","ALT","DP","AD","GT","AF","peak_chr","peak_start","peak_end","peak_name","peak_score","SNV_ID")

m<-merge(heterozygote_SNVs_in_peaks,combined_counts,by.x="SNV_ID", by.y="row.names")

m$nTotal<-m$nA+m$nC+m$nG+m$nT

 dim(subset(m,nTotal>=6))

# count the number of reference and variant allele counts
for(i in 1:dim(m)[1] ) {
  if(m$REF[i] == "A" ) {
    m$nRef[i] = m$nA[i]
} else if (m$REF[i] == "C") {
  m$nRef[i] = m$nC[i]
} else if (m$REF[i] == "G") {
  m$nRef[i] = m$nG[i]
}
else if (m$REF[i] == "T") {
  m$nRef[i] = m$nT[i]
}
}

for(i in 1:dim(m)[1] ) {
  if(m$ALT[i] == "A" ) {
    m$nAlt[i] = m$nA[i]
} else if (m$ALT[i] == "C") {
  m$nAlt[i] = m$nC[i]
} else if (m$ALT[i] == "G") {
  m$nAlt[i] = m$nG[i]
}
else if (m$ALT[i] == "T") {
  m$nAlt[i] = m$nT[i]
}
}

# add info if it's a winning or losing situation

for(i in 1:dim(m)[1] ) {
    if (m$GT[i] == "1|0" && m$nAlt[i] > m$nRef[i]) { m$WIN[i] ="HAP1_WINS" }
    else if (m$GT[i] == "1|0" && m$nAlt[i] < m$nRef[i]) { m$WIN[i] ="HAP1_LOSES" }
    else if (m$GT[i] == "1|0" && m$nAlt[i] == m$nRef[i]) { m$WIN[i] ="DRAW" }
    else if (m$GT[i] == "0|1" && m$nAlt[i] > m$nRef[i]) { m$WIN[i] ="HAP2_WINS" }
    else if (m$GT[i] == "0|1" && m$nAlt[i] < m$nRef[i]) { m$WIN[i] ="HAP2_LOSES" }
    else if (m$GT[i] == "0|1" && m$nAlt[i] == m$nRef[i]) { m$WIN[i] ="DRAW" }}
m_min6<-subset(m,nTotal>=6)

for (i in 1:dim(m_min6)[1]) {
m_min6$pval[i]<-binom.test(m_min6$nAlt[i],m_min6$nTotal[i])$p.value }

m_min6$padj<-p.adjust(m_min6$pval, method="fdr")

write.table(file="m_min6.txt", m_min6, row.names=F, sep="\t", quote=F)
save(file="m_min6.RData",m_min6)


m_min10<-subset(m,nTotal>=10)

for (i in 1:dim(m_min10)[1]) {
m_min10$pval[i]<-binom.test(m_min10$nAlt[i],m_min10$nTotal[i])$p.value }

m_min10$padj<-p.adjust(m_min10$pval, method="fdr")

print("Number of significant events at min-read count 6 (padj<0.05)")
dim(subset(m_min6,padj<=0.05))

print("Number of significant events at min-read count 6 (padj<0.2)")
dim(subset(m_min6,padj<=0.2))

save(file="m_min10.RData",m_min10)

write.table(file="m_min10.txt", m_min10, row.names=F, sep="\t", quote=F)


# plot allelic-ratios
m_min6$VAF<-m_min6$nAlt/m_min6$nTotal
m_min10$VAF<-m_min10$nAlt/m_min10$nTotal

png(file="allelic_ratio_m_min6.png",width=800, height=800, res=150, type="cairo")
hist(m_min6$VAF, breaks = 17, col=rgb(0.5,0.5,0.5,0.5),xlab="VAF", main="allelic ratios with m_min 6")
hist(subset(m_min6,padj<=0.05)$VAF, breaks = 17, col=rgb(1,0,0,0.5),add=T)
dev.off()

png(file="allelic_ratio_m_min10.png",width=800, height=800, res=150, type="cairo")
hist(m_min10$VAF, breaks = 17, col=rgb(0.5,0.5,0.5,0.5),xlab="VAF", main="allelic ratios with m_min 10")
hist(subset(m_min10,padj<=0.05)$VAF, breaks = 17, col=rgb(1,0,0,0.5),add=T)
dev.off()
