# 
# #1. If you have separate vcf files for each chromosome, first get the union of calls made on the same samples (i.e. Merge all VCFs files) using GATK. It is easier to work on one file for Principal componenet analysis.
# #Here, reference.fasta should be the same file that was used to create vcf file (if same reference.fasta is not available, hg19.fasta can be corrected by keeping the same chromosomes and contigs included in the vcf file). Check GATK command for CombineVariants for details. This also requires the use of picard tool to create dictionary and index file (http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference) before executing the commands below. Then execute the GATK command to merge multiple VCF files as shown below:
# java -jar /home/an/Desktop/qc_paper/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
# -T CombineVariants \
# -R /home/an/Desktop/qc_paper/reference2.fasta \
# --variant:chr1 AOGC-Genotyping.output.chr1.recalibrated.filtered.vcf \
# --variant:chr2 AOGC-Genotyping.output.chr2.recalibrated.filtered.vcf \
# --variant:chr3 AOGC-Genotyping.output.chr3.recalibrated.filtered.vcf \
# --variant:chr4 AOGC-Genotyping.output.chr4.recalibrated.filtered.vcf \
# --variant:chr5 AOGC-Genotyping.output.chr5.recalibrated.filtered.vcf \
# --variant:chr6 AOGC-Genotyping.output.chr6.recalibrated.filtered.vcf \
# --variant:chr7 AOGC-Genotyping.output.chr7.recalibrated.filtered.vcf \
# --variant:chr8 AOGC-Genotyping.output.chr8.recalibrated.filtered.vcf \
# --variant:chr9 AOGC-Genotyping.output.chr9.recalibrated.filtered.vcf \
# --variant:chr10 AOGC-Genotyping.output.chr10.recalibrated.filtered.vcf \
# --variant:chr11 AOGC-Genotyping.output.chr11.recalibrated.filtered.vcf \
# --variant:chr12 AOGC-Genotyping.output.chr12.recalibrated.filtered.vcf \
# --variant:chr13 AOGC-Genotyping.output.chr13.recalibrated.filtered.vcf \
# --variant:chr14 AOGC-Genotyping.output.chr14.recalibrated.filtered.vcf \
# --variant:chr15 AOGC-Genotyping.output.chr15.recalibrated.filtered.vcf \
# --variant:chr16 AOGC-Genotyping.output.chr16.recalibrated.filtered.vcf \
# --variant:chr17 AOGC-Genotyping.output.chr17.recalibrated.filtered.vcf \
# --variant:chr18 AOGC-Genotyping.output.chr18.recalibrated.filtered.vcf \
# --variant:chr19 AOGC-Genotyping.output.chr19.recalibrated.filtered.vcf \
# --variant:chr20 AOGC-Genotyping.output.chr20.recalibrated.filtered.vcf \
# --variant:chr21 AOGC-Genotyping.output.chr21.recalibrated.filtered.vcf \
# --variant:chr22 AOGC-Genotyping.output.chr22.recalibrated.filtered.vcf \
# --variant:chrX AOGC-Genotyping.output.chrX.recalibrated.filtered.vcf \
# --variant:chrY AOGC-Genotyping.output.chrY.recalibrated.filtered.vcf \
# --variant:chrM AOGC-Genotyping.output.chrM.recalibrated.filtered.vcf \
# -o AOGC-Genotyping.output.chrALL.recalibrated.filtered.vcf \
# -genotypeMergeOptions PRIORITIZE \
# -priority chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM
# 
# ##2. Now using vcftools, use merged vcf file to write plink files (.ped and .map files) 
# vcftools --vcf AOGC-Genotyping.output.chrALL.recalibrated.filtered.vcf --remove-filtered-geno-all --plink  --out chrALL

##3. The .map file thus created could have rs ids which needs to be fixed. The map files shouls have same ids as in HapMap file. Here we use 'chr' and 'position' to match with HapMap SNPs. The rscript below can be used to fix the rs id position in .map file: 
newqcvcf <- read.delim("/home/an/Desktop/qc_paper/chrALL.map",header=FALSE,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
#newqcvcf <- read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/_Collaborations/Neupane_Leo_Collaboration/UQCCG/all_bams_analysis/vcf_create/SNPs/per_chromosome/gwas/qcvcffile_copy.txt",header=FALSE,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
colnames(newqcvcf)<-c("chr","rsid","centi","position")
newqcvcf[1:10,]
newqcvcf[,"chr"] <- paste0("chr",newqcvcf[,"chr"])
newqcvcf[,"rsid"] <- paste(newqcvcf[,"chr"],newqcvcf[,"position"],sep=":")
newqcvcf[1:10,]
newqcvcf[,"chr"] <- gsub("chr", "",newqcvcf[,"chr"])
newqcvcf[1:10,]
write.table(newqcvcf,"/home/an/Desktop/qc_paper/chrALL.map",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE) #This is the corrected map file 

# ##4. Now, work with plink tools (http://zzz.bwh.harvard.edu/plink/download.shtml) to convert plink files into bed format. These files will be spiked with the bed files for HapMap data.
# ##see http://zzz.bwh.harvard.edu/plink/reference.shtml and http://zzz.bwh.harvard.edu/plink/data.shtml for details
# plink --file chrALL --make-bed --out qcvcf
# #Next MAF and missingness filters are set
# plink --bfile qcvcf  --maf 0.02 --geno 0.1 --make-bed --out qcvcf.f   
# #Check frequencies for the selected SNPs and founders from your data
# plink --bfile qcvcf.f  --freq  
# #Check frequencies for the selected SNPs and founders from HapMap data
# plink --bfile ALL_eth_650Y_forward_hg19_Best  --freq
# 
# #########################################################
##5. Now, run this rscript to get the list of common SNPs in HapMap and your data
code.dir <- "/home/an/Desktop/qc_paper" #set directory that contains both .bim files 
setwd(code.dir)
bim.file.list <- c("ALL_eth_650Y_forward_hg19_Best.bim","qcvcf.f.bim") ## list all bim files to get common snp basis

num.vars<-6; skip.lines <- 0
all.snps <- {}
# i <- 2
combined.bim <- {}
for(i in 1:length(bim.file.list)){

a.bim <- try(scan(bim.file.list[i],what=character(num.vars),skip=skip.lines,fill=TRUE))
num.lines <- length(a.bim)/(num.vars)
dim(a.bim) <- c(num.vars,num.lines)
a.bim <- t(a.bim)
if(sum(grepl("^RS",a.bim[,2]))>0){
  print(paste("warning",bim.file.list[i]," contains RS NOT rs ids use -sed -i '/RS/rs/' file.bim . Converting...",sep=" "))
  a.bim[,2] <- gsub("^RS","rs",a.bim[,2])
      }
all.snps <- c(all.snps,a.bim[,2])
a.bim[1:5,]
if(is.null(dim(combined.bim))){ combined.bim <- a.bim}else{combined.bim <- rbind(combined.bim,a.bim)}
}

counts <- tapply(all.snps,all.snps,length)
in.common <- counts==length(bim.file.list)
print(paste(sum(in.common)," snps in common for all bims"))

colnames(combined.bim) <- c("chr","SNP","CM","start","A1","A2")
combined.bim.common <- combined.bim[ combined.bim[,"SNP"] %in% names(counts[in.common]) ,]
dim(combined.bim.common)

order.by <- order(combined.bim.common[,"SNP"])
combined.bim.common <- combined.bim.common[order.by,]
combined.bim.common[1:6,]

unique.snps <- unique(combined.bim.common[,"SNP"])
dups <- duplicated(combined.bim.common[,"SNP"])
dups[1:5]

table.common <- cbind(combined.bim.common[!dups,],combined.bim.common[dups,c("A1","A2")])
table.common[1:5,]
test <- apply(table.common[,5:8],1,function(x) length(unique(x)) )
test[1:10]
names(test) <- table.common[,"SNP"]

bad.alleles <- names(test)[test>2]

keep <- names(counts[in.common])[!(names(counts[in.common]) %in% bad.alleles)]
length(keep)
combined.bim[,"SNP"]
write.table(keep,"common_snps_Chip.650.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

# #################################################### End of rscript to get the list of common SNPs
# 
# ##6. Now, use plink tool to extract coommon snps
# plink --bfile ALL_eth_650Y_forward_hg19_Best  --extract common_snps_Chip.650.txt --make-bed --out  ALL_eth_650Y_forward_hg19_Best.C
# 
# plink --bfile qcvcf.f  --extract common_snps_Chip.650.txt --make-bed --out  qcvcf.f.C
# 
# plink --bfile  ALL_eth_650Y_forward_hg19_Best.C --bmerge qcvcf.f.C.bed qcvcf.f.C.bim qcvcf.f.C.fam --allow-no-sex --make-bed  --out data_WTCCC_f_650 ### merger here
# ##8. extract long range LD positions using plink tool
# plink --bfile data_WTCCC_f_650   --exclude Price2008_hg19_paul.txt --range --make-bed --out data_WTCCC_f_650 --noweb
# 
# ##9. Perform principal component analysis. 
# #You can Use any one of PCA tools as demonstrated below:
# #Use shellfish.py to run PCA (download from here: http://www.stats.ox.ac.uk/~davison/software/shellfish/shellfish.php)
# #Note: all data_WTCCC_f_650 files created in step #8 should be in /shellfish directory
# /shellfish/shellfish.py --pca --numpcs 10 --maxprocs 8 --file  data_WTCCC_f_650 --out qcvcf --ignore-sge  
# 
# ##or use plink2 (http://www.cog-genomics.org/plink/2.0) to run PCA
# plink2 --bfile data_WTCCC_f_650 --pca approx 10 --out qcvcf
# 
# #or use eigenstrat
# ##Note: whichever method you use for PCA as mentioned above, I have written a snippet of code to process the PCA file in 'Cleaned_plot_650_strat_with_dist_with_euclidian_distance_4d.r', and you have to select PCA type in the code. You can also use your choice of PCA tool and customize the rscript 'Cleaned_plot_650_strat_with_dist_with_euclidian_distance_4d.r' to PCA data. 
# 
# ##10. now use plot_650_strat_with_dist.r script provided