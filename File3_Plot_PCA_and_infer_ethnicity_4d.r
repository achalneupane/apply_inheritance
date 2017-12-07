##This script is used to plot PCA and also to infer ethnicity of samples. The ethnicity of samples are parsed in "RELATED" IBD output file obtained after running IBS calculation (ApplyinheritanceQC) code provided.

##Choose either "shellfish" or "plink2" or "other" tools you have used to generate pca. Can use any other PCA tools provided it takes plink format files. You may also need to re-build your data structure for 'vec' below if you decide to use some other tools for pca.
#used.pca.tool <- "other" #may need further processing of data below
used.pca.tool <- "shellfish"  
#used.pca.tool <- "plink2"
#Choose directory where "data_WTCCC_f_650.fam", "qcvcf.evecs", including all the required files are present.
ann.dir <- "/home/an/Desktop/qc_paper"
#Choose pca file correctly. For example, "qcvcf.evecs" is from shellfish; "qcvcf.evecs_pca_eigenvectors.txt" is from plink2
plot.pca.file <- "qcvcf.evecs" # "qcvcf.evecs_pca_eigenvectors.txt"
#fam.file<-"exomeChip.GWAS.HBM.650Y.strat.FINAL.fam"  ### fam files used with PCA analysis from plink2 or shellfish or any other tools used.
fam.file <- "data_WTCCC_f_650.fam"

##Choose directory that contains "ann.650Y.file" (i.e. "Sample_Detail_All_ethnicity.txt"), "Generalized_Ethnicity_definition.csv" and "qcvcf.evecs". If "plot.ann.file" is available should be in this directory as well
ann.dir <- "/home/an/Desktop/qc_paper" 
ann.650Y.file <- "Sample_Detail_All_ethnicity.txt"
plot.ann.file <- {} # leave as {} if set to default otherwise provide tab delimited file with column names "Sample","status" information
####### Read IBS file with maximum number of SNPs and with "RELATED" extension in file name which was obtained after running IBD code. 
are.related <- read.delim("/home/an/Desktop/qc_paper/521932_unfiltered_84574_filtered__SNPs_Vcf_Merge.GQ-20.Min.ALT-10_AOGC-Genotyping.output.ARE_RELATED.Mon_Sep_05_2016.genetic_QC.txt",header=T,sep="",fill=TRUE,stringsAsFactors=FALSE)

path.save <- "/home/an/Desktop/qc_paper"  ##choose path where you want to save the output files such as pca plot and inferred ethnicity data

color.with <- "Ethnicity"

eig1 <- "e.1"
eig2 <- "e.2"
eig3 <- "e.3"

################ get PCA loaded with 650Y samples #######

if(exists("vec")){rm(vec)}

## If you decide to you some other tool for pca you may have to work on this bit
if(used.pca.tool == "other"){
  print ("WARNING!! You used other tool for PCA, you may have to check again or re-structre your 'vec' to process further!!")
  vec<-read.delim(paste(ann.dir,plot.pca.file,sep="/"), header=T,sep="",fill=TRUE,stringsAsFactors=FALSE)  
  rownames(vec)<-unlist(lapply(strsplit(rownames(vec), split=":"),function(x) x[1]))   # fix colnum names
  colnames(vec)<-c(paste("e.",1:10,sep=""), "status")   #,"origin","color","points")
  
  }else if(used.pca.tool == "plink2"){
  print ("USING PCA from plink2")
  vec <- read.delim(paste(ann.dir,plot.pca.file,sep="/"), header=T,sep="",fill=TRUE,stringsAsFactors=FALSE)
  rownames(vec) <- vec[,"X.FID"]
  vec <- vec[,c( "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
  
  }else if (used.pca.tool == "shellfish") {
  print ("USING PCA from shellfish")
  vec<-read.delim(paste(ann.dir,plot.pca.file,sep="/"),header=F,sep="",fill=TRUE,stringsAsFactors=FALSE)  
  vec<-t(vec)
  }

  fam<-read.delim(paste(ann.dir,fam.file,sep="/"),header=F,sep="",fill=TRUE,stringsAsFactors=FALSE)
  
  fam[fam[,6]==2,6] <- "Control"
  fam[fam[,6]==1,6] <- "Case"
  
  if(dim(fam)[1] != dim(vec)[1]){
  print("Error fam and evec file had different number of samples that MUST BE the same if shellfish was used")}
  
  if(sum(fam[,1] != rownames(vec_plink))> 0){
  print("Error fam and evec file had different number of samples that MUST BE the same if shellfish was used")}
  vec <- as.matrix(vec)
  rownames(vec) <- fam[,1]   
  vec <- cbind(vec,fam[,6])
  colnames(vec) <- c(paste("e.",1:10,sep=""),"status")
  
#   > vec[1:3,]
#   e.1        e.2         e.3        e.4        e.5        e.6         e.7         e.8         e.9       
#   1percent  "0.020551" "-0.021417" "0.0048"   "0.001146" "0.091746" "-0.033781" "-0.227411" "-0.251145" "0.054517"
#   5percent  "0.020973" "-0.021318" "0.004452" "0.001047" "0.094519" "-0.034128" "-0.230977" "-0.252499" "0.055655"
#   10percent "0.021189" "-0.020907" "0.004242" "0.001162" "0.094721" "-0.033001" "-0.227096" "-0.248857" "0.054731"
#   e.10       status
#   1percent  "0.009821" "-9"  
#   5percent  "0.010094" "-9"  
#   10percent "0.00962"  "-9"  
  
  pca.for.logistic <- cbind(rownames(vec),rownames(vec),vec[,c(paste("e.",1:10,sep=""))])
  colnames(pca.for.logistic) <- c("FID","IID",c(paste("PCA",1:10,sep="")))
  paste(plot.pca.file,"pca_eigenvectors.txt",sep="_")
  getwd()
  setwd(ann.dir)
  write.table(pca.for.logistic,file=paste(plot.pca.file,"pca_eigenvectors.txt",sep="_"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

####################

if(!is.null(plot.ann.file)){
  plot.ann <- read.delim(paste(ann.dir,plot.ann.file,sep="/"),header=T,sep="",fill=TRUE) ### annotation file pca.evec samplea not in 650Y
}else{
  if(!("status" %in% colnames(vec))){ vec[,"status"]<-"unknown"}
  plot.ann<-cbind(rownames(vec),as.character(vec[,"status"]))
  colnames(plot.ann) <- c("Sample","status")
}
####################

ann.650Y <- read.delim(paste(ann.dir,ann.650Y.file,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ### annotation file for 650 Y
ann.650Y[1:5,]
dim(ann.650Y)

## Now combine sample annotation file this has the 650Y and the pca.evec samples.

dim(plot.ann)
dim(ann.650Y) 


sample.ann <- merge(ann.650Y,plot.ann,by.x="Sample",by.y="Sample",all=TRUE,sort=FALSE)
dim(sample.ann)

rownames(sample.ann) <- sample.ann[,"Sample"]


##### both the same size, put in order
dim(sample.ann)
dim(vec)

posns <- match(rownames(sample.ann),rownames(vec))
missing <- is.na(posns)
if(sum(missing)>0){print("ERROR missing annotations")}
vec <- vec[posns[!missing],]
sample.ann <- sample.ann[!missing,]

dim(vec)


sample.ann[1:5,]
dim(sample.ann)

#################################################
sample.ann[is.na(sample.ann[,color.with]),color.with] <- as.character(sample.ann[is.na(sample.ann[,color.with]),"status"])
sample.ann[is.na(sample.ann[,color.with]),color.with] <- "Data"

sort(tapply(sample.ann[,"Ethnicity"],sample.ann[,"Ethnicity"],length))
unique(sample.ann[,"Ethnicity"])  ### Here -9 indicates cases. Their ethnicities need to be inferred. 

# [1] "Adygei"        "Balochi"       "Bantu"         "Basque"        "Bedouin"       "Biaka_Pygmies" "Brahui"       
# [8] "Burusho"       "Cambodian"     "Colombian"     "Dai"           "Daur"          "Druze"         "French"       
# [15] "Han"           "Hazara"        "Hezhen"        "Italian"       "Japanese"      "Kalash"        "Karitiana"    
# [22] "Lahu"          "Makrani"       "Mandenka"      "Maya"          "Mbuti_Pygmies" "Melanesian"    "Miaozu"       
# [29] "Mongola"       "Mozabite"      "Naxi"          "Orcadian"      "Oroqen"        "Palestinian"   "Papuan"       
# [36] "Pathan"        "Pima"          "Russian"       "San"           "Sardinian"     "She"           "Sindhi"       
# [43] "Surui"         "Tu"            "Tujia"         "Tuscan"        "Uygur"         "Xibo"          "Yakut"        
# [50] "Yizu"          "Yoruba"        "-9"           


tapply(rownames(sample.ann),sample.ann[,color.with],length) # counts for individual from different ethnicities
unique(ann.650Y[,color.with])


############################## Get orginal colors for 650Y ###################
classes<-levels(as.factor(ann.650Y[,color.with]))  #order ok with the above
color.set<-rainbow(length(classes))          #get auto colors
color.set[10]<-"yellow4"
color.set[20]<-"forestgreen"
color.set[17]<-"lightblue"

pch.set<-c(0:25,33:45,47:58) # 650Y ethnicity
pch.set[21] <- pch.set[5]
pch.set[35] <- pch.set[2] 
pch.set[36] <- pch.set[5]
pch.set[28] <- pch.set[4]
pch.set[31] <- pch.set[3]
pch.set[34] <- pch.set[6]
pch.set[32] <- pch.set[16]
pch.set[33] <- pch.set[1]
pch.set[30] <- pch.set[9]
pch.set[38:51] <- pch.set[1:14]

classes.full <- levels(as.factor(sample.ann[,color.with]))
classes.extra <- classes.full[!(classes.full %in% classes)]

classes <- c(classes.extra,classes)

gray.steps <- seq(from=0,to=0.75,by=(1-0.25)/length(classes.extra))
color.extra <- gray(gray.steps[1:length(classes.extra)])
pch.extra <- rep(20,times=length(classes.extra))

color.set <- c(color.extra,color.set)
pch.set <- c(pch.extra,pch.set) # 650Y ethnicity (HapMap sets)

names(color.set) <- classes
names(pch.set) <- classes

color.array <- color.set[sample.ann[,color.with]]
pch.array <- pch.set[sample.ann[,color.with]]
length(color.array)
length(pch.array)

range.eig1 <- range(as.numeric(vec[,eig1]))
range.eig2 <- range(as.numeric(vec[,eig2]))

include <- c(1:dim(vec)[1])

#save.image("plot_test_ethnicity_650_strat_before_par_line.RData") ##save in case you need the raw data file

par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.1,1,0),mar=c(5,5,1,0)+0.1)
bp <- plot(as.numeric(vec[include,eig1]),as.numeric(vec[include,eig2]),pch=pch.array, cex=1,lwd=1.0,col = color.array[include], xlim=range.eig1,ylim=range.eig2,main="",xlab="Eigenvector1", ylab="Eigenvector1",cex.lab=1,cex.axis=1)
legend("topleft",legend= classes,col= color.set,pch=pch.set,bty="n",ncol = 3, cex=0.7, xjust= 0.5, x.intersp= 2.0, y.intersp=0.6,pt.cex=0.7)
##Save plots if needed
tiff(paste0(used.pca.tool,"PCA_PLOT.tiff"), width = 600, height = 600)
bp <- plot(as.numeric(vec[include,eig1]),as.numeric(vec[include,eig2]),pch=pch.array, cex=1,lwd=1.0,col = color.array[include], xlim=range.eig1,ylim=range.eig2,main="",xlab="Eigenvector1", ylab="Eigenvector1",cex.lab=1,cex.axis=1)
legend("topleft",legend= classes,col= color.set,pch=pch.set,bty="n",ncol = 3, cex=1.3, xjust= 0.5, x.intersp= 2.0, y.intersp= 0.7,pt.cex=0.7)
dev.off()
######################################## Using 4 Principal components (PCs), we can infer relationship between the individuals based on multidimensional euclidean distance

myref <- sample.ann
myref[,"Ethnicity"][myref[,"Ethnicity"]=="Case"|myref[,"Ethnicity"]=="Control"| myref[,"Ethnicity"]=="-9"] <- "unknown"
uIX = which(myref[,"Ethnicity"] == "unknown")

mydf <- vec
PCs <- c("e.1", "e.2", "e.3", "e.4") #number of principal components
mydf <- mydf[,PCs]  #selecting the number of PCs for calculating euclidean distance between each sample pair

#Calculating euclidean distance between each sample pair
dMat <- as.matrix(dist(mydf, method = "euclidean", diag = FALSE, upper = FALSE, p = 2))

dMat[uIX, uIX] <- Inf
nn <- apply(dMat, 1, order)[2, ]
myref[,"Ethnicity"][uIX] <- myref[,"Ethnicity"][nn[uIX]]
rownames(myref) <- myref[,"Sample"] # "myref" has inferred ethnicity for cases
#######End of euclidean distance calculation 

#######Now transfer the ethnicity to the IBD output file obtained after running IBD code. 
##If samples have leading 0's (sometimes samples names 890 are labelled as 00890 which need to be corrected), this will remove any leading 0's in sample name
are.related[,"sample_A"] <- sub("^0+([0-9]+)$", "\\1", are.related[,"sample_A"])
are.related[,"sample_B"] <- sub("^0+([0-9]+)$", "\\1", are.related[,"sample_B"])


posns <- match(are.related[,"sample_A"],myref[,"Sample"])
posns2 <- match(are.related[,"sample_B"],myref[,"Sample"])

matched.ethnicity <- cbind(are.related,Ethnicity_Sample_A=myref[posns,"Ethnicity"], Ethnicity_Sample_B=myref[posns2,"Ethnicity"])
#Writing the the file with specific ethnicity (non-generalized ethnicity straight from HapMap data) information and IBS values (IBS>= 0.06) for cases.
write.table(matched.ethnicity,file=paste0(path.save,"/","Non-generalized_inferred_ethnicity_PCA_with_",used.pca.tool,"_", length(PCs),"D_PCA.txt_RELATED"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

####################### Now, read generalized ethnicity ("Generalized_Ethnicity_definition.csv") file provided. 

generalized.ethnicity <- read.delim(paste0(ann.dir,"/Generalized_Ethnicity_definition.csv"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
posns <- match(matched.ethnicity[,"Ethnicity_Sample_A"],generalized.ethnicity[,"Ethnicity"])
sum(is.na(posns)) #has to be zero
matched.ethnicity <- cbind(matched.ethnicity, generalized.ethnicity_A=generalized.ethnicity [posns,"Generalised"])
posns2 <- match(matched.ethnicity[,"Ethnicity_Sample_B"],generalized.ethnicity[,"Ethnicity"])
sum(is.na(posns2)) #has to be zero
matched.ethnicity <- cbind(matched.ethnicity, generalized.ethnicity_B=generalized.ethnicity [posns2,"Generalised"])
#Writing the the very file with ethnicity information and IBS values (IBS>= 0.06) used in QC analysis. 
write.table(matched.ethnicity, file=paste0(path.save,"/","Generalized_ethnicity_","inferred_ethnicity_PCA_with_",used.pca.tool,"_", length(PCs),"D_PCA_RELATED.txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

####################################################################################END!!
####################################################################################END!!
####################################################################################END!!
####################################################################################END!!
####################################################################################END!!
####################################################################################END!!
####################################################################################END!!