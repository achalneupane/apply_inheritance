#This code should be used to process VCF files to calculate IBS, first using 1000 unfiltered SNPs, then doubling the number of SNPs until all SNPs are used. The final output file (i.e. IBS calculated using maximum number of SNPs; the total number of filtered and unfiltered SNPs used is indicated in the output file names ) should be used for the IBS relationship analysis. 

##Need the following packages
#source("https://bioconductor.org/biocLite.R")
#biocLite()
myPackages <- c("biomaRt", "dplyr", "stringi", "GenomicFeatures", "Rsamtools", "foreach", "doMC", "HardyWeinberg")

#install missing packages
my.installed.packages <- installed.packages()
available.packages <- myPackages %in% my.installed.packages
missing.packages <- myPackages[!available.packages]

if (length(missing.packages) > 0){
  print (paste0("Installing missing packages: ",missing.packages))  
  source("https://bioconductor.org/biocLite.R")
  biocLite((missing.packages), suppressUpdates=TRUE)
  }else {
  print ("ALL packages required for File1_ApplyinheritanceQC.r are installed!")
}

#force load all the packages required; try until all packages are loaded successfully! This will avoid any issues while working with multithreading.
tryCount <- 0    

while( !all(myPackages %in% (.packages())) ){
  
  try(require(biomaRt))
  try(require(dplyr))
  try(require(stringi))
  try(require(GenomicFeatures))
  try(require(Rsamtools))
  try(require(foreach))
  try(require(doMC))
  try(require(HardyWeinberg))
  tryCount <- tryCount + 1
  
  if( !all(myPackages %in% (.packages()))  ){
    cat(paste0("Failure: ", tryCount, "\n"))
    cat("Failed to load: ")
    cat(myPackages[ !myPackages %in% (.packages()) ])
    cat("\n")
  } else {
    print(paste0("Package loading successful!"))
  }
  
  Sys.sleep(5)
  
}

#path.save is where all the output files are saved. Change as needed***
path.save <- "/home/an/Desktop/qc_paper/subset_vcf/southasian/SNPs/IBD_output"

## Directory where subroutines 'annotate_SNPs_subroutines.r' and 'hwe.r' are located
code.dir <- "/home/an/Desktop/qc_paper/qcpaper_dtata"

#Project.dir is where subdirectory named 'SNPs' is present (i.e. /home/an/Desktop/qc_paper/subset_vcf/southasian/SNPs). Within /SNPs, all vcf files to be used are present. Change as needed***
project.dir <- "/home/an/Desktop/qc_paper/subset_vcf/southasian"

project.name <-  "AOGC-Genotyping.output" #This is the name that should be used before ".chr1." named for each vcf file used in the analysis. Here, example vcf files named:
#AOGC-Genotyping.output.chr1.recalibrated.filtered.vcf
#AOGC-Genotyping.output.chr2.recalibrated.filtered.vcf
#AOGC-Genotyping.output.chrX.recalibrated.filtered.vcf


the.sample.sheet <- "no_file"  # We will create samplesheet below when no sample sheet is available

###To avoid this error: Error in nchar(as.character(match.string)) : could not find symbol "keepNA" in environment of the generic function
setMethod("nchar", "ANY", base::nchar)
lengths <- function(x) vapply(x,length,1L)

###############################################################

skip.annovar.run <- FALSE # if FALSE it just redoes the files annovar summarization IF not already run
update.annovar.annotations <- TRUE # set TRUE will re-read VCF file and force annovar to always run
GATK.SB <- TRUE
genome.build <- "hg19"
dbSNP.build <- "131"
vcf.type="v4" # "annovar" "v4" "v3" "plink_assoc"

bam.extension <- ".ReCal.sort.bam"
combined.extension <- ""
snp.extension <- ".recalibrated.filtered.vcf"  # "_snps.raw.vcf"
indel.extension <- ".indelFiltered.vcf"    # "_DINDEL.raw.vcf"  
small.extension <- "_All_VARIANTS.raw.vcf"  # "_All_VARIANTS.raw.vcf"
variant.types <- c("snp") ##MUST have a extension type for each indel defined
names(variant.types) <- c("v4") ### define data type for when reading below
use.variants <- c("snp")
################

##Adjust the required number of cores ('num.cores')
num.cores <- 8 
num.cores.ori <- num.cores
registerDoMC(cores = num.cores)
if(!exists("read.in.all.vcf")){read.in.all.vcf <- FALSE} # if FALSE will use Max.reads to read in file in chunks
if(!exists("max.reads")){max.reads <- 1000 } # 0- get ALL IN ONE GO  :
## NA get all
## 415 sample 700000 (7 cores) < 20 Gb

###############10000 -1000 samples about 23Gb RAM' 12hrs
###############ALL   -58 samples about  23Gb 30 mins

if(!exists("project.name.alternative")){project.name.alternative <- project.name}

## Removed by Mhairi 9th Feb 2016 as should be set in params at top of script
###project.dir <- paste(UQCCG.data,project,sep="/")

snp.dir <- paste(project.dir,"SNPs",sep="/")
bam.dir <- paste(project.dir,"BAM",sep="/") # bam directory where it contains the sample sheet else disregard this. 
small.dir <- paste(project.dir,"SNPs",sep="/")
indel.dir <- paste(project.dir,"DINDELs",sep="/")
analysis.dir <- paste(project.dir,"Analysis",sep="/")

###Approximate time for processing:
###############10000 -1000 samples about 23Gb RAM' 12hrs
###############ALL   -58 samples about  23Gb 30 mins

############################################################ Now set up the filter options; leave as is

if(read.in.all.vcf){max.reads <- 0} # max.reads of 0 will cause it to read in the whole file
missing.threshold <- 0.80 # less.than  0.80  #  missing.threshold <- 1.0
inbreeding.threshold <- 0.75 # less than must be (1-inbreeding.threshold)*100 PERCENT NOT 0/1 OR MISSING # inbreeding.threshold <- 1
hwe.threshold <- 1e-4
GQ.threshold <- 20
use.variants <- c("snp") ### THIS MUST BE LOWER CASE SNP ELSE GENOTYPE FILTERING WILL FAIL

QC.dir <- paste(project.dir,"QC",sep="/")
if(!exists("annotate.dir")){annotate.dir <- paste(project.dir,"Annotate",sep="/")}
if(!exists("analysis.dir")){analysis.dir <- annotate.dir}

############## Accessing code files and subroutines
core.ann <- c("chr","start","end","REF","ALT","TYPE") 
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")
################# Setting up quality filtering; leave as is
global.quality.labs <- c("QUAL","QD","FS","HRun","SB","FILTER","FILTER","TYPE","STR") ### these become the "good.qual" filter
global.quality.names <- c("QUAL","QD","FS","HRun","SB","FILTER_PASS","FILTER_100","flat","notREPEAT")
global.quality.cut <- c(50,0.5,60,5,1,"PASS","99.90to100.00","flat","0")
global.quality.type <- c("numeric","numeric","numeric","numeric","numeric","factor","factor","factor","factor")
global.quality.dirn <- c("greater","greater","less","less","less","exact","ends_with","not_ends_with","exact")

pass.filter <- c("FILTER_PASS","QUAL","QD","FS","HRun","SB","flat","notREPEAT")  ## all these are required to be TRUE  (can have FILTER AND FILTER_100


names(global.quality.cut) <- global.quality.labs
names(global.quality.dirn) <- global.quality.labs
names(global.quality.type) <- global.quality.names

global.quality.cut
global.quality.dirn
global.quality.type
####################################################################################END: set up the filter options

variant.types <- variant.types[variant.types %in% use.variants]

i <-  1
print(variant.types)
for(i in 1:length(variant.types)){  ### CHOOSE BY project.name * extension
  files <- dir( eval(as.name( paste(variant.types[i],"dir",sep=".") )) )
  print((variant.types[i]))
  print(files)
  the.extension <- paste(eval(as.name( paste(variant.types[i],"extension",sep=".") )),"$",sep="")
  the.extension <- paste(combined.extension,the.extension,sep="")
  the.files <- files[grepl(the.extension ,files)]
  the.files <- the.files[grepl(paste("^",project.name,sep=""),the.files)]
  
  names(the.files) <- gsub(the.extension,"",the.files)
  assign( paste(variant.types[i],"files",sep="."),value=the.files)
}

snp.files
snp.dir

######

## To check if the snp.files are ok when running job in HPC cluster. This will print out the list of all vcf files ("Used_snpfiles.txt") to be analysed.
if(length(snp.files)>0){
  check.snp.files <- as.data.frame(snp.files)
  write.table(check.snp.files,file= paste(path.save,"Used_snpfiles.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}

########################### QC analysis only use SNPs
sample.files <- snp.files
names(sample.files) <- rep("snp",times=length(sample.files))

a.single.test <- FALSE ## need to be FALSE to overwrite Aij through each looping process
# isamp <- 1 ; a.single.test <- TRUE

######################################################
#Note: This is where the loop begins; the purpose of this loop is to get IBD output starting from first 1000 input SNPs, using the doubled number of input SNPs, then ultimaltely using all SNPs in the VCF files. IBD values using maximum number of SNPs should be considered for the analysis. The output file will indicate total number of SNPs (both filtered and unfiltered) used
#############################

total.nos.snps.indel.unfiltered <- {}
total.nos.filtered.genotypes <- {}
doubled.SNPs <- {}
isamp <- 1
for (isamp in 1:length(sample.files)){

#   #Stop switch to check the generated data
#   if(isamp==2){
#     print ("Stopping to check the generated data")
#     stop 
#   }
  
print(paste("Doing: ",sample.files[isamp],sep=""))

################## get initial read information
start.data <- prepare.for.Vcf.file.read(sample.files[isamp])
for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}
###########################################



setwd(eval(as.name(paste(names(sample.files)[isamp],"dir",sep="."))))
con <- file(sample.files[isamp], open="r")  # close(con)
num.lines <- 1 # so does while llop at least once
reads <- max.reads  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter <-  -1
while (num.lines >0){
counter <- counter+1
print(paste("Reading chunk, counter number",counter,sep=":"))

if(counter==0){
indels <- try(scan(con,what=character(num.vars),skip=(reads*counter)+skip.lines,nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}else{
indels <- try(scan(con,what=character(num.vars),nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}

num.lines <- length(indels)/(num.vars)
print(num.lines)
if(num.lines==0){next}
dim(indels) <- c(num.vars,num.lines)
indels <- t(indels)
colnames(indels) <- column.labels

#} ###Test for indel

if(dim(indels)[1]<100){next} ### not worth the effort... not enough data

nos.snps.indel <- dim(indels)[1]
total.nos.snps.indel.unfiltered <- sum(total.nos.snps.indel.unfiltered,nos.snps.indel)
print(paste0("Total raw/unfiltered SNPs in indel:",total.nos.snps.indel.unfiltered))

format.posn <- match("FORMAT",colnames(indels))
samples.processing <- length(column.labels)-format.posn # number of samples in the vcf file
samples.order.in.ALL <- column.labels[(format.posn+1):length(column.labels)] # samlple labels
the.samples <- samples.order.in.ALL

#assigning filter for min allele frequency 
if(length(the.samples)>=100){
allele.min <- max(c(10,as.integer(samples.processing/40))) 
}
if(length(the.samples)>40 & length(the.samples)< 100 ){
allele.min <- max(c(8,as.integer(samples.processing/10)))
}
if(length(the.samples)<=40 & length(the.samples)>=15  ){
allele.min <- 2
}
if(length(the.samples)<15  ){
allele.min <- 1
}
                                        # typically 10 
#allele.min/length(the.samples) # 10/40
print(paste0("allele.min:",allele.min))
####################################### FINISHED Read in data DO PROCESSIng below
###################################################################################################

indels <-   suppressWarnings(process.GATK.indels(indels,sample.files[isamp],vcf.type,format.types,info.types,info.class,num.cores.ori)) ## names(the.sample.file) must be the mutation type

rownames(indels) <- build.key(indels,core.ann,add.chr.label=FALSE)
######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
####################################################################################################

all.sample.labels <- colnames(indels)
all.sample.labels <- all.sample.labels[grep("\\.GT$",all.sample.labels)]
all.sample.labels <- gsub(".GT$","",all.sample.labels)
all.sample.labels[1:10]
############################################33

if(!grepl("^chr",indels[1,"chr"])){
key.indels <- build.key(indels,core.ann,add.chr.label=TRUE)
indels[,"chr"] <- paste("chr",indels[,"chr"],sep="") ## REQUIRE "chr in chromosome labels below
rownames(indels) <- key.indels
}else{key.indels <- build.key(indels,core.ann)      }

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

the.samples <- all.sample.labels
the.samples
print("Get filtered genoypes QC")
## all.genotypes <- filtered.genotype.summary(indels,the.samples,prefix="",suffix=".ALL",20,0.2,0.80,0.05,0.95,10,5,GQ.threshold) # most often used

num.bits <- num.cores.ori ### can make num.bits larger and keep num cores small to reduce memeory
while((dim(indels)[1] %% num.bits)< 2){num.bits <- num.bits+1} ### go don't get matrix issues
print(paste0("Filtered.genotypes num.bits: ",num.bits))
all.genotypes <- foreach(indels.bit=iter(indels,by='row',chunksize=as.integer(dim(indels)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(indels.bit,the.samples,prefix="",suffix=".ALL",20,0.2,0.80,0.05,0.95,10,5,GQ.threshold)


dim(all.genotypes)
targets <- c("ALL") #c("NMD","ex.Control","AOGC")
names(targets) <- targets
targets
#use.samples <- the.samples[pheno[,"SampleProject"]==targets[it]]
summary.geno.group <- genotype.summary(as.matrix(all.genotypes[,paste0(the.samples,".GT")]))
colnames(summary.geno.group) <- paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets),sep=".")


position.filter.full <- position.quality.filter(indels,global.quality.cut,global.quality.type,global.quality.dirn)
pass.filter.have <- pass.filter[pass.filter %in% colnames(position.filter.full)]

## position.filter.full[1:5,]

## names(all.genotypes)
dim(all.genotypes)
#referencing the paper (Joseph Powell et al. [november 2010], "Reconciling the analysis of IBD and IBS in complex trait studies"), we can now use filtered genotypes to calculate Aij
#http://www.nature.com/nrg/journal/v11/n11/extref/nrg2865-s1.pdf

print("start QC")
all.genotypes[all.genotypes=="NA"] <- NA
all.genotypes[all.genotypes=="0/0"] <- 0
all.genotypes[all.genotypes=="0/1"] <- 1
all.genotypes[all.genotypes=="1/1"] <- 2

p <- as.numeric(summary.geno.group[,"MAF.ALL"])
p[is.na(p)] <- 0

hw.p.ok <- getHWE(summary.geno.group[,"GENO.ALL"])
hw.p.ok <- hw.p.ok >  hwe.threshold


#sum(hw.p.ok)
#test[1:500]

## Hg18 pseudoautosomal regions (PARs)
## chrY:1-2709520 and chrY:57443438-57772954 
## chrX:1-2709520 and chrX:154584238-154913754


## Hg19 pseudoautosomal regions (PARs)
## chrY:10001-2649520 and chrY:59034050-59363566 
## chrX:60001-2699520 and chrX:154931044-155260560

#mouse PAR: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC311143/figure/F6/  Fxy to tel:  166350000  End: 166,404,117
if(genome.build=="mm9"){
  parX <- (indels[,"chr"]=="chrX" & (indels[,"start"]> 166350000 & indels[,"start"]< 166650000) )
  print("Pseudo-Autosomal region on mouse mm9 coordinates")
}else{
parX <- (indels[,"chr"]=="chrX" & (indels[,"start"]> 60001 & indels[,"start"]< 2699520) ) | (indels[,"chr"]=="chrX" & (indels[,"start"]> 154931044 & indels[,"start"]< 155260560) )
print("Pseudo-Autosomal region on human hg19 coordinates")
}

on.X <- indels[,"chr"]=="chrX"

print(paste0("on.x; ",sum(on.X)))
summary.geno.group[1:5,]

start.duplicated <- duplicated(indels[,"start"]) # don't need cause used TYPE==snp==exact
is.unwound.geno <- grepl("snp:\\d+$",indels[,"TYPE"]) | grepl("indel:\\d+$",indels[,"TYPE"])

missing.rate <- (as.numeric(summary.geno.group[,"MISSING.Alleles.ALL"])/2)/length(the.samples)
missing.ok <- missing.rate < missing.threshold  ## genotype is call not NA

inbreeding.rate <- ((as.numeric(summary.geno.group[,"MISSING.Alleles.ALL"])/2) + as.numeric(summary.geno.group[,"ALT_HETRO.ALL"]))/length(the.samples)
inbreeding.ok <-  (inbreeding.rate < inbreeding.threshold) | on.X

REF.length <- nchar(as.character(indels[,"REF"]))
ALT.length <- nchar(as.character(indels[,"ALT"]))

large.indel <- REF.length>1 | ALT.length>1

are.repeats <- identify.repeats(indels,di.run.max=2,homo.run.max=4)

length(large.indel) 
length(are.repeats)
sum(are.repeats)
#rownames(indels)[are.repeats][1:20]

#################### Now searching repeats by looking  forward

chk.in.repeat <- large.indel & !are.repeats
if(sum(chk.in.repeat)>1){
are.sub.repeat <- indentify.IN.repeat(indels[chk.in.repeat,],looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats <- key.indels[chk.in.repeat][are.sub.repeat]
are.in.repeats.forward <- key.indels %in% remove.repeats
}else{
  are.in.repeats.forward <- rep(FALSE,times=length(REF.length))
}
  
#remove.repeats[1:20]
sum(are.in.repeats.forward)


###################### Now searching repeats by looking backward

sum(chk.in.repeat)
chk.in.repeat <- large.indel & !are.repeats & !are.in.repeats.forward
if(sum(chk.in.repeat)>1){
are.sub.repeat <- indentify.IN.repeat(indels[chk.in.repeat,],looking="back",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats <- key.indels[chk.in.repeat][are.sub.repeat]
are.in.repeats.back <-  key.indels %in% remove.repeats
}else{
  are.in.repeats.back <- rep(FALSE,times=length(REF.length))
}
  
sum(are.in.repeats.back)

are.in.repeats <-  are.in.repeats.back | are.in.repeats.forward

length(are.in.repeats)
sum(are.in.repeats)

##########################################################################
##########################################################################

maf.fil <- (p >= allele.min/length(the.samples)) & (p <= (1-allele.min/length(the.samples))) # maf.fil <- (p >= 0.01 & p <= 0.1)
maf.fil[is.na(maf.fil)] <-  FALSE

all.true <- length(pass.filter.have) ### only use snps and higly filtered data in the QC: position.filter.full[1:5,pass.filter]
position.filter.QC <- apply(position.filter.full[,pass.filter.have],1,function(x) sum(x)==all.true) # sum(position.filter)

position.filter.QC <- position.filter.QC & maf.fil &  missing.ok  & !is.unwound.geno & !are.in.repeats & !are.repeats & (hw.p.ok | on.X) & !start.duplicated & !is.unwound.geno  # hw and x chromosome don't mix!



sum(maf.fil)
sum(missing.ok)
sum(inbreeding.ok)
sum(!is.unwound.geno)
sum(!are.in.repeats)
sum(!are.repeats)
sum(hw.p.ok)
print(sum(position.filter.QC))

############################################# Now working on chunks to run the process fast

indels <- indels[position.filter.QC,]
summary.geno.group <- summary.geno.group[position.filter.QC,]
all.genotypes <- all.genotypes[position.filter.QC,]

fil.nos.genotypes <- dim(all.genotypes)[1]
total.nos.filtered.genotypes <- sum(total.nos.filtered.genotypes,fil.nos.genotypes)
print(paste0("Total number of filtered genotype/SNPs from all.genotypes:",total.nos.filtered.genotypes))

p <- p[position.filter.QC]
on.X <- on.X[position.filter.QC]
parX <- parX[position.filter.QC]
position.filter.QC <- position.filter.QC[position.filter.QC] # position.filter.QC now lost

############################################################ subset to make fast

print(paste("The dimension of indel is",dim(indels)[1],dim(indels)[2],sep=":"))
print(paste("The working directory is",getwd(),sep=":"))
## getwd()
# save.indel <- rbind(save.indel,indels)
#print(paste0("the dimension of save.indel:",dim(save.indel)[1],":",dim(save.indel)[2]))

#############################################################

if(sum(position.filter.QC)>1){
num.bits <- num.cores.ori
while((dim(all.genotypes)[1] %% num.bits)< 2){num.bits <- num.bits+1} ### to avoid getting matrix issues
num.bits
if(dim(all.genotypes)[1]<200){num.bits <- 1; print(paste("Only",dim(all.genotypes)[1],"Filtered genotypes present",sep=" ")) }
  print(paste0("do Aij calculation num.bit:",num.bits))
Aij <- foreach(genotypes.bit=iter(as.matrix(c(1:(dim(all.genotypes)[1]))),by='row',chunksize=as.integer(dim(all.genotypes)[1]/num.bits) ), .combine='+', .multicombine=FALSE,.inorder=FALSE) %dopar% genetic.sample.QC.accumilate(all.genotypes[as.integer(genotypes.bit),],p[as.integer(genotypes.bit)],position.filter.QC[as.integer(genotypes.bit)],the.samples,on.X[as.integer(genotypes.bit)],parX[as.integer(genotypes.bit)])

}else{
  print("WARNING  0 position left after filtering skipping Aij calculation")
}
print(paste0("num.cores and num.bits",num.cores," : ",num.bits))
if((counter==0 & isamp==1) | a.single.test ){
Aij.total <- Aij
a.single.test <- FALSE
print("first assign")
plink.files.all <- paste(gsub(".vcf$","",sample.files[isamp]),"plink",counter,sep="_")
}else{
Aij.total <- Aij.total+Aij
plink.files.all <- c(plink.files.all,paste(gsub(".vcf$","",sample.files[isamp]),"plink",counter,sep="_"))
print("Add Aij")
}



#####################

if((isamp==1) & (counter==0)& total.nos.snps.indel.unfiltered>=500){
 
  related <- sum.QC.matrix(Aij.total,the.samples) 
  mat.Aij.total <- as.matrix(Aij.total)
  snp.columns <-  colnames(mat.Aij.total)[grepl("_C$",colnames(mat.Aij.total))]
  mat.Aij.total  <- mat.Aij.total[,colnames(mat.Aij.total)%in% snp.columns]
  mat.Aij.total  <- Aij.total[rownames(mat.Aij.total)%in% snp.columns,]
  i1 <- lower.tri(mat.Aij.total, diag=TRUE)
  i2 <- which(i1, arr.ind=TRUE)
  df.Aij  <- data.frame(sampleA = colnames(mat.Aij.total)[i2[,1]], 
                      sampleB = colnames(mat.Aij.total)[i2[,2]], SNPs = mat.Aij.total[i1])
  df.Aij[,"sampleA"] <-  sub("^0+([0-9]+)$", "\\1", df.Aij[,"sampleA"])
  df.Aij[,"sampleB"] <-  sub("^0+([0-9]+)$", "\\1", df.Aij[,"sampleB"])
  #df.Aij  <- cbind(key.Aij= paste(gsub("_C","",df.Aij[,"sampleA"]),gsub("_C","",df.Aij[,"sampleB"]),sep=":"),df.Aij)
  df.Aij  <- cbind(key.Aij= paste(gsub("_C","",df.Aij[,"sampleB"]),gsub("_C","",df.Aij[,"sampleA"]),sep=":"),df.Aij) #May 10 changed
  df.Aij <- as.matrix(df.Aij)
 
  ###if the the.samples (i.e. the column names in the indel has leading 0 [sometimes sample names such as "890" are labelled as 00890]), this will fix it. Sample names must match exactly!
  the.samples <- sub("^0+([0-9]+)$", "\\1", the.samples)
  related[,"sample_A"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_A"])
  related[,"sample_B"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_B"])
  
  related <- cbind(key.related= paste(related[,"sample_A"],related [,"sample_B"],sep=":"),related)
  related <- merge(related,df.Aij,by.x="key.related",by.y= "key.Aij")
  
  
  ###########################################################
  
  the.samples <- sub("^0+([0-9]+)$", "\\1", the.samples)
  related[,"sample_A"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_A"])
  related[,"sample_B"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_B"])

  ################################Now creating sample sheet with same sample names in vcf file
  
  if(exists("the.sample.sheet") & the.sample.sheet!="no_file"){
    
    xx <- try(sample.sheet.full <- read.delim(paste(bam.dir,the.sample.sheet,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE))
    if(inherits(xx, "try-error")){
      sample.sheet.full <- read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
    }else{  ## no sample sheet exists
      #sample.sheet.full <- unique(related[,"sample_A"])
      sample.sheet.full <- the.samples
      
      sample.sheet.full <- cbind(sample.sheet.full,sample.sheet.full,sample.sheet.full)
      colnames(sample.sheet.full) <- c("SampleProject","ParticipantCode","Sex")
      sample.sheet.full[,"SampleProject"] <- project.name
      sample.sheet.full[,"Sex"] <- 9
    }
    if(!("ParticipantCode" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
      if(sum(grepl("^Participant",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Participant",colnames(sample.sheet.full) )] <- "ParticipantCode"}
    }
    
    if(!("SampleProject" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
      if(sum(grepl("^Sample Project",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Sample Project",colnames(sample.sheet.full) )] <- "SampleProject"}
    }
    
    if(!("SampleProject" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
      if(sum(grepl("^Sample.Project",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Sample.Project",colnames(sample.sheet.full) )] <- "SampleProject"}
    }
  } else{   ## to create a sample.sheet.full object if  sample sheet exists
    sample.sheet.full <- unique(related[,"sample_A"])
    sample.sheet.full <- cbind(sample.sheet.full,sample.sheet.full,sample.sheet.full)
    colnames(sample.sheet.full) <- c("SampleProject","ParticipantCode","Sex")
    sample.sheet.full[,"SampleProject"] <- project.name
    sample.sheet.full[,"Sex"] <- 9
  } #exists("the.sample.sheet"
  
   #########################################################################333
  samples.not.in.samplesheet <-  the.samples[!(the.samples %in% sample.sheet.full[,"ParticipantCode"])]
  if(length(samples.not.in.samplesheet)>0){print(paste("WARNING Samples missing from SampleSheet: ",toString(samples.not.in.samplesheet),sep=""))}
  posns <- match(related[,"sample_A"],sample.sheet.full[,"ParticipantCode"])
  diagonal <- related[,"sample_A"]==related[,"sample_B"]
  Ident.Samples <- diagonal
  related.sheet <- cbind(related,Ident.Samples,sample.sheet.full[posns,])
  #added later 
  colnames(related.sheet)[grepl("IBS",colnames(related.sheet))] <- paste0("IBS_from:",total.nos.snps.indel.unfiltered,"SNPS")
  colnames(related.sheet)[grepl("Num_Good_SNPs_A",colnames(related.sheet))] <- paste0("Num_Good_SNPs_A:",total.nos.snps.indel.unfiltered,"SNPS") 
  ########

  setwd(path.save)
  out.file.name <- paste0(path.save,paste(total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered_","_SNPs_","Vcf_Merge.","GQ-",GQ.threshold,".Min.ALT-",allele.min,"_",project.name.alternative,".",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep=""))
  write.table(related.sheet,out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
  all.related.sheet.ones  <- related.sheet
  all.related.sheet.ones <-  cbind(key.related.sheet = paste(all.related.sheet.ones[,"sample_A"],all.related.sheet.ones[,"sample_B"],sep=":"), all.related.sheet.ones)        
  
  diagonal <- related.sheet[,"sample_A"]==related.sheet[,"sample_B"]
  bad.sex <- (related.sheet[,"sex_Predicted"]!=related.sheet[,"Sex"] | related.sheet[,"sex_Predicted"]==9  ) & diagonal & !is.na(related.sheet[,"Sex"])
  are.related <- as.numeric(as.character(related[,"IBS"]))>0.06 & is.finite(as.numeric(as.character(related[,"IBS"])))
  out.file.name <- paste0(path.save,paste(total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered_","_SNPs_","Vcf_Merge.","GQ-",GQ.threshold,".Min.ALT-",allele.min,"_",project.name.alternative,".ARE_RELATED.",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep=""))
  write.table(related.sheet[are.related | bad.sex,],out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
  related.ones <- related.sheet[are.related | bad.sex,]
  related.ones <-  cbind(key.related.ones= paste(related.ones[,"sample_A"], related.ones[,"sample_B"],sep=":"), related.ones)                            
  doubled.SNPs  <- total.nos.snps.indel.unfiltered*2 #doubling the number of SNPs
  save.image(paste0(path.save,total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered.Rdata"))

  }else if (total.nos.snps.indel.unfiltered >= doubled.SNPs){
    related <- sum.QC.matrix(Aij.total,the.samples) 
    mat.Aij.total <- as.matrix(Aij.total)
    snp.columns <-  colnames(mat.Aij.total)[grepl("_C$",colnames(mat.Aij.total))]
    mat.Aij.total  <- mat.Aij.total[,colnames(mat.Aij.total)%in% snp.columns]
    mat.Aij.total  <- Aij.total[rownames(mat.Aij.total)%in% snp.columns,]
    i1 <- lower.tri(mat.Aij.total, diag=TRUE)
    i2 <- which(i1, arr.ind=TRUE)
    df.Aij  <- data.frame(sampleA = colnames(mat.Aij.total)[i2[,1]], 
                      sampleB = colnames(mat.Aij.total)[i2[,2]], SNPs = mat.Aij.total[i1])
    df.Aij[,"sampleA"] <-  sub("^0+([0-9]+)$", "\\1", df.Aij[,"sampleA"])
    df.Aij[,"sampleB"] <-  sub("^0+([0-9]+)$", "\\1", df.Aij[,"sampleB"])
    df.Aij  <- cbind(key.Aij= paste(gsub("_C","",df.Aij[,"sampleB"]),gsub("_C","",df.Aij[,"sampleA"]),sep=":"),df.Aij)  #May 10 changed
    df.Aij <- as.matrix(df.Aij)
    related <- sum.QC.matrix(Aij.total,the.samples) 
    ###if the the.samples (i.e. the column names in the indel has leading 0 [sometimes sample names such as "890" are labelled as 00890]), this will fix it. Sample names must match exactly!
    the.samples <- sub("^0+([0-9]+)$", "\\1", the.samples)
    related[,"sample_A"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_A"])
    related[,"sample_B"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_B"])
    related <- cbind(key.related= paste(related[,"sample_A"],related [,"sample_B"],sep=":"),related)
    related <- merge(related,df.Aij,by.x="key.related",by.y= "key.Aij")
    the.samples <- sub("^0+([0-9]+)$", "\\1", the.samples)
    related[,"sample_A"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_A"])
    related[,"sample_B"] <-  sub("^0+([0-9]+)$", "\\1", related[,"sample_B"])
    posns <- match(related[,"sample_A"],sample.sheet.full[,"ParticipantCode"])
    diagonal <- related[,"sample_A"]==related[,"sample_B"]
    Ident.Samples <- diagonal
    related.sheet <- cbind(related,Ident.Samples,sample.sheet.full[posns,])
    related.sheet <-  as.matrix(related.sheet)
  
  ################################# Now saving output file
    out.file.name <- paste0 (path.save,paste(total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered_","_SNPs_","Vcf_Merge.","GQ-",GQ.threshold,".Min.ALT-",allele.min,"_",project.name.alternative,".",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep=""))
    write.table(related.sheet,out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
    unrelated.related.sheet <- cbind(key.related.sheet=paste(related.sheet[,"sample_A"],related.sheet[,"sample_B"],sep=":"),related.sheet)
    unrelated.related.sheet <- unrelated.related.sheet[,c("key.related.sheet","IBS","Num_Good_SNPs_A")]
    colnames(unrelated.related.sheet) <- c("key.related.sheet",paste0("IBS_from:",total.nos.snps.indel.unfiltered,"SNPS"),paste0("Num_Good_SNPs_A:",total.nos.snps.indel.unfiltered,"SNPS") )
    all.related.sheet.ones <-  merge(unrelated.related.sheet, all.related.sheet.ones, by="key.related.sheet", all = TRUE)
    out.file.name <- paste0(path.save,paste(total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered_","_Accumulated_SNPs_","Vcf_Merge.","GQ-",GQ.threshold,".Min.ALT-",allele.min,"_",project.name.alternative,".",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep=""))
#   #uncomment the line below to get accumulated IBD file starting from 1000 SNPs (both related and unrelated)
#   write.table(all.related.sheet.ones,paste0(out.file.name),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
    diagonal <- related.sheet[,"sample_A"]==related.sheet[,"sample_B"]
    bad.sex <- (related.sheet[,"sex_Predicted"]!=related.sheet[,"Sex"] | related.sheet[,"sex_Predicted"]==9  ) & diagonal & !is.na(related.sheet[,"Sex"])
    are.related <- as.numeric(as.character(related[,"IBS"]))>0.06 & is.finite(as.numeric(as.character(related[,"IBS"])))
    out.file.name <- paste0(path.save,paste(total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered_","_SNPs_","Vcf_Merge.","GQ-",GQ.threshold,".Min.ALT-",allele.min,"_",project.name.alternative,".ARE_RELATED.",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep=""))
    write.table(related.sheet[are.related | bad.sex,],out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
    new.related.ones  <- related.sheet[are.related | bad.sex,]
    new.related.ones <- cbind(key.related.ones=paste(new.related.ones[,"sample_A"],new.related.ones[,"sample_B"],sep=":"), new.related.ones)
    new.related.ones <- new.related.ones[,c("key.related.ones","IBS","Num_Good_SNPs_A")]
    colnames(new.related.ones) <- c("key.related.ones",paste0("IBS_from:",total.nos.snps.indel.unfiltered,"SNPS"), paste0("Num_Good_SNPs_A:",total.nos.snps.indel.unfiltered,"SNPS") )
    related.ones  <- merge(new.related.ones, related.ones, by="key.related.ones", all = TRUE)
    out.file.name <- paste0(path.save,paste(total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered_","_Accumulated_SNPs_","Vcf_Merge.","GQ-",GQ.threshold,".Min.ALT-",allele.min,"_",project.name.alternative,".ARE_RELATED.",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep=""))
#   #uncomment two lines below to get accumulated related IBD file starting from 1000 SNPs
#   write.table(related.ones,paste0(out.file.name),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#   save.image(paste0(path.save,total.nos.snps.indel.unfiltered,"_unfiltered_",total.nos.filtered.genotypes,"_filtered.Rdata"))
    doubled.SNPs  <- total.nos.snps.indel.unfiltered*2 #doubling the number of SNPs

  }

###################

  } ## while loop over data chunks


close(con)

} ## loop over isamp  -sample files which are all snps

#Note: This will generate two sets of accumulated and individual output files as shown below:
# 1. 'Accumulated_SNPs_Vcf_Merge.GQ-20.Min.ALT-2_AOGC-Genotyping.output.ARE_RELATED'
# 2. 'Accumulated_SNPs_Vcf_Merge.GQ-20.Min.ALT-2_AOGC-Genotyping.output'
# 3. 'SNPs_Vcf_Merge.GQ-20.Min.ALT-2_AOGC-Genotyping.output.ARE_RELATED'
# 4. 'SNPs_Vcf_Merge.GQ-20.Min.ALT-2_AOGC-Genotyping.output'
#File 1 has the IBD using all SNPs; and includes IBD >= 0.06
#File 2 has the IBD using all SNPs; and includes IBD = all values
#File 3 has the IBD using indicated number of SNPs; and includes IBD >= 0.06
#File 4 has the IBD using indicated number of SNPs; and includes IBD = all values
#File 'Accumulated_SNPs_Vcf_Merge.GQ-20.Min.ALT-2_AOGC-Genotyping.output.ARE_RELATED' with maximum number of filtered and unfiltered SNPs should be considered for IBD analysis. This is the output file that is generated after using all good quality SNPs.

######################### END!! ####################################################
######################### END!! ####################################################
######################### END!! ####################################################
######################### END!! ####################################################
######################### END!! ####################################################
######################### END!! ####################################################
######################### END!! ####################################################