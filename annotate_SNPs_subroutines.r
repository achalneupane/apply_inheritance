#This code will process VCF files to calculate IBS, first using 1000 unfiltered SNPs then doubling the number of SNPs until all SNPs are used up. The final output file (i.e. IBS calculated using maximum number of SNPs; the total number of filtered and unfiltered SNPs used is indicated in the output file names ) should be used for the IBS relationship analysis. 
#Processing GATK file structure
process.GATK.indels<-function(indels,the.sample.file,vcf.type,format.types,info.types,info.class,num.cores){
  #   print(dim(indels))
  column.labels <- colnames(indels)
  #######################Process the INFO section common to all indels
  ###
  #print(num.cores)
  ###
  if(vcf.type=="v3"){
    hetero<-indels[,dim(indels)[2]]=="0/1"
    indels[!hetero,"INFO"]<- paste("AB=1",indels[!hetero,"INFO"],sep=";")
  }
  
  ########### to process using lists  use this method timing:  52.91s    0.05  ( {   } )
  
  
  ## system.time(
  ##  {
  
  ####multicore way WORKS! safe but slower  74.410  16.490  21.155 compare to  11.110   0.010  11.124 
  ## info<-strsplit(indels[,"INFO"],split=";")   
  ## key.indel<-build.key(indels,c("chr","POS","REF","ALT"))
  ## names(info)<-key.indel
  ## info<-mclapply(info,process.info.list,mc.cores=num.cores)
  ## posns<-match(key.indel,names(info))
  ## info<-info[posns]
  ## num.info<-length(info.types)
  ## len.info<-length(info)
  ## info<-unlist(info)
  ## dim(info)<-c(num.info,len.info)  
  ## info<-t(info) # info[1:5,]
  ## colnames(info)<-info.types
  
  ########vector.method WORKS in # 11.110   0.010  11.124  indels[53,"INFO"]
  info<-indels[,"INFO"]
  #if(!skip.info.prosessing){  # skipping info processing causes error in fix multi alleles since that used AC anf AF
  
  flag.items <- info.types[info.class=="Flag"]
  for(iflag in 1:length(flag.items)){
    info<-gsub(paste(flag.items[iflag],";",sep=""),paste(flag.items[iflag],"=1;",sep=""),info) # is flag present just adds a one
    info<-gsub(paste(flag.items[iflag],"$",sep=""),paste(flag.items[iflag],"=1",sep=""),info) # if flag is the last element in teh string
  }
  
  
  info<-gsub("=",";",info) # get read to split
  info<-strsplit(info,split=";") # unique(unlist(lapply(info,length)))
  info.length<-unlist(lapply(info,length)) # table(info.length)
  
  #abs(c(1,1,2.1,-2.6,0.001)-trunc(c(1,1,2.1,-2.6,0.001))
  
  if(sum(info.length %% 2)!=0){print("Error in processing INFO unpaired values in process.GATK.indels")}
  
  
  
  info<-unlist(info)
  wanted<-seq(from=2, to=length(info),by=2)
  wanted.labels<-seq(from=1, to=length(info),by=2)
  
  info.labels<-info[wanted.labels]
  info.num<-info[wanted]
  info.length<-info.length/2
  info.index<-rep(1:length(indels[,"INFO"]),times=info.length)
  info<-matrix(data=0,nrow=length(info.length),ncol=length(info.types))
  colnames(info)<-info.types
  
  # i.info.col<-1
  for(i.info.col in 1:length(info.types)){
    # print(info.types[i.info.col])
    the.col<-info.types[i.info.col]
    posns.in.flatten.info<-grep(paste("^",the.col,"$",sep=""),info.labels,fixed=FALSE) ## can have AC or AC1 etc correct make 01/21/2013
    missing<-is.na(as.numeric(posns.in.flatten.info))
    sum(missing)
    posns<-info.index[posns.in.flatten.info]
    info[posns,the.col]<-info.num[as.numeric(posns.in.flatten.info)]
  }
  # numerical data wanted
  #} # skip info 
  
  ## }) # system time test
  posns
  sum(is.na(info.num[posns.in.flatten.info]))
  sum(info.num[posns.in.flatten.info]=="NA")
  test<-info.num[posns.in.flatten.info]
  
  chk<-is.na(as.numeric(test))
  info.num[posns.in.flatten.info][chk][1:10]
  grep("NA",as.character(chk))[1:10]
  sum(is.na(as.numeric(test)))
  sum(test=="NA")
  info[posns,the.col]<-test
  info[1:5,]
  ## to.test<-15
  ## a.na<-is.na(info.ori[,to.test]) # | (is.na(info[,to.test]) | is.na(info.ori[,to.test]))
  ## sum(a.na)
  ## sum(  is.na(info[,to.test]))
  ## sum(  is.na(info.ori[,to.test]))
  ## sum(  info[!a.na,to.test]!=info.ori[!a.na,to.test]) 
  
  ####  Alternate mulitcore fails
  ##   library(foreach)  ## could not get this to work well
  ##    library(doMC)
  ## registerDoMC(cores=7)
  ## foreach(ifi=1:length(info),.combine=c,.multicombine=TRUE,.inorder=TRUE) %dopar% process.info.list(info[[ifi]],info.types)
  ###
  
  
  ########### to process using arrays use this method timing:  73.070s  ## may be abl to make multicore  
  ## system.time(
  ## {
  ## info<- apply(as.matrix(indels[,"INFO"]),1,process.info,info.types)
  ## info<-t(info) # info[1:5,]
  ## colnames(info)<-info.types
  ## }
  ## )
  # info[1:5,] # info contains the common metrics for that location for all samples 
  ############################################
  
  #################### Process the FORMAT sectionfor each sample
  ### assumes samples are listed after FORMAT COLUMN
  info.sample<-{}
  format.posn<-match("FORMAT",colnames(indels))
  
  samples.processing<-length(column.labels)-format.posn # number of sample in the vcf file
  samples.order.in.ALL<-column.labels[(format.posn+1):length(column.labels)] # samlple labels
  
  #system.time({
  format.string<-indels[1,format.posn]
  complex.format<-sum(indels[,format.posn]!=format.string) !=0  # TRUE if all format are not the same order
  format.string<-unlist(strsplit(format.string,split=":")) # used for column names below
  
  if(complex.format){
    format.string<-indels[,format.posn]
    print("Complex Format found in vcf file")
  }else{
    format.types<-format.string ## have same format labels but make sure in the correct order
  } #use all possible attributes the function will then give the correct order
  
  
  
  
  #  test2<-info.parse(indels[,c(format.posn+1)],format.types,format.string,complex.format) # protoy
  # iter(indels[,c((format.posn+1):(format.posn+samples.processing))], by='col')
  
  
  ########Multicore run with foreach ## Done columnwise and Format the same for each column
  print(paste("START processing",samples.processing,"samples", " in, ", dim(indels), " , file ",the.sample.file,sep=" "))  
  info.sample<-foreach(info2=iter(indels[,c((format.posn+1):(format.posn+samples.processing))], by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% info.parse(info2,format.types,format.string,complex.format)
  ###
  print(paste("Finished processing",samples.processing,"samples",sep=" "))  
  
  ###construct the column names
  
  if(vcf.type=="v3") {  
    #  if(grepl(combined.extension,the.sample.file)){
    sample.class<-as.character( t(matrix(data=rep(samples.order.in.ALL,times=length(c("GT"))),nrow=length(samples.order.in.ALL),ncol=length(c("GT")),byrow=FALSE)) )
    colnames(info.sample)<-paste(sample.class,c("GT"),sep=".")
    #}
    
    #if(!grepl(combined.extension,the.sample.file)){ colnames(info.sample)<-c("GT")  }
    
  }else{ # modern vcf 4.0 format
    
    #if(grepl(combined.extension,the.sample.file)){
    sample.class<-as.character( t(matrix(data=rep(samples.order.in.ALL,times=length(format.types)),nrow=length(samples.order.in.ALL),ncol=length(format.types),byrow=FALSE)) )
    colnames(info.sample)<-paste(sample.class,format.types,sep=".")
    #}
    
    #if(!grepl(combined.extension,the.sample.file)){ colnames(info.sample)<-format.types  }
    
  }
  
  #################### FINISHED processing individual filters   indels<-indels.ori[c(1000,grep("1/2",indels.ori[,10])[1:3], grep("/3",indels.ori[,10])) ,]
  TYPE<-rep(names(the.sample.file),times=dim(indels)[1])
  filter.posn<-match("FILTER",colnames(indels))
  
  indels<-cbind(subset(indels,select=1:filter.posn),TYPE,info,info.sample) ## use subset here incase indels has only 1 element
  
  
  rm(TYPE)
  rm(info)
  rm(info.sample)
  #colnames(data)<-c( "chr","start","end","strand",c("ID","REF","ALT","QUAL","FILTER"),"TYPE",colnames(info),colnames(info.sample))
  ############################### MAKE the final data file
  print("Have data now fix multi-alleles")
  
  # the.gen<-expand.labels.to.samples(c("GT"),samples.order.in.ALL)
  ############################## multi-allele genotypes need to fix
  # poly position typicaly have no SB or HRun given
  
  #####################################################################uses just indels from this point ###################################
  
  ########Multicore run with foreach  #### WARNING PROBLEM IS LESS THAN num
  if(dim(indels)[1]<1000){num.bits<-1}else{num.bits<-num.cores}  #### WARNING PROBLEM IS LESS indels than num.cores
  print(paste("START processing muli-alleles",samples.processing,"samples", " in, ", dim(indels),sep=" "))  
  dim(indels)
  ## indels<-foreach(indels.bit=iter(indels, by='row',chunksize=as.integer(dim(indels)[1]/num.bits)), .combine=rbind,.multicombine=TRUE,.inorder=TRUE) %dopar% correct.multi.alleles.flattern(indels.bit,samples.order.in.ALL)
  #test<-correct.multi.alleles(indels,samples.order.in.ALL)
  ###
  
  indels<-foreach(indels.bit=iter(indels, by='row',chunksize=as.integer(dim(indels)[1]/num.bits)), .combine=rbind,.multicombine=TRUE,.inorder=TRUE) %dopar% correct.multi.alleles(indels.bit,samples.order.in.ALL)
  ###
  print(paste("Finished processing muli-alleles",samples.processing,"samples",sep=" "))
  
  
  ########################### PUT ALLELES IN ANNOVAR FORMAT from convert2annovar.pl line 1083########
  
  indels[is.na(indels[,"REF"]),"REF"]<-"NA" ## see a REF =N ALT= NA R confuses with "NA" with NA - kills substr below
  indels[is.na(indels[,"ALT"]),"ALT"]<-"NA"
  
  
  ref.length<-nchar(as.character(indels[,"REF"]))
  alt.length<-nchar(as.character(indels[,"ALT"]))
  is.snp<-(ref.length==1 & alt.length==1)
  
  # sum(!is.snp)
  
  POS.end<-as.numeric(indels[,"POS"])
  del<-ref.length > alt.length
  ins<-(ref.length <= alt.length) & !is.snp
  ## indels[del,][1:5,]
  ## POS.end[del][1:5]
  ### deletion or block substitution
  head<-substr(as.character(indels[del,"REF"]),1,alt.length[del])
  head.is.mut<-(head==as.character(indels[del,"ALT"]))
  indels[del,"REF"][head.is.mut]<-substr(as.character(indels[del,"REF"][head.is.mut]),(alt.length[del][head.is.mut]+1),ref.length[del][head.is.mut])
  indels[del,"ALT"][head.is.mut]<-"-"
  indels[del,"POS"][head.is.mut]<-as.numeric(indels[del,"POS"][head.is.mut]) + nchar(as.character(head[head.is.mut]))
  POS.end[del]<-POS.end[del]+ref.length[del]-1  # same for both head is mut and not head is mut
  
  ## indels
  ## POS.end
  ### insertion or block substitution
  head<-substr(as.character(indels[ins,"ALT"]),1,ref.length[ins])
  head.is.ref<-(head==as.character(indels[ins,"REF"]))
  indels[ins,"ALT"][head.is.ref]<-substr(as.character(indels[ins,"ALT"][head.is.ref]),(ref.length[ins][head.is.ref]+1),alt.length[ins][head.is.ref])
  indels[ins,"REF"][head.is.ref]<-"-"
  indels[ins,"POS"][head.is.ref]<-as.numeric(indels[ins,"POS"][head.is.ref]) + ref.length[ins][head.is.ref]-1
  POS.end[ins]<-POS.end[ins]+ref.length[ins]-1
  ########################################################################################################
  #$$$  all.alleles
  
  
  POS.end<-as.integer(POS.end) ### else can get 8e+7 type errors
  begining.cols<-c("chr","POS")
  strand<-rep("+",times=dim(indels)[1])
  ## indels<-cbind(indels[,begining.cols],POS.end,strand,indels[,colnames(indels)[!(colnames(indels) %in% begining.cols)]])
  indels<-cbind(subset(indels,select=begining.cols),POS.end,strand,subset(indels,select=colnames(indels)[!(colnames(indels) %in% begining.cols)])) ##use incase indels 1 element long
  
  
  colnames(indels)[colnames(indels)=="POS.end"]<-"end"
  colnames(indels)[colnames(indels)=="POS"]<-"start"
  
  rm(POS.end)
  rm(strand)
  indels[indels[,"chr"]=="chrMT","chr"]<-"chrM"  ### annotation algorithms expect MT
  indels
  
} # end function  process.GATK.indels

#Build keys
build.key<-function(table,key.cols,add.chr.label=FALSE){
  options(scipen=300)
  if(is.null(dim(table))){table<-as.matrix(table)} # in casea vector sent
  if(length(key.cols)<1){print("FAIL no keys columns specified");key<-1:dim(table)[1]}else{
    for (i in 1:length(key.cols)){
      if(i==1){key<-table[,key.cols[i]]}else{
        key<-paste(key,table[,key.cols[i]],sep=":")}
    }}
  if(add.chr.label){key<-paste("chr",key,sep="")}
  key}
#End of build.keys

#genotype summary
genotype.summary<-function(x){
  if(is.null(dim(x))){x<-matrix(data=x,nrow=1,ncol=length(x))}
  hetero.geno<-x=="0/1"
  ref.geno<- x=="0/0"
  alt.geno<-x=="1/1"
  alt.homo.counts<-apply(alt.geno,1,function(xx) sum(xx,na.rm=TRUE))
  hetero.counts<-apply(hetero.geno,1,function(xx) sum(xx,na.rm=TRUE))
  ref.homo.counts<-apply(ref.geno,1,function(xx) sum(xx,na.rm=TRUE))
  alt.counts<- alt.homo.counts*2+hetero.counts
  ref.counts<- ref.homo.counts*2+hetero.counts
  total.alleles<-ref.counts+alt.counts
  maf<-alt.counts/total.alleles
  genos.sum<-paste(alt.homo.counts,hetero.counts,ref.homo.counts,sep=",")
  cbind(signif(maf, digits = 6),alt.counts,ref.counts, total.alleles,(2*dim(x)[2]-total.alleles),alt.homo.counts,hetero.counts,genos.sum)
}


genotype.summary.with.apply<-function(x){
  counts<-tapply(x,x,length)
  missing.counts<-counts["NA"]
  genos<-counts[c("1/1","0/1","0/0")]
  genos[is.na(genos)]<-0 # case does not have that genotype
  alt.counts<-genos[1]*2+genos[2]
  ref.counts<-genos[2]+genos[3]*2
  total.alleles<-ref.counts+alt.counts
  maf<-alt.counts/total.alleles
  genos.sum<-paste(genos,collapse=",")
  c(signif(maf, digits = 6),alt.counts,ref.counts,total.alleles,missing.counts,genos[1],genos[2],genos.sum)
}

#Position Filtering

position.quality.filter<-function(indels,quality.cut,quality.type,quality.dirn){
  
  ##### Missing or NA values are ALWAYS set to FALSE
  ## quality.cut
  ## quality.type
  ## quality.dirn
  ################################### only test attributes that are avaialable
  avail<-names(quality.cut) %in% colnames(indels)
  print(paste("Attributes not in indels file:",toString(names(quality.cut)[!avail]),sep=" "))
  quality.cut<-quality.cut[avail]
  quality.type<-quality.type[avail]
  quality.dirn<-quality.dirn[avail]
  
  
  quality.names<- names(quality.cut) ## qaulity.names is the column name in indel indels[1:10,quality.names[j]]
  quality.measure<-names(quality.type) ## is the filter name :  X1000753W.Aff.DP_Low
  ## quality.thresh<-data.frame(key=key.indels,stringsAsFactors=FALSE)
  quality.thresh<-matrix(data=TRUE,nrow=dim(indels)[1],ncol=length(quality.cut))
  colnames(quality.thresh)<-quality.measure
  #quality.thresh[1:5,]
  
  for(j in 1:length(quality.cut)){
    if(quality.type[j]=="factor"){
      
      if(quality.dirn[j]=="exact"){
        quality.thresh[,quality.measure[j]] <- indels[,quality.names[j]]==quality.cut[j]
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## "J"=NA -> NA so NEEDED
      }
      
      if(quality.dirn[j]=="starts_with"){
        quality.thresh[,quality.measure[j]] <- grepl(paste("^",quality.cut[j],sep=""),indels[,quality.names[j]]) 
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  #  grepl("J",NA) -> FALSE so no really needed
      }
      
      if(quality.dirn[j]=="not_starts_with"){
        quality.thresh[,quality.measure[j]] <- grepl(paste("^",quality.cut[j],sep=""),indels[,quality.names[j]]) 
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  #  grepl("J",NA) -> FALSE so no really needed
      }
      
      if(quality.dirn[j]=="ends_with"){
        quality.thresh[,quality.measure[j]] <- grepl(paste(quality.cut[j],"$",sep=""),indels[,quality.names[j]]) 
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
      }
      
      if(quality.dirn[j]=="not_ends_with"){
        quality.thresh[,quality.measure[j]] <- !grepl(paste(quality.cut[j],"$",sep=""),indels[,quality.names[j]]) 
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
      }
      
      if(quality.dirn[j]=="contains"){
        quality.thresh[,quality.measure[j]] <- grepl(quality.cut[j],indels[,quality.names[j]]) 
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
      }
      
      if(quality.dirn[j]=="not_contains"){
        quality.thresh[,quality.measure[j]] <-  !grepl(quality.cut[j],indels[,quality.names[j]]) 
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE # grepl("J",NA) -> FALSE so no really needed
        quality.thresh[is.na(indels[,quality.names[j]]),quality.measure[j]] <- FALSE  #resquires additional check NA->FALSE ; FALSE->TRUE!!
      }
      
      #  testing agaisnst more than one item c("CAT","DOG") ;-quality.cut[j] would need to be a list for this to work  properly
      ## if(quality.dirn[j]=="contains"){
      ##                     quality.thresh[,quality.measure[j]] <- indels[,quality.names[j]] %in% quality.cut[j] #  quality.cut[j] would need to be alist for this to work
      ##                     quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE #"J" %in% NA -> FALSE so not really needed
      ##                     }
      
      ##     if(quality.dirn[j]=="not_contains"){
      ##                     quality.thresh[,quality.measure[j]] <- !(indels[,quality.names[j]] %in% quality.cut[j])
      ##                     quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE #"J" %in% NA -> FALSE so not really needed
      ##                     quality.thresh[is.na(indels[,quality.names[j]]),quality.measure[j]] <- FALSE  #resquires additional check NA->FALSE ; FALSE->TRUE!!
      ##                     }
      
    }     ## factor character comparison NA gets FASLE automaticaly
    if(quality.type[j]=="numeric"){
      
      if(quality.dirn[j]=="greater"){
        quality.thresh[,quality.measure[j]] <- as.numeric( as.character(indels[,quality.names[j] ]) ) >= as.numeric(as.character(quality.cut[j])) ## get warning here of NAs intro by coersion
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
      }
      if(quality.dirn[j]=="less"){
        quality.thresh[,quality.measure[j]] <- as.numeric( as.character(indels[,quality.names[j] ]) ) <  as.numeric(as.character(quality.cut[j])) ## get warning here of NAs intro by coersion
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
      }
      if(quality.dirn[j]=="exact"){
        quality.thresh[,quality.measure[j]] <- as.numeric( as.character(indels[,quality.names[j] ]) ) ==  as.numeric(as.character(quality.cut[j])) ## get warning here of NAs intro by coersion
        quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
      }
      
      
    }
    #  print(paste("passing filter",names(quality.type)[j],":",sum(quality.thresh[,quality.measure[j]] ),sep=" "))
  }
  
  ## sum(test[[paste("All",names(quality.type)[j],sep="::")]] != quality.thresh[,names(quality.type)[j]])
  quality.thresh
  
} # position.quality.filter subroutine end

#quality.thresh[1:5,]



# chr2:209113113:209113113
## indels.bit<-indels[1:4442,]

correct.multi.alleles<-function(indels.bit,samples.order.in.ALL){ ## add a flattern genotype as well to expanded indels
  
  #####################################################################uses just indels from this point ###################################
  ######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
  ## cols.wanted<-colnames(indels.bit) %in% core.ann | grepl(".GT$",colnames(indels.bit))
  ## print(class(indels.bit))
  ## print(dim(indels.bit))
  ##    indels.bit[6607:6614,cols.wanted]
  add.flattern<-TRUE
  
  # print(dim(indels.bit))
  the.cols<-colnames(indels.bit)
  has.comma.ref<-grepl(",",indels.bit[,"REF"])
  if(sum(has.comma.ref)>0){print("ERROR comma in REF allele")}
  
  ###### expand indels.bit to make room for new data
  
  alt.list<-strsplit(indels.bit[,"ALT"],split=",")
  number.of.alleles<-unlist(lapply(alt.list,length))
  #  number.of.alleles<- ((number.of.alleles+1)*(number.of.alleles+2)/2)- (number.of.alleles+1)  # extra +1 to count REf allele
  if(add.flattern){number.of.alleles[number.of.alleles>1]<-number.of.alleles[number.of.alleles>1]+1} ### so can have the collapsed one too
  
  flat.index<-rep(1:length(number.of.alleles),times=number.of.alleles)
  indels.bit.ori.size<-dim(indels.bit)[1]
  indels.bit<-indels.bit[flat.index,]
  
  if(is.null(dim(indels.bit))){ # rare case id indels.bit is only one element long so above line makes a col
    the.cols<-names(indels.bit)
    indels.bit<-matrix(data=indels.bit,nrow=1,ncol=length(indels.bit))
    colnames(indels.bit)<-the.cols
  }
  ################### fix the genotypes
  
  
  if(dim(indels.bit)[1]> indels.bit.ori.size){
    null.genotype<-{} # posns in indels.bit where there are no genotypes  
    ## to.fix<-  c("AD","GT","AD")
    fix.posns<-tapply(flat.index,flat.index,length)
    fix.posns<-fix.posns[fix.posns>1]
    location.in.data<-match(names(fix.posns),flat.index)
    names(fix.posns)<-location.in.data
    # fix.posns
    ##  ## Names(fix.posn) give the index location in data to strat and
    
    ## for(ifix in 1:length(to.fix)){ ##loop over items to fix then put into cols.wanted
    ##   a.string<-to.fix[ifix]
    ##   if(a.string %in% format.types){ cols.wanted<-paste(samples.order.in.ALL,a.string,sep=".")}
    ##   if(a.string %in% info.types ){cols.wanted<-a.string}
    ##   cols.wanted indels.bit[10847:10859,1:6]
    
    # print(paste("fixing",length(fix.posns),"poly morphic sites",sep=" "))
    ## ipos<-8  grep("4799",names(fix.posns)) grep("rs141819080",indels.bit[,"ID"])
    for(ipos in 1:length(fix.posns)){
      # print(ipos)
      a.posn<-as.integer(names(fix.posns[ipos]))
      num<-fix.posns[ipos]
      
      alleles <- unlist(strsplit(c(paste(indels.bit[a.posn,"REF"],indels.bit[a.posn,"ALT"], sep=",")),','))
      names(alleles)<-0:(length(alleles)-1)
      ## alleles
      ACs <- unlist(strsplit(indels.bit[a.posn,"AC"],','))
      names(ACs)<-1:(length(alleles)-1) # one AC for each ALT allele
      AFs <- unlist(strsplit(indels.bit[a.posn,"AF"],','))
      names(AFs)<-1:(length(alleles)-1) # one AF for each ALT allele
      
      
      
      ############### do AD for all samples are this is the sample for each posn, just replicate once
      ############### sort out the ref and alt alleles 1/1 could be from 0,1 OR 1,2
      ###############  use the PL to decide and test
      
      cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
      the.genotypes<-indels.bit[a.posn,cols.wanted]
      the.alleles<-the.genotypes
      the.alleles[the.alleles=="NA" | is.na(the.alleles) ]<-"0/0" #just so next step does not overlop the.genotypes==NA are filled with 0,0,0 AD
      the.alleles<-strsplit(the.alleles,split="/")
      
      ref.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[1]))) # +1 # COLUMN POSN IN AD.LIST
      ## ref.allele.ori<-ref.allele
      alt.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[2]))) # +1 # COLUMN POSN IN AD.LIST
      
      ##################
      ## ref.allele
      ## ref.allele.ori
      ##  alt.allele
      ## a.posn
      ## indels.bit[a.posn,cols.wanted]
      
      
      ############################### define the PL array so can use to get most likely alleles for ties
      
      PL.combo<-matrix(data=NA,nrow=length(alleles),ncol=length(alleles))
      ############ set up the PL structure
      for(iref in 1:length(alleles)){
        for(ialt in iref:length(alleles)){
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          PL.combo[iref,ialt]<-paste(iref-1,ialt-1,sep=":") 
        }}
      PL.combo<-as.character(PL.combo)
      PL.combo<-PL.combo[!is.na(PL.combo)]
      ## PL.combo ### the v4 samtools 
      
      cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
      the.PL<- indels.bit[a.posn,cols.wanted]
      the.PL[the.PL=="NA" | is.na(the.PL)]<-paste(rep("3000",times=length(PL.combo)),collapse=",")
      pl.list<-strsplit(the.PL,split=",")
      pl.list<-pl.list[cols.wanted]
      pl.list.length<-length(pl.list)
      pl.list<-as.integer(unlist(pl.list))
      dim(pl.list)<-c(length(PL.combo), pl.list.length)
      pl.list<-t(pl.list)
      colnames(pl.list)<-PL.combo ## p.list and array of sample by allele combos 
      
      
      
      
      ######### fix location where have   1/1 or 2/2 or 3/3 etc choose two most probable REF allele 0 if a tie..
      ### no loger need this the ref allele is always 0 now
      ref.allele<-rep(0,times=length(ref.allele))
      
      ## ties<-which(ref.allele==alt.allele)
      ## pl.ties<-pl.list[ties,] ## P
      ## if(length(ties)>0){  ## just in case is no ties which causes an error
      ##   # ir<-1
      ## for(ir in 1:length(ties)){
      ##   if(is.null(dim(pl.ties))){pl.ties<-as.matrix(pl.ties);pl.ties<-t(pl.ties)} ## in case inly one ties and get a vector!
      ##   best.pl<-sort(pl.ties[ir,]) ## this sort works as 0:1 reported before 1:2 is they are tied (see PL.combo order)
      ##   the.alleles<-{}
      ##   while(length(the.alleles)<2){
      ##     the.alleles<-sort(as.integer(unique(unlist(strsplit(paste(names(best.pl[1:2]),collapse=":"),split=":")))))
      ##   }
      ##   ref.allele[ties[ir]]<-the.alleles[1]
      ##   alt.allele[ties[ir]]<-the.alleles[2]
      
      ## }}
      ###################################################
      ## ref.allele ### these contain the best choese of alleles fo the genotypes listed
      ## alt.allele  ### these contain the best choese of alleles fo the genotypes listed 
      ############### do AD for all samples are this is the sample for each posn, just replicate once
      
      
      
      
      #        cols.wanted.AD<-gsub(".GT$",".AD",cols.wanted)
      #        cols.wanted.DP<-gsub(".GT$",".DP",cols.wanted)
      #        for(i in 1:length(cols.wanted)){
      #        test<-(indels[,cols.wanted[i]]=="0/0" & indels[,cols.wanted.DP[i]]==".")
      #            print(sum(test))
      #     }
      #        indels[test,c(colnames(indels)[1:10],cols.wanted[i],cols.wanted.AD[i],cols.wanted.DP[i])][1:5,]
      #            indels[test,c(colnames(indels)[1:10],cols.wanted[i],cols.wanted.AD[i],cols.wanted.DP[i])]
      #        
      #        [1:5,]
      
      #####get The AD 
      cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
      the.AN<- indels.bit[a.posn,cols.wanted]
      ad.list<-strsplit(the.AN,split=",")
      AD.alleles<-max(as.integer(unique(unlist(lapply(ad.list,length)))))
      
      the.AN[the.AN=="NA" | is.na(the.AN) | the.AN=="."]<-paste(rep("0",times= AD.alleles),collapse=",")
      ## Mhairi to do - replace length(alleles) with nrows<-length(ad.list[[1]])
      
      ad.list<-strsplit(the.AN,split=",")
      ad.list<-ad.list[cols.wanted]
      ad.list.length<-length(ad.list)
      ad.list<-as.integer(unlist(ad.list))
      
      dim(ad.list)<-c(AD.alleles, ad.list.length)
      ad.list<-t(ad.list)
      ad.list<-ad.list[,1:length(alleles)]
      
      ref.counts<-rep(NA,times=dim(ad.list)[1])
      alt.counts<-rep(NA,times=dim(ad.list)[1])
      for(ir in 1:dim(ad.list)[1]){
        ref.counts[ir]<-ad.list[ir,(ref.allele[ir]+1)] ## ad list is ordered in assending allele
        alt.counts[ir]<-ad.list[ir,(alt.allele[ir]+1)]
      }
      the.AN<-paste(ref.counts,alt.counts,sep=",")
      names(the.AN)<-cols.wanted
      
      #################################################################
      
      #####get The PL
      cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
      the.PL<-rep(NA,times=dim(pl.list)[1])
      for(ir in 1:dim(pl.list)[1]){
        pl.wanted<-c(paste(ref.allele[ir],ref.allele[ir],sep=":"),
                     paste(ref.allele[ir],alt.allele[ir],sep=":"),
                     paste(alt.allele[ir],alt.allele[ir],sep=":"))
        
        the.PL[ir]<-paste(pl.list[ir,pl.wanted],collapse=":")
      }
      names(the.PL)<-cols.wanted
      
      #################################################################
      
      ## the.genotypes
      ## ACs
      ## AFs
      ## the.AN
      ## the.PL
      ## ref.allele
      ## ref.allele.ori
      ## alt.allele
      ## ad.list[1:5,]
      ## pl.list[1:5,]
      
      ## print(alleles)
      
      
      
      ###########################################################
      iplace<-1 #################### now fix all positions that start with a.posn iplace goes from 1:num
      
      if(add.flattern){
        iref<-0
        ialt<-1
        ## for(iref in 0:(length(alleles)-1)){
        ##   for(ialt in iref:(length(alleles)-1)){
        ##     if(iref==ialt){next}
        i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
        i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
        ##    print(paste(i.ref,i.alt,sep=":"))
        ##   }}
        
        a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change
        
        ## print(a.place)
        ## indels.bit[a.place,cols.wanted]
        ## test[a.place,cols.wanted]
        ## indels.bit[a.place,1:10]
        
        
        indels.bit[a.place,"REF"]<-alleles[i.ref]
        indels.bit[a.place,"ALT"]<-alleles[i.alt]
        indels.bit[a.place,"AC"]<-sum(as.numeric(ACs),na.rm=TRUE) # 
        indels.bit[a.place,"AF"]<-sum(as.numeric(AFs),na.rm=TRUE) # AFs[i.alt]
        
        if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
        if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
        if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
        if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
        
        indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],"flat",sep=":")
        
        ########################################  FIX GT ################################
        
        new.gt<-the.genotypes
        target.GT<-!is.na(new.gt) | new.gt!="NA"
        
        
        new.gt[!target.GT]<-"NA"
        if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
        
        
        ######### here convert all possible alternative alleles to one ( 2 stage so can handle 1/10 11/15 etc:
        to.one<-paste("^[",paste(names(alleles)[!(names(alleles) %in% c(i.ref))],collapse=" "),"]/",sep="") ## everying to aly except  0
        new.gt<-gsub(to.one,"1/",new.gt)  #convect all left alleles to zero
        to.one<-paste("/[",paste(names(alleles)[!(names(alleles) %in% c(i.ref))],collapse=" "),"]$",sep="") ## everying to aly except  0
        new.gt<-gsub(to.one,"/1",new.gt)  #convect all right alleles to zero
        
        new.gt<-gsub(i.alt,"1",new.gt)
        new.gt[new.gt=="1/0"]<-"0/1"
        
        
        cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
        indels.bit[a.place,cols.wanted]<-new.gt
        
        ##################################################################################
        
        ########################################  FIX AD ################################
        
        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
        
        ##################################################################################
        
        ########################################  FIX PL ################################
        
        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
        
        ##################################################################################
        ################################# step one out to get the others
        
        
        
        iplace<-2 #################### now fix all positions that start with a.posn iplace goes from 1:num
        
      } # end add flattern if added flatten iplace is 2 otherwise is 1
      
      ## iref<-0
      ## ialt<-2
      for(iref in 0:0){  ## iref==0 always as now always have the reference alleles
        for(ialt in iref:(length(alleles)-1)){
          if(iref==ialt){next}
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          ## print(paste(i.ref,i.alt,sep=":"))
          ##  }}
          
          a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change
          
          ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:35]
          
          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-ACs[i.alt]
          indels.bit[a.place,"AF"]<-AFs[i.alt]
          
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],sep=":")
          
          ########################################  FIX GT ################################
          
          new.gt<-the.genotypes
          missing.genotype<-is.na(new.gt) | new.gt=="NA"
          target.GT<- !missing.genotype
          
          ## new.gt[!target.GT & !missing.genotype]<-"0/0"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          
          
          ## new.gt<-gsub(i.alt,"1",new.gt) ## set the 0/2 to 0/1 --- fails in case where have 2/3 for i.alt=2
          
          ## new.gt<-gsub(i.ref,"0",new.gt) ## set the other alleles to reference now refernce alwasy zero
          ## Need a 2 step process for alleles like 0/10 11/11 etc 
          to.zero<-paste("^[",paste(names(alleles)[!(names(alleles) %in% c(i.alt))],collapse=" "),"]/",sep="") ## everying to aly except  0
          new.gt<-gsub(to.zero,"0/",new.gt)  #convect all left alleles to zero
          to.zero<-paste("/[",paste(names(alleles)[!(names(alleles) %in% c(i.alt))],collapse=" "),"]$",sep="") ## everying to aly except  0
          new.gt<-gsub(to.zero,"/0",new.gt)  #convect all right alleles to zero
          
          new.gt<-gsub(i.alt,"1",new.gt)
          new.gt[new.gt=="1/0"]<-"0/1"
          
          
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt
          
          ##################################################################################
          
          ########################################  FIX AD ################################
          
          cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
          new.AD<-the.AN
          new.AD[!target.GT]<-"NA"
          indels.bit[a.place,cols.wanted]<-new.AD
          
          ##################################################################################
          
          ########################################  FIX PL ################################
          
          cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
          new.PL<-the.PL
          new.PL[!target.GT]<-"NA"
          indels.bit[a.place,cols.wanted]<-new.PL
          
          ##################################################################################
          
          
          iplace<-iplace+1
        }} #iref ialt iplace fill in same genotypes for location
      
    } #ipos new location
    
    #} #ifix
    if(length(null.genotype)>1){  indels.bit<-indels.bit[-1*null.genotype,]}
  } # indels.bit.ori.size has to do work and fix alleles
  
  
  
  #####################################################################uses just indels.bit from this point ###################################
  
  indels.bit
  
  
  
}  ## end function correct.muli.alleles

############## map 1 to many relationship
############## map 1 to many relationship



indentify.IN.repeat<-function(indel,looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19"){
  # requires "chr","start","end","REF","ALT","TYPE" columns
  
  require(genome,character.only=TRUE)
  require("Biostrings")
  require("IRanges")
  #Hsapiens
  
  ## test<-indel[chk.in.repeat,1:10]
  ## test["chr5:172385267:172385267:-:AA:indel",]
  ## indel<-a.indel[chk.in.repeat,]
  ## print(bases.about)
  
  
  if( (bases.about< (homo.run.max+1)) | (bases.about< (di.run.max*2)) ){bases.about<-max(c((homo.run.max+1),(di.run.max*2)))} # look to look past the length of the repeart run
  ## chk.in.repeat<-large.indel & !are.repeats
  ## sum(chk.in.repeat)
  
  
  if(is.null(rownames(indel))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;key.small<-build.key(indel,core.ann,add.chr.label=FALSE)}else{key.small<-rownames(indel)}
  the.chrom<-indel[,"chr"]
  
  ### look forward
  if(looking=="forward"){
    starts<-as.numeric(indel[,"end"]) 
    ends<-as.numeric(indel[,"end"]) +   bases.about
  }else{ ## looking back
    starts<-as.numeric(indel[,"start"]) -  bases.about 
    ends<-as.numeric(indel[,"start"])   
  }
  ##### stay within the existing region 
  ## starts<-pmax.int(starts,peaks[,"start"])
  ## ends<-pmin.int(ends,peaks[,"end"])
  ############################ found if I did not do this I got a single extra region  15496:  16 89084207 89084461  Krtap16-7
  
  
  large.indels.test<-getSeq(Hsapiens, the.chrom, starts, ends) # large.indels.test<-all.genomic
  ## length(all.genomic)
  ## system.time(x<- XStringViews(all.genomic, "DNAString"))
  ## x.labels<-paste(the.chrom,starts,ends,sep=":")
  ## names(x)<-x.labels
  
  
  
  di<-dinucleotideFrequency(large.indels.test,step=1)
  widths<-width(large.indels.test)
  #homo.run<-pmax(di[,colnames(di) %in% c("AA","TT","GG","CC")])
  homo.run<-pmax(di[,"AA"],di[,"TT"],di[,"GG"],di[,"CC"])
  homo.run[homo.run>0]<-homo.run[homo.run>0]+1  
  #di.run<-pmax(di[,colnames(di) %in% c("AC","AG","AT","CG","CT","GT")])
  di.run<-pmax(di[,"AC"],di[,"AG"],di[,"AT"],di[,"CG"],di[,"CT"],di[,"GT"])
  di.run<-di.run*2 # number of bases in repeats
  
  ## di.run[1:50]
  ## homo.run[1:50]
  ## widths[1:50]
  
  repeats<-(di.run==widths & (di.run/2)>=di.run.max) | (homo.run>=(widths-1) &  homo.run>=homo.run.max)
  
  
  ## repeats[1:5]
  ## grep(TRUE,repeats)[1:5]
  ## large.indels.test[grep(TRUE,repeats)[1:5]]
  ## length(repeats)
  ## sum(repeats)
  ## key.small[repeats][1:5]
  ## homo.run[repeats][1:5]
  ## di.run[repeats][1:5]
  
  return(repeats)
  
  ## remove.repeats<-key.small[repeats]
  ## length(remove.repeats)
  ## remove.repeats[1:20]
  
  ## are.in.repeats.forward<- key %in% remove.repeats
  
  ## are.in.repeats.forward[1:20]
  ## sum(are.in.repeats.forward)
  
} # indentify.IN.repeat

###############################################################
## indel<-a.indel
## di.run.max=3
## homo.run.max=5


## a.indel[grep("REF",a.indel[,"REF"])-1,1:10]

identify.repeats<-function(indel,di.run.max=3,homo.run.max=5){
  # requires "chr","start","end","REF","ALT","TYPE" columns 
  require(Biostrings)
  ## poly.morphic<-grepl(":\\d*$",key,perl=T)
  ## key[poly.morphic][1:50]
  if(is.null(rownames(indel))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;key<-build.key(indel,core.ann,add.chr.label=FALSE)}else{key<-rownames(indel)}
  REF.length<-nchar(as.character(indel[,"REF"]))
  ALT.length<-nchar(as.character(indel[,"ALT"]))
  
  large.indel<-REF.length>1 | ALT.length>1
  
  #sum(large.indel)
  
  the.large.indels<-indel[,"REF"][large.indel]
  use.alt<-ALT.length[large.indel]>REF.length[large.indel]
  #the.large.indels[60:65]
  the.large.indels[use.alt]<-indel[,"ALT"][large.indel][use.alt]
  
  
  ## the.large.indels[grep("E",the.large.indels)][1:5]
  
  
  large.indels.test<-DNAStringSet(the.large.indels)
  #large.indels.test[1:50]
  di<-dinucleotideFrequency(large.indels.test,step=1)
  widths<-width(large.indels.test)
  #homo.run<-pmax(di[,colnames(di) %in% c("AA","TT","GG","CC")])
  homo.run<-pmax(di[,"AA"],di[,"TT"],di[,"GG"],di[,"CC"])
  homo.run[homo.run>0]<-homo.run[homo.run>0]+1  
  #di.run<-pmax(di[,colnames(di) %in% c("AC","AG","AT","CG","CT","GT")])
  di.run<-pmax(di[,"AC"],di[,"AG"],di[,"AT"],di[,"CG"],di[,"CT"],di[,"GT"])
  di.run<-di.run*2 # number of bases in repeats
  
  ## di.run[1:50]
  ## homo.run[1:50]
  ## widths[1:50]
  
  repeats<-(di.run==widths & (di.run/2)>=di.run.max) | (homo.run>=(widths-1) &  homo.run>=homo.run.max)
  
  ## length(repeats)
  ## sum(repeats)
  ## the.large.indels[repeats][1:5]
  ## homo.run[repeats][1:5]
  ## di.run[repeats][1:5]
  
  remove.repeats<-key[large.indel][repeats]
  ## length(remove.repeats)
  ## remove.repeats[1:50]
  
  are.repeats<- key %in% remove.repeats
  ## sum(are.repeats)
  
  ## length(large.indel) 
  ## length(are.repeats)
  dim(indel)
  
  return(are.repeats)
}



#genotypes, fam,"pheno",num.samples)
logistic.run<- function(gdata, pdata,pcolumn,num.samp) {
  ## gdata <- genotypes
  ## pdata <- fam
  ## pcolumn <- "pheno"
  ## num.samp <- num.samples
  
  
  ## print(dim(gdata))
  ## print(rownames(gdata))
  n <- dim(gdata)[1]
  ##   if(is.null(n)){ # only one element
  ##     gdata<-as.matrix(gdata)
  ##     if(dim(gdata)[1]!=1){gdata<-t(gdata)}
  ## }
  #  n <- names(gdata)
  markers <- length(n) - 1
  
  all.missing<-apply(gdata,1,function(x) sum(is.na(as.numeric(x))))
  all.missing<-grep(dim(gdata)[2],all.missing,fixed=TRUE)
  do.rows<-1:dim(gdata)[1]
  if(length(all.missing)>0){do.rows<-do.rows[!(do.rows %in% all.missing)]}
  ## all.missing[1:5]
  ## sum(missing)
  
  if(is.null(rownames(gdata))){rownames(gdata)<-1:dim(gdata)[1]}
  
  results <- matrix(data=NA,nrow=nrow(gdata),ncol=7)
  #results <- matrix(nrow=ncol(gdata),ncol=5)
  colnames(results) <- c("SNP","REFfreq","ngeno","Estimate", "Std Error", "Z", "pval" )
  #  i<-2
  for( i in 1:length(do.rows)) {
    print(i)
    
    #   for( i in 2:length(n)) {  
    # pdata[,"missing"]<-as.numeric(gdata[i,])
    ## for aogc & hbm use both
    #  fit1 <- glm(form1,data=data,family=binomial(logit)) # logistic for case/control
    #  fit <- lm(gdata[i,] ~ pdata[,"PCA1"]+ pdata[,"PCA2"] + pdata[,"PCA3"] + pdata[,"PCA4"]+pdata[,pcolumn])
    fit <- lm(as.numeric(pdata[,pcolumn])~as.numeric(gdata[do.rows[i],]) ) # same results as plink
    #  fit <- lm(as.numeric(gdata[do.rows[i],]) ~ as.numeric(pdata[,pcolumn]) )
    a <- summary(fit)
    
    ### get coeff for pheno1 only
    #     Estimate <- a$coefficients[2,1]
    #     Std_error <- a$coefficients[2,2]
    #     T <- a$coefficients[2,3]
    #     pval <- a$coefficients[6,1:4]
    the.geno<-as.numeric(gdata[do.rows[i],])
    the.geno<-the.geno[!is.na(the.geno)]
    p= sum(the.geno)/(2*length(the.geno))
    # p= 1- sum(as.numeric(gdata[do.rows[i],]),na.rm=TRUE)/(2*num.samp)
    
    if( dim(a$coefficients)[1]==1){ # no slope (monomorphic typically)  so get P if intercept - bad! as dim(a$coefficients)[1]==1 NOT 2 (slope and intercept)
      results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),rep(NA,times=4))
    }else{
      # results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
      results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients["as.numeric(gdata[do.rows[i], ])",1:4]) ## just in case ordr changes
      
    }
    
  }
  
  if(length(all.missing)>0){results[all.missing,"SNP"]<-rownames(gdata)[all.missing]}
  #rownames(results) <- seq(from = 1, to = markers)
  #colnames(results) <- c("Estimate", "Std Error", "Z", "pval" )
  
  return(results)
  
}   
 
#Calculation of Aij

###################################################################


genetic.sample.QC.accumilate<-function(genotypes,p,position.filter,all.sample.labels,on.X,parX){
  
  # position.filter<-position.filter.QC  
  # genotypes<-all.genotypes$"genotypes"
  # print(dim(genotypes))
  # print(length(p))
  
  if(is.null(dim(genotypes))){  # just in case not a matrix or data.frame
    ncol<-length(genotypes)
    nrow<-1
    dim(genotypes)<-c(nrow,ncol)
  }
  
  Aij<-matrix(data=0,nrow=(length(all.sample.labels)+length(all.sample.labels)),ncol=(length(all.sample.labels)+2+length(all.sample.labels)))
  colnames(Aij)<-c(all.sample.labels,"x.hetro","x.homo",paste(all.sample.labels,"C",sep="_"))
  rownames(Aij)<-c(all.sample.labels,paste(all.sample.labels,"C",sep="_"))
  #Aij
  print(dim(genotypes))
  for(i in 1:length(all.sample.labels)){  # columns
    for(j in i:length(all.sample.labels)){  # rows
      
      ii<-i+2+length(all.sample.labels) # columns
      jj<-j+length(all.sample.labels)  # rows
      
      ## print(paste(i,j,sep=":"))
      ## print(paste(rownames(Aij)[j],colnames(Aij)[i],sep=" : "))
      ## print(paste(rownames(Aij)[jj],colnames(Aij)[ii],sep=" : "))
      #   test=paste(rownames(Aij)[j],colnames(Aij)[i],sep=" : ") 
      
      present<-!(is.na(genotypes[,i]) | is.na(genotypes[,j]))        
      ok<-present & position.filter
      if(sum(ok)<2){next} # not worth it move on to avoid matrix issues
      
      if(i==j){
        ## tapply(indels[ok,"chr"],indels[ok,"chr"],length)
        Aij[j,"x.hetro"]<-sum(genotypes[ok & on.X & !parX,i]==1)
        Aij[j,"x.homo"]<-sum(genotypes[ok & on.X & !parX,i]==2)
        
        
        #         ####AN:PCA use for confounding biases: ###########Start of pca biases correction  
        #         cp<-as.data.frame(genotypes)
        #         cc<-t(cp)
        #         print(paste0("dim of cc (if i=j in Aij) cc:",dim(cc)[1]," rows and ", dim(cc)[2], " columns "))
        #         cc<-as.data.frame(cc)
        #         
        #         cc<-cbind(sample.cc=rownames(cc),cc)
        #         #pca<-pca[,-1]
        #        suppressWarnings( if (!grepl(".GT$",pca[,"ParticipantCode"])==TRUE){
        #         pca[,"ParticipantCode"]<-paste0(pca[,"ParticipantCode"],".GT")
        #         })
        #         
        #         posns<-match(cc[,"sample.cc"],pca[,"ParticipantCode"])
        #         new.pca<-cbind(cc,pca[posns,])
        #         new.pca<-new.pca[,!colnames(new.pca)%in%c("ParticipantCode")]
        #         suppressWarnings(new.pca[,!(colnames(new.pca)%in%c("ParticipantCode"))] <- sapply(new.pca[,!(colnames(new.pca)%in%c("ParticipantCode"))], function (x){as.numeric(as.character(x))}))
        #         
        #         #colnames(new.pca)[!colnames(new.pca)%in%c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")]
        #           colnames(new.pca)<-gsub(":","_",colnames(new.pca))  
        #         
        #         ##if there is a column with NAs only
        #         lst <- lapply(names(new.pca)[!colnames(new.pca)%in%c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")], function(x)  {
        #           x1 <- new.pca[[x]]
        #           if(all(is.na(x1))){
        #             NA } else lm(formula(paste(x, '~', paste0("PCA", 1:4, collapse="+"))), 
        #                          data= new.pca)$fitted.values
        #         })
        #         
        #         colnames(new.pca)<- gsub("_",":",colnames(new.pca))
        #         new.pca<-new.pca[,colnames(new.pca)[!colnames(new.pca)%in%c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")]]
        #         genotypes<-t(new.pca)
        #         genotypes<-as.matrix(genotypes)
        ###########end of pca biases correction  #############
        Aijk<- (as.numeric(genotypes[ok,j])^2 -(1.0 + 2*p[ok])*as.numeric(genotypes[ok,j]) + 2*p[ok]^2)/(2*p[ok]*(1-p[ok]))
        
        
        Aij[j,i]<- sum( Aijk ,na.rm=TRUE)  #sum(Aijk,na.rm=TRUE)
        Aij[jj,ii]<- sum(!is.na(Aijk))
      }else{
        
        #         ####AN:PCA use for confounding biases: ###########Start of pca biases correction  
        #         cp<-as.data.frame(genotypes)
        #         cc<-t(cp)
        #         print(paste0("dim of cc (if i!==j in Aij) cc:",dim(cc)[1]," rows and ", dim(cc)[2], " columns "))
        #         cc<-as.data.frame(cc)
        #         
        #         cc<-cbind(sample.cc=rownames(cc),cc)
        #         #pca<-pca[,-1]
        #         suppressWarnings( if (!grepl(".GT$",pca[,"ParticipantCode"])==TRUE){
        #           pca[,"ParticipantCode"]<-paste0(pca[,"ParticipantCode"],".GT")
        #         })
        #         
        #         
        #         posns<-match(cc[,"sample.cc"],pca[,"ParticipantCode"])
        #         new.pca<-cbind(cc,pca[posns,])
        #         new.pca<-new.pca[,!colnames(new.pca)%in%c("ParticipantCode")]
        #         suppressWarnings(new.pca[,!(colnames(new.pca)%in%c("ParticipantCode"))] <- sapply(new.pca[,!(colnames(new.pca)%in%c("ParticipantCode"))], function (x){as.numeric(as.character(x))}))
        #         
        #         #colnames(new.pca)[!colnames(new.pca)%in%c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")]
        #         colnames(new.pca)<-gsub(":","_",colnames(new.pca))  
        #         
        #         ##if there is a column with NAs only
        #         lst <- lapply(names(new.pca)[!colnames(new.pca)%in%c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")], function(x)  {
        #           x1 <- new.pca[[x]]
        #           if(all(is.na(x1))){
        #             NA } else lm(formula(paste(x, '~', paste0("PCA", 1:4, collapse="+"))), 
        #                          data= new.pca)$fitted.values
        #         })
        #         
        #         colnames(new.pca)<- gsub("_",":",colnames(new.pca))
        #         new.pca<-new.pca[,colnames(new.pca)[!colnames(new.pca)%in%c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")]]
        #         genotypes<-t(new.pca)
        #         genotypes<-as.matrix(genotypes)
        #         ###########end of pca biases correction  #############
        #####
        Aijk<- (as.numeric(genotypes[ok,i])-2*p[ok])*(as.numeric(genotypes[ok,j])-2*p[ok])/(2*p[ok]*(1-p[ok]))
        Aij[j,i]<- sum( Aijk ,na.rm=TRUE)  #sum(Aijk,na.rm=TRUE)
        Aij[jj,ii]<- sum(!is.na(Aijk))
      }
      
    }}
  
  Aij
}
###################################################################



sum.QC.matrix<-function(Aij,all.sample.labels){
  
  wanted.att<-c("sample_A","sample_B","IBS","sex_Predicted","sex.code.Predicted","Het_chrX_A","ChrX-Het:Homo","Num_Good_SNPs_A")
  related<-matrix(data=NA,nrow=(length(all.sample.labels)*(length(all.sample.labels)+1)/2),ncol=length(wanted.att))
  colnames(related)<-wanted.att # set the columns of related; "sample_A","sample_B","IBS","sex_Predicted" can't be modified they are used in "Apply inheritance QC_MULTI_NEW.r "
  
  k<-1
  for(i in 1:length(all.sample.labels)){  # columns
    for(j in i:length(all.sample.labels)){  # rows
      
      ii<-i+2+length(all.sample.labels) # columns
      jj<-j+length(all.sample.labels)  # rows
      
      if(i==j){
        ## tapply(indels[ok,"chr"],indels[ok,"chr"],length)
        x.het<- Aij[j,"x.hetro"]/( Aij[j,"x.hetro"]+ Aij[j,"x.homo"])
        sex=2;sex.code="F";het.count<-paste(Aij[j,"x.hetro"],Aij[j,"x.homo"],sep=":")
        if(!is.finite(x.het)){sex=9;sex.code="Unknown"}else{
          if(x.het<0.3){sex=1;sex.code="M"}
          if(x.het<0.3){sex=1;sex.code="M"}
        }
        
        
        A<-  1.0+ Aij[j,i]/Aij[jj,ii]
        #  Aij<-  1.0+(1-(1.0/sum(!is.na(Aijk))/var(Aijk,na.rm=TRUE)))*(Aij-1.0)
        
      }else{
        
        A<-  Aij[j,i]/Aij[jj,ii]
        #  Aij<-  (1-(1.0/sum(!is.na(Aijk))/var(Aijk,na.rm=TRUE)))*Aij
      }
      
      ## sum(ok)
      ## Aijc("sample_A","sample_B","IBS","sex_Predicted","sex.code_A","Het_chrX_A","ChrX-Het:Homo,""Num_Good_SNPs_A","Num_POLYmorph_SNPs_A")
      
      related[k,] <- c(all.sample.labels[i],all.sample.labels[j],signif(A,digits=3),sex,sex.code,signif(x.het,digits=3),het.count,Aij[jj,ii])
      k<-k+1
    }}
  
  related
}

#Filtering genotypes
filtered.genotype <- function(indels,the.samples, prefix = "", suffix = "", het.read.thresh = 20, het.low = 0.2, het.high = 0.80, het.low.ref = 0.05, het.high.alt = 0.95, depth.het.low = 7, depth.homo.low = 5, min.gq = 20) {
  
  ########################### ALL ###############
  ## applied filters to genotypes and then makes a summary
  ## indels<-  a.indel #[41458:41460,]
  ## the.samples<-gsub(".GT$","",the.samples)
  ## prefix<-""
  ## suffix<-""
  ## het.read.thresh<-20  # depth must exceed this to apply filer on genotypes
  ## het.low<-0.20
  ## het.high<-0.80
  ## het.low.ref<-0.05
  ## het.high.alt<-0.95
  ## depth.het.low<-10
  ## depth.homo.low<-5
  ## min.gq<-20
  ## ## the.sample
  ## s<-fam[,1]
  #the.samples
  #if(prefix!=""){prefix<-paste(prefix,".",sep="")}
  
  # print(dim(indels))
  
  if (is.null(rownames(indels))) {
    core.ann <- c("chr", "start", "end", "REF", "ALT", "TYPE") 
    rownames(indels) <- build.key(indels, core.ann, add.chr.label = FALSE)
  }
  
  key.indels <- rownames(indels)
  
  
  allele.depths.group <- indels[, paste(the.samples, "AD", sep = ".")]
  summary.het.indiv <- apply(as.matrix(allele.depths.group), 2, allele.summary.individuals)
  colnames(summary.het.indiv) <- gsub(".AD$", ".HET", colnames(summary.het.indiv))
  
  
  if (sum(paste(the.samples, "DP", sep = ".") %in% colnames(indels)) == length(the.samples)) {
    allele.depths <- as.matrix(indels[, paste(the.samples, "DP", sep = ".")])
  } else {
    allele.depths <- apply(as.matrix(allele.depths.group), 2, allele.DP.from.AN)
    colnames(allele.depths) <- gsub(".AD", ".DP", colnames(allele.depths))
  }
  
  ncol <- length(the.samples)
  nrow <- dim(indels)[1]
  allele.depths <- as.numeric(allele.depths)
  allele.depths[is.na(allele.depths)] <- 0
  
  if (is.null(dim(allele.depths))) {
    dim(allele.depths) <- c(nrow, ncol)
    colnames(allele.depths) <- paste(the.samples, "DP", sep = ".")
  }
  
  genotypes <- as.matrix(indels[, paste(the.samples, "GT", sep = ".")])
  
  if (sum(colnames(summary.het.indiv) != paste(the.samples, "HET", sep = ".")) != 0) {
    summary.het.indiv <- summary.het.indiv[, paste(the.samples, "HET", sep = ".")]
  }
  
  if (sum(colnames(allele.depths) != paste(the.samples, "DP", sep = ".")) != 0) {
    summary.het.indiv <- summary.het.indiv[, paste(the.samples, "HET", sep = ".")]
  }
  
  ## tapply(genotypes[,1],genotypes[,1],length)
  passing.genos <- (
    (genotypes == "1/1" & (( summary.het.indiv > het.high.alt) | allele.depths < het.read.thresh) & allele.depths > depth.homo.low &  grepl("^snp", indels[, "TYPE"])) |
      (genotypes == "0/1" & ((summary.het.indiv >= het.low & summary.het.indiv <= het.high) | allele.depths < het.read.thresh)  &  allele.depths > depth.het.low &  grepl("^snp", indels[, "TYPE"])) |
      (genotypes == "0/0" & ((summary.het.indiv < het.low.ref) | allele.depths < het.read.thresh ) & allele.depths > depth.homo.low  &  grepl("^snp", indels[, "TYPE"]))  |
      (genotypes == "1/1" & allele.depths > depth.homo.low & grepl("^indel", indels[, "TYPE"])) |
      (genotypes == "0/1" & allele.depths > depth.het.low  & grepl("^indel", indels[, "TYPE"])) |
      (genotypes == "0/0" & allele.depths > depth.homo.low & grepl("^indel", indels[, "TYPE"])) |
      (genotypes == "1/1" & allele.depths > depth.homo.low & grepl("^CREST", indels[, "TYPE"])) |
      (genotypes == "0/1" & allele.depths > depth.het.low  & grepl("^CREST", indels[, "TYPE"])) |
      (genotypes == "0/0" & allele.depths > depth.homo.low & grepl("^CREST", indels[, "TYPE"]))
  )
  
  genotypes[!passing.genos] <- NA  ### exclude failing genotypes from the calculation
  
  # Apply GQ filtering. I'm (jje) not sure if this is the most efficient way to do this.
  
  GQ.column.names <- paste0(the.samples, ".GQ")
  # print(min.gq)
  if (all(GQ.column.names %in% colnames(indels))) {
    gq <- indels[, GQ.column.names]
    gq[gq == "NA"] <- NA
    # storage.mode(gq) <- "integer" ### this was failing as it was not chaning the storge mode.
    gq[gq <= min.gq] <- NA
    genotypes[is.na(gq)] <- NA
  } else {
    warning("missing GQ columns. GQ filtering NOT applied")
  }
  
  return(genotypes)
}

#Prepare VCF file to read
# sample.files.ori<-sample.files
## sample.files<- sample.files[isamp]
prepare.for.Vcf.file.read<-function(sample.files){
  print(paste("Doing samples",sample.files,sep=" "))
  setwd(eval(as.name(paste(names(sample.files),"dir",sep=".")))) ## different directory for each indel type
  
  ######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format
  
  ### get location of header:
  chromo1<-try(scan(sample.files,what=character(),n=5000,sep="\n",skip=0,fill=TRUE,na.strings="",quote="\"")) ## find the start of the vcf file
  skip.lines<-grep("^#CHROM",chromo1)
  if(length(skip.lines)>1){print("ERROR multiple chrom lables found");skip.lines<-skip.lines[1]}
  skip.lines<-skip.lines
  
  options(show.error.messages = TRUE)
  column.labels<-read.delim(sample.files,header=F,nrows=1,skip=(skip.lines-1),sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"",colClasses="character")
  num.vars<-dim(column.labels)[2]
  
  vcf.head<-read.delim(sample.files,header=F,nrows=(skip.lines-1),skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,quote="\"",colClasses="character")
  info.flag.string<-"^##INFO=<ID="
  info.types<-grep(info.flag.string,vcf.head[,1])
  info.class<-extract.value.from.format(vcf.head[info.types,1],"Type=")
  info.description<-extract.value.from.format(vcf.head[info.types,1],"Description=")  
  info.types<-extract.value.from.format(vcf.head[info.types,1],"ID=")
  
  
  format.flag.string<-"^##FORMAT=<ID="
  format.types<-grep(format.flag.string,vcf.head[,1])
  format.class<-extract.value.from.format(vcf.head[format.types,1],"Type=")
  format.description<-extract.value.from.format(vcf.head[format.types,1],"Description=")  
  format.types<-extract.value.from.format(vcf.head[format.types,1],"ID=")
  
  
  column.labels[match("#CHROM",column.labels)]<-"chr"
  
  out.list<-list(column.labels,skip.lines,num.vars,info.types,info.class,info.description,format.types,format.class,format.description)
  names(out.list)<-c("column.labels","skip.lines","num.vars","info.types","info.class","info.description","format.types","format.class","format.description")
  out.list
}

extract.value.from.format<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
  match.string.general<-paste(match.string,"[a-zA-Z0-9_\\ \\.\\-]*",sep="") #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
  #match.string.general<-paste(match.string,"[\\w]*",sep="")
  match<- regexec(match.string.general,info)
  start<-unlist(match)
  length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
  end<-start+length-1 # minus one to get last position
  start<-start+nchar(as.character(match.string))
  substr(info,start=start,stop=end)
}
info.parse<-function(info2,format.types,format.string,complex.format=FALSE){
  #print(complex.format)
  if(vcf.type=="v3") {
    no.reads<-grepl("./.",info2)  # info2=="./." # no reads so no information sum(no.reads) ## "./.:.:1" oct GATK recal versions 2.7-2-g6bda569
    info2[no.reads]<-"NA/NA"
    info2<-as.matrix(info2)
    colnames(info2)<-c("GT")
    
  }else{
    
    #info2<-apply(info2,1,process.format) # original way
    
    ####multicore way
    if(complex.format){  ## if complex format the the format column is not always the same
      
      ## info2<-cbind(format.string,info2)
      ## info2<-apply(info2,1,list)
      ## key.indel<-1:dim(info2)[1]
      ## names(info2)<-key.indel
      ## info2<-mclapply(info2,process.format.list,mc.cores=1)
      ## posns<-match(key.indel,names(info2))
      ## info2<-info2[posns]
      
      ############use same trick used to process info in process.GATK.indels
      
      format.string<-strsplit(format.string,split=":") # unique(unlist(lapply(info,length)))
      format.length<-unlist(lapply(format.string,length))
      
      
      no.reads<-(grepl("./.",info2,fixed=TRUE) | (info2=="."))  #  (info2=="./.") | (info2==".") # no reads so no information sum(no.reads) changes Nov 4 2013 seen: "./.:.:1"
      info2[no.reads]<-unlist(lapply(format.length[no.reads],function(x){ paste(rep(NA,times=x),collapse=":")}))
      
      
      # WARNING info2 The format string can have MORE ITEMS than in sample string in merged VCF files
      ### fix this below make sure have the same length
      format.num<-strsplit(info2,split=":")
      format.num.length<-unlist(lapply(format.num,length))
      
      ## sum(format.num.length !=format.length) #
      if(sum(format.num.length >format.length)){print("ERROR ERROR format and genotype mismatch")}
      to.fix<-format.num.length !=format.length ## these have all got more formats than lengths
      if(sum(to.fix)>0){
        to.strip<-format.length[to.fix] -format.num.length[to.fix]
        to.strip.types<-tapply(to.strip,to.strip,length)
        for(is in 1:length(to.strip.types)){
          
          the.strip.type<-as.integer(names(to.strip.types)[is])
          to.fix.subset<-to.strip==the.strip.type
          new.format.string<-lapply(format.string[to.fix][to.fix.subset],function(x,the.strip.type){x[1:(length(x)-1)]})
          format.string[to.fix][to.fix.subset]<-new.format.string
          format.length[to.fix][to.fix.subset]<-format.num.length[to.fix][to.fix.subset]
        } #strip.types many need to take off more than one element
      } # to.fix required
      
      ########################
      
      
      format.labels<-unlist(format.string)
      format.num<-unlist(strsplit(info2,split=":"))
      
      format.index<-rep(1:length(info2),times=format.length)
      info2<-matrix(data=NA,nrow=length(info2),ncol=length(format.types))
      colnames(info2)<-format.types
      
      # i.info.col<-5
      for(i.info.col in 1:length(format.types)){
        the.col<-format.types[i.info.col]
        posns.in.flatten.info<-grep(paste("^",the.col,"$",sep=""),format.labels,fixed=FALSE) ## can have AC or AC1 etc correct make 01/21/2013
        posns<-format.index[posns.in.flatten.info]
        info2[posns,the.col]<-format.num[posns.in.flatten.info]
      }
      
      
    }else{ # not complex format
      no.reads<-(grepl("./.",info2,fixed=TRUE) | (info2=="."))  #   (info2=="./.") | (info2==".") # no reads so no information sum(no.reads) changes Nov 4 2013 seen: "./.:.:1"
      info2[no.reads]<-gsub(", ",":",toString(rep("NA",times=length(format.types)))) #  this will be the longest required 
      info2<-strsplit(info2,split=":")  # fastest  3.890   0.010   3.894
      num.info2<-length(format.types)
      len.info2<-length(info2)
      info2<-unlist(info2)
      dim(info2)<-c(num.info2,len.info2)  
      info2<-t(info2)
      colnames(info2)<-format.types
    } ## not complex format
    
  } ## not v3 VCF
  info2
} ### end sub info.process
# sum(info2.ori[,1] != info2[,1])

## x<-as.matrix(allele.depths.group)[,1]
allele.summary.individuals<-function(x){  ### work on COLUMNS
  #  print(x)
  missing<-is.na(x) | x=="NA,NA" | x=="NA" | !grepl(",",x) # not an NA must have a comma
  het<-rep(NA,times=length(x))
  x<-x[!missing]
  x<-as.numeric(unlist(strsplit(x,split=",")))
  dim(x)<-c(2,sum(!missing)) 
  ref.hom<-x[2,]==0
  alt.hom<-x[1,]==0
  x<-x[2,]/(x[2,]+x[1,])
  
  het[!missing]<-x
  het[!missing][ref.hom]<-0
  het[!missing][alt.hom]<-1
  het
}