library("GLAD")
library(zoo)
library(Rsge)

#################
## Constants
#################

## Compute log2 R ratio
## Add pseudo-signal to prevent -Inf and NaNs
DOUBLE_LOSS_THRESHOLD = -3
GAIN_THRESHOLD = 1
NORMAL_AB_LOWER_BAF_THRESHOLD = 0.33
NORMAL_AB_UPPER_BAF_THRESHOLD = 0.66
PSEUDO = 0.000001
CORFACTOR = 0.53

##MIN_REQUIRED_HET_DENSITY = 1/20000 ## One het SNP every 10 kb
MIN_REQUIRED_HET_SNPS = 15 ## At least 5 SNPs confirm an LOH event

# Method ID 1
# Take all samples from one individual
# Take all breakpoints of events from one indvidual
# Partition the genome based on ALL breakpoints
# This procedure is somewhat immune to having any kind of events
# even if only single, double copy loss is passed, it will partition on that,
# if AAB,BBA,AABB,A,- events are passed, it will partition on that

##############
## Functions
##############

match.conditions = function(x, low.bounds, up.bounds){
  res = NULL
  for(i in 1:length(x)){
    if(length(low.bounds)==length(up.bounds)){
      res=c(res,c((x[i] >= low.bounds) & (x[i] <= up.bounds)))
    } else {
      stop("Vector lengths must match")
    }
  }
  as.numeric(res)
}


blank.homs.log2r = function(daf.data,log2r.data){
  indxs = which(!is.na(daf.data))
  ## Check for chromosome Y, which has all its SNPs having DAF of 0.0
  if(length(daf.data)>=length(indxs)+2){
    log2r.data[indxs] = NA ## Blank all heterozygous SNPs
  }
  log2r.data
}

blank.hets.log2r = function(daf.data,log2r.data){
  indxs = which(is.na(daf.data))    
  ## Check for chromosome Y, which has all its SNPs having DAF of 0.0
  if(length(daf.data)>=length(indxs)+2){
    log2r.data[indxs] = NA ## Blank all homozygous SNPs
  }
  log2r.data
}

## Haplotype phasing function

phase.region = function(x,vec1){
  ##THRES = 0.8
  flipper.1 = c("A","AA","AAA","AAB")
  flipper.2 = c("B","BB","BBB","BBA")
  start.pos = as.numeric(x["start_snp_id"])
  stop.pos = as.numeric(x["stop_snp_id"])
  res1 = vec1[start.pos:stop.pos,"baf"]>CORFACTOR ## CORFACTOR because of dye-bias
  ##res2 = vec2[start.pos:stop.pos,"baf"]>0.5
  res3 = ab.snps[start.pos:stop.pos]%in%c(0,1)
  phased = ab.snps[start.pos:stop.pos]
  if(sum(res3)==0){
    x["state"]     # If there are no informative snps, return the same state
  } else {
    ##perc.conc = sum(res1[res3]==res2[res3])/sum(res3)   
    if(sum(phased[res3]==res1[res3]) > sum(phased[res3]!=res1[res3])){
      flipper.2[which(flipper.1==x["state"])] # Majority of matches, call it the "B" haplotypes
    } else {
      x["state"]    # Majority of opposites, call it the "A" haplotype
    }
  }
}


###############
## Main Program
###############

# Connect to mySQL database
library(RMySQL)
passwd = as.character(read.table("/etc/keys/.mysqlpass",as.is=T)[1,1])
con = dbConnect(MySQL(), user="root", password=passwd, dbname="sbep", host="localhost")

## Load snps
snps_all = dbGetQuery(con, "select chromosome, position, snp_id from snps order by snp_id")
## Filter and replace X, Y with 23, 24
snps_all = snps_all[which(snps_all$chromosome %in% c(c(1:22),"X","Y")),]
snps_all[which(snps_all$chromosome=="X"),"chromosome"] = 23
snps_all[which(snps_all$chromosome=="Y"),"chromosome"] = 24

snps_annotation = dbGetQuery(con, "select * from snps_annotation")
## Fix X=23 Y=24
snps_annotation[23,1]=23
snps_annotation[24,1]=24


## #######################
## For features calling
## Get snps
snps = dbGetQuery(con, "select chromosome, position, snp_id from snps order by chromosome, position")
## Make bins
## Replace chromosome X and Y by 23 and 24
snps[snps[,"chromosome"]=="X","chromosome"] = 23
snps[snps[,"chromosome"]=="Y","chromosome"] = 24
## Remove all extra snps, those unmapped on chromosomes (position is 0)
## and those mapped on MT, UNKNOWN, XY chromosomes.
snps = snps[snps[,"chromosome"] %in% c(1:24),]
snps = snps[snps[,"position"]>0,]

## End of features calling
## ########################

## We did 672 already
patients = dbGetQuery(con, "select distinct(patient_id) from acs_paired_samples where patient_id!=672")

## ########################
## Start of GLAD processing
## ########################
for(pind in 1:nrow(patients)){
  
  patient_id = patients[pind,1]  
  samples = dbGetQuery(con, paste("select * from acs_paired_samples",
    "where patient_id=",patient_id))

  ## Previous reference id will be used to avoid loading reference samples every time per patient
  prev_reference_id = 0
  
  for(i in 1:nrow(samples)){

    sample_id = samples[i,"sample_id"]
    reference_id = samples[i,"reference_id"]
    print(paste("Processing pair ",i," ",sample_id," ",reference_id))

    ## Check if the sample is slotted into MEMORY
    check = dbGetQuery(con,paste("select * from data_memory where sample_id=",sample_id,"limit 1"))
    if(length(check)>0){
      sample_data=dbGetQuery(con,paste("select r, baf, snp_id from data_memory where sample_id=",
        sample_id,"order by snp_id"))
    } else {
      sample_data=dbGetQuery(con,paste("select r, baf, snp_id from data where sample_id=",
        sample_id,"order by snp_id"))
    }
    
    
    if(prev_reference_id!=reference_id){
      check = dbGetQuery(con,paste("select * from data_memory where sample_id=",reference_id,"limit 1"))
      if(length(check)>0){
        reference_data=dbGetQuery(con,paste("select r, baf, snp_id from data_memory where sample_id=",
          reference_id,"order by snp_id"))
      } else {
        reference_data=dbGetQuery(con,paste("select r, baf, snp_id from data where sample_id=",
          reference_id,"order by snp_id"))
      }
    }
    prev_reference_id = reference_id
    
    ## Get the bad probes from the reference and remove them from both reference and sample
    snp_ids.to.filter=dbGetQuery(con,paste("select snp_id from sample_filtered_snps where sample_id=",
      reference_id,"order by snp_id"))
    
    ## Get rid of some 
    sample_data = sample_data[which(!sample_data[,"snp_id"]%in%snp_ids.to.filter[,1]),]
    reference_data = reference_data[which(!reference_data[,"snp_id"]%in%snp_ids.to.filter[,1]),]
    
    ## ###############
    ## Pre-processing
    ## ###############
    
    ## 1. Remove all probes that have low signal in the reference sample for the patient
    ## These are CNV probes that the patient lacks, or genomic microregions that are missing
    snps = snps_all[which(!snps_all[,"snp_id"]%in%snp_ids.to.filter[,1]),]
    
    
    ## 2. Take the Log2 (R sample / R reference) and
    ##    take the difference between B allele frequencies of sample and reference
    ##    or "DAF" - difference in allelic frequency
    ## Shift reference BAF data to get it centered at 0.5
    ## Normalize datasets, center at 0.55

    if(prev_reference_id!=reference_id){
      reference_data[reference_data[,"baf"]<=CORFACTOR,"baf"]=reference_data[reference_data[,"baf"]<=CORFACTOR,"baf"]*0.5/CORFACTOR
      reference_data[reference_data[,"baf"]>CORFACTOR,"baf"]=0.5+0.5*(reference_data[reference_data[,"baf"]>CORFACTOR,"baf"]-CORFACTOR)/(1-CORFACTOR)
    }
    sample_data[sample_data[,"baf"]<=CORFACTOR,"baf"]=sample_data[sample_data[,"baf"]<=CORFACTOR,"baf"]*0.5/CORFACTOR
    sample_data[sample_data[,"baf"]>CORFACTOR,"baf"]=0.5+0.5*(sample_data[sample_data[,"baf"]>CORFACTOR,"baf"]-CORFACTOR)/(1-CORFACTOR)
    
      
    data = cbind(snps$snp_id, snps,
      log2((sample_data[,"r"]+PSEUDO)/(reference_data[,"r"]+PSEUDO)),
      abs(sample_data[,"baf"]-0.5)*2)

    colnames(data) = c("snp_id", "Chromosome","PosBase","PosOrder","LogRatio","DAF")
    
    ## 3. Threshold the Log2R for double loss
    data[which(data$LogRatio < DOUBLE_LOSS_THRESHOLD), "LogRatio"] = DOUBLE_LOSS_THRESHOLD
    
    ## 4. Threshold, so that all homozygous SNPs in the reference are blanked, or have DAF = 0
    data[which(reference_data[,"baf"]<NORMAL_AB_LOWER_BAF_THRESHOLD | reference_data[,"baf"]>NORMAL_AB_UPPER_BAF_THRESHOLD), "DAF"]=NA
    
    ## Combined log2r and daf
    ## If log2r exceeds threshold for gain, add DAF to Log2R, otherwise, subtract DAF from Log2R
    ## This will boost signal  
    
    hets.logr.profile = data[,c("snp_id", "Chromosome","PosBase","PosOrder","LogRatio")]
    homs.logr.profile = data[,c("snp_id", "Chromosome","PosBase","PosOrder","LogRatio")]
    hets.daf.profile = data[,c("snp_id", "Chromosome","PosBase","PosOrder","DAF")]
    colnames(hets.daf.profile) = c("snp_id", "Chromosome","PosBase","PosOrder","LogRatio")
    
    ## Compute a n-SNP moving average to detect stable copy gains
    for(chr in 1:24){
      c.start = snps_annotation[snps_annotation[,"chromosome"]==chr,"start"]
      c.p = snps_annotation[snps_annotation[,"chromosome"]==chr,"p"]
      c.q = snps_annotation[snps_annotation[,"chromosome"]==chr,"q"]
      c.stop = snps_annotation[snps_annotation[,"chromosome"]==chr,"stop"]
      if(c.start!=c.p){
        indxs = which(data$Chromosome==chr & data$snp_id>=c.start & data$snp_id<=c.p)
        if(chr%in%c(23,24)){
          hets.logr.profile[indxs,"LogRatio"] = data[indxs,"LogRatio"]
          homs.logr.profile[indxs,"LogRatio"] = data[indxs,"LogRatio"]
          hets.daf.profile[indxs,"LogRatio"] = rep(NA,length(indxs))
        } else {
          hets.logr.profile[indxs,"LogRatio"] = blank.hets.log2r(data[indxs,"DAF"],data[indxs,"LogRatio"])
          homs.logr.profile[indxs,"LogRatio"] = blank.homs.log2r(data[indxs,"DAF"],data[indxs,"LogRatio"])
        }
      }
      indxs = which(data$Chromosome==chr & data$snp_id>=c.q & data$snp_id<=c.stop)
      if(chr%in%c(23,24)){
        hets.logr.profile[indxs,"LogRatio"] = data[indxs,"LogRatio"]
        homs.logr.profile[indxs,"LogRatio"] = data[indxs,"LogRatio"]
        hets.daf.profile[indxs,"LogRatio"] = rep(NA,length(indxs))
      } else {
        hets.logr.profile[indxs,"LogRatio"] = blank.hets.log2r(data[indxs,"DAF"],data[indxs,"LogRatio"])
        homs.logr.profile[indxs,"LogRatio"] = blank.homs.log2r(data[indxs,"DAF"],data[indxs,"LogRatio"])
      }

      ## Affix centromeres to have Log2R of 0 and DAF of 0
      if(c.p!=c.q-1){
        indxs = which(data$Chromosome==chr & data$snp_id>c.p & data$snp_id<c.q)
        hets.logr.profile[indxs,"LogRatio"] = NA
        homs.logr.profile[indxs,"LogRatio"] = NA
        hets.daf.profile[indxs,"LogRatio"] = NA
      }
      if(c.stop<max(data[data$Chromosome==chr,"snp_id"])){
        indxs = which(data$Chromosome==chr & data$snp_id>c.stop)
        hets.logr.profile[indxs,"LogRatio"] = NA
        homs.logr.profile[indxs,"LogRatio"] = NA
        hets.daf.profile[indxs,"LogRatio"] = NA
      }
    }

    ## Post-threshold the Log2r for double loss
    indxs = which(hets.logr.profile$LogRatio < DOUBLE_LOSS_THRESHOLD)
    if(any(indxs)){
      hets.logr.profile[indxs, "LogRatio"] = DOUBLE_LOSS_THRESHOLD
    }
    indxs = which(homs.logr.profile$LogRatio < DOUBLE_LOSS_THRESHOLD)
    if(any(indxs)){      
      homs.logr.profile[indxs, "LogRatio"] = DOUBLE_LOSS_THRESHOLD
    }

        
    ## Run GLAD on the 3 profiles
    ## The FDR parameter 0<q<0.5 controls the FDR of breakpoints in the data.
    ## Low values of q will reduce the false-positives at the possible cost of
    ## increasing the false-negatives, and vice versa.
    ## Warning - GLAD changes the values in the PosOrder column
    ## resulting in loosing the association between PosOrder and snp_id
    

    ## hets.logr.glad = daglad(cg, smoothfunc="lawsglad", lkern="Exponential", model="Gaussian",
    ##           qlambda=0.999,  bandwidth=1, base=FALSE, round=1.5,
    ##           lambdabreak=8, lambdaclusterGen=40, param=c(d=6), alpha=0.001, msize=5,
    ##           method="centroid", nmin=1, nmax=8,
    ##           amplicon=1, deletion=-1, deltaN=0.1,  forceGL=c(-0.1,0.1), nbsigma=3,
    ##           MinBkpWeight=0.35, CheckBkpPos=TRUE)
    haar.start = 3
    haar.end = 9
    haar.fdr = 0.00000001
        
    cg = as.profileCGH(hets.logr.profile)    
    hets.logr.glad = daglad(cg, mediancenter=F,  genomestep=F, smoothfunc="haarseg", 
      haarStartLevel=haar.start, haarEndLevel=haar.end, breaksFdrQ=haar.fdr,
      verbose=F, OnlySmoothing=T )
    ## PosOrder changed after GLAD, but snp_id stayed
    if(!is.null(hets.logr.glad$profileValuesNA)){
      n = nrow(hets.logr.glad$profileValuesNA)
      temp = cbind(hets.logr.glad$profileValuesNA,rep(NA,n),rep(NA,n),rep(NA,n),rep(0,n))
      colnames(temp) = colnames(hets.logr.glad$profileValues)
      hets.logr.glad = rbind(hets.logr.glad$profileValues, temp)
      hets.logr.glad = hets.logr.glad[order(hets.logr.glad$snp_id),]
    } else {
      hets.logr.glad = hets.logr.glad$profileValues
    }

    cg = as.profileCGH(homs.logr.profile)
    homs.logr.glad = daglad(cg, mediancenter=F,  genomestep=F, smoothfunc="haarseg", 
      haarStartLevel=haar.start, haarEndLevel=haar.end, breaksFdrQ=haar.fdr,
      verbose=F, OnlySmoothing=T )
    ## PosOrder changed after GLAD, but snp_id stayed
    if(!is.null(homs.logr.glad$profileValuesNA)){
      n = nrow(homs.logr.glad$profileValuesNA)
      temp = cbind(homs.logr.glad$profileValuesNA,rep(NA,n),rep(NA,n),rep(NA,n),rep(0,n))
      colnames(temp) = colnames(homs.logr.glad$profileValues)
      homs.logr.glad = rbind(homs.logr.glad$profileValues, temp)
      homs.logr.glad = homs.logr.glad[order(homs.logr.glad$snp_id),]
    } else {
      homs.logr.glad = homs.logr.glad$profileValues
    }
    
    cg = as.profileCGH(hets.daf.profile)
    hets.daf.glad = daglad(cg, mediancenter=F,  genomestep=F, smoothfunc="haarseg", 
      haarStartLevel=haar.start, haarEndLevel=haar.end, breaksFdrQ=haar.fdr,      
      verbose=F, OnlySmoothing=T )
    ## PosOrder changed after GLAD, but snp_id stayed
    if(!is.null(hets.daf.glad$profileValuesNA)){
      n = nrow(hets.daf.glad$profileValuesNA)
      temp = cbind(hets.daf.glad$profileValuesNA,rep(NA,n),rep(NA,n),rep(NA,n),rep(0,n))
      colnames(temp) = colnames(hets.daf.glad$profileValues)
      hets.daf.glad = rbind(hets.daf.glad$profileValues, temp)
      hets.daf.glad = hets.daf.glad[order(hets.daf.glad$snp_id),]
    } else {
      hets.daf.glad = hets.daf.glad$profileValues
    }
    rm(temp)
    gc()


    ## Take all break points post GLAD smoothing
    ## Make sure breakpoints is a vector of snp_ids
        
    ## Implant the start stop of p, q arms of chromosomes
    ## Go break point by break point to create features
    ## Features overlap by 1 SNP
    
    ## In GLAD, if the cnai spans snps N..M
    ## the breakpoint always are set to snps N-1 and M
    ## Therefore, the features have overlapping start and end SNPs
    ## e.g. chr1    1 .. N-1 N .. M M+1 ... Z-1 Z ...F F+1 ... E
    ## will have glad breakpoints N-1, M, Z-1, F
    ## which will get translated to features  1..N-1  N-1..M M..Z-1 Z-1..F F..E

    ## Combine the breakpoints from all 3 profiles
    breakpoints = hets.logr.glad[hets.logr.glad$Breakpoints==1 |
                                 homs.logr.glad$Breakpoints==1 |
                                 hets.daf.glad$Breakpoints==1,"snp_id"]
    
    all.bins = NULL
    for(chr in 1:24){
      bins.chr = NULL
      ## run p arm, should handle acrocentric chromosomes   
      bins.p.arm = cbind(c(snps_annotation[chr,"start"],
        breakpoints[breakpoints>snps_annotation[chr,"start"] & breakpoints<snps_annotation[chr,"p"]]),
        c(breakpoints[breakpoints>snps_annotation[chr,"start"] & breakpoints<snps_annotation[chr,"p"]],
          snps_annotation[chr,"p"]))
      ## run q arm
      bins.q.arm = cbind(c(snps_annotation[chr,"q"],
        breakpoints[breakpoints>snps_annotation[chr,"q"] & breakpoints<snps_annotation[chr,"stop"]]),
        c(breakpoints[breakpoints>snps_annotation[chr,"q"] & breakpoints<snps_annotation[chr,"stop"]],
          snps_annotation[chr,"stop"]))
      bins.chr = rbind(bins.p.arm,bins.q.arm)
      if(nrow(bins.chr)>0){
        bins.chr = cbind(rep(chr, nrow(bins.chr)),bins.chr)
      }
      all.bins = rbind(all.bins, bins.chr)
    }
    colnames(all.bins) = c("Chromosome","Start","Stop")
    
    get.first.non.na.element = function(arr){
      for(arrind in 1:length(arr)){
        if(!is.na(arr[arrind])){
          return(arr[arrind])
        }
      }
      return(NA)
    }
    
    ## Molecular states table, hardcoded
    ## State, upper bound on Log2r, lower bound on Log2R, upper bound 
    sample_features = NULL
    for(j in 1:nrow(all.bins)){
      if(j %% 100 == 0 ) {
        print(paste(j,"features created out of ",nrow(all.bins)))
      }
      ## A segment spans N..M
      ## GLAD's breakpoints are set to N-1 and M
      ##
      seg.start = which(hets.logr.glad$snp_id == all.bins[j,"Start"]+1)
      seg.stop  = which(hets.logr.glad$snp_id == all.bins[j,"Stop"])
      if(length(seg.start)==0){
        seg.start = max(which(hets.logr.glad$snp_id < all.bins[j,"Start"]+1))
      }
      if(length(seg.stop)==0){
        seg.stop = max(which(hets.logr.glad$snp_id < all.bins[j,"Stop"]))
      }
      het.logr = get.first.non.na.element(hets.logr.glad[seg.start:seg.stop,"Smoothing"])
      hom.logr = get.first.non.na.element(homs.logr.glad[seg.start:seg.stop,"Smoothing"])
      mbaf = get.first.non.na.element(hets.daf.glad[seg.start:seg.stop,"Smoothing"])
      ## Count the number of heterozygous snps support for this feature
      num.het.snps = sum(!is.na(hets.logr.glad[seg.start:seg.stop,"LogRatio"]))
      ## Size of this feature in base pairs
      feature.size = hets.logr.glad[seg.stop,"PosBase"]-hets.logr.glad[seg.start,"PosBase"]+1

      ## Determine the molecular state of this feature  
      sample_features = rbind(sample_features,cbind(
        NA,
        sample_id,
        reference_id,
        all.bins[j,"Chromosome"],
        hets.logr.glad[seg.start,"PosBase"],
        hets.logr.glad[seg.stop,"PosBase"],
        hets.logr.glad[seg.start,"snp_id"],
        hets.logr.glad[seg.stop,"snp_id"],
        round(het.logr,digits=4),
        round(hom.logr,digits=4),
        round(mbaf,digits=4),
        num.het.snps,
        feature.size))
     }

    sample_features = as.data.frame(sample_features)
    colnames(sample_features)=c("feature_id","sample_id","reference_id","chromosome",
              "start_position","stop_position", "start_snp_id", "stop_snp_id",
              "het_logr", "hom_logr", "mbaf","num_het_snps","size")
    row.names(sample_features) = NULL
    ##sample_features[sample_features[,"mean_daf"]=="NaN","mbaf"] = NA
    ## Save to sbepdb
    rs = dbGetQuery(con, paste("delete from sample_features where sample_id=",sample_id,"and reference_id=",reference_id))
    dbWriteTable(con, "sample_features", sample_features, row.names=F, append=TRUE)
            
  }  ## End of looping the samples per patient
}
## ######################
## End of GLAD processing
## ######################



## #################################
## Start of patient feature state calling
## #################################

for(pind in 1:nrow(patients)){
##for(pind in 12:12){
  
  patient_id = patients[pind,1]
  print(paste("Patient",patient_id))
  samples = dbGetQuery(con, paste("select * from acs_paired_samples",
    "where patient_id=",patient_id))

  ## #########################
  ## Start of state assignment
  ## #########################
##  for(kkk in 1:1){  
  
  ## Get state calling thresholds
  sthres = dbGetQuery(con, "select * from state_calling_thresholds")

  for(j in 1:nrow(samples)){
    print(paste("Molecular state assignment sample",j,"done"))
    sample_id = samples[j,"sample_id"]
    reference_id = samples[j,"reference_id"]
    ## Assign molecular states to features
    ## Cartesian product of states and features:
    features = dbGetQuery(con, paste("select * from sample_features where sample_id=",
      sample_id,"and reference_id=",reference_id)) 
    sf = expand.grid(sthres[,"state_id"],features[,"feature_id"])
    fs.het.logr = match.conditions(features[,"het_logr"],sthres[,"het_logr_lower"],sthres[,"het_logr_upper"])
    fs.hom.logr = match.conditions(features[,"hom_logr"],sthres[,"hom_logr_lower"],sthres[,"hom_logr_upper"])
    fs.mbaf     = match.conditions(features[,"mbaf"],    sthres[,"mbaf_lower"],sthres[,"mbaf_upper"])
    sid = features[match(sf[,2],features$feature_id),"sample_id"]
    rid = features[match(sf[,2],features$feature_id),"reference_id"]
    fsp = cbind(sid,rid,sf[,2],sf[,1],sthres[sf[,1],"state"],fs.het.logr,fs.hom.logr,fs.mbaf)
    colnames(fsp) = c("sample_id","reference_id","feature_id","state_id","state",
              "het_logr_prob","hom_logr_prob","mbaf_prob")
    fsp = as.data.frame(fsp)    
    rs = dbGetQuery(con, paste("delete from feature_state_probabilities where sample_id=",sample_id,"and reference_id=",reference_id))
    dbWriteTable(con, "feature_state_probabilities", fsp, row.names=F, append=T)
  }
  
  ## ##################################
  ## Start of creating patient features
  ## ##################################
  breakpoints = NULL
  for(j in 1:nrow(samples)){
    print(paste("Creating features sample",j,"done"))
    sample_id = samples[j,"sample_id"]
    reference_id = samples[j,"reference_id"]
    d=dbGetQuery(con,paste("select * from feature_state_probabilities f inner join sample_features s
                          on s.feature_id=f.feature_id
                          where s.sample_id=",sample_id," and s.reference_id=",reference_id))
    ## state_AA = d[d$state_id==3 & d$het_logr_prob==1 & d$hom_logr_prob==1 & d$mbaf_prob==1 &
    ##              d$num_het_snps > MIN_REQUIRED_HET_SNPS &
    ##              d$num_het_snps/d$size > MIN_REQUIRED_HET_DENSITY, "feature_id"]
        state_AA = d[d$state_id==3 & (d$het_logr_prob==1 | d$hom_logr_prob==1) & d$mbaf_prob==1 &
                 d$num_het_snps > MIN_REQUIRED_HET_SNPS, "feature_id"]

    state_gain = d[d$state_id%in%c(4,5,7,8) & ( d$het_logr_prob==1 | d$hom_logr_prob==1 )& d$mbaf_prob==1,"feature_id"]
    state_A = d[d$state_id==2 & (d$hom_logr_prob==1 | d$het_logr_prob==1) & d$mbaf_prob==1, "feature_id"]
    state_0 = d[d$state_id==1 & d$het_logr_prob==1 & d$hom_logr_prob==1 & d$mbaf_prob==1, "feature_id"]
    state_0_mixed = d[d$state_id==2 & d$het_logr < -1 & d$hom_logr < -1, "feature_id"]
    
    orig.feats = sort(unique(d$feature_id))
    ## The default state is AB for all features
    feats = rep("AB",length(orig.feats)) 
    feats[match(state_AA,orig.feats)] = "AA"
    feats[match(state_gain,orig.feats)] = "AAB"
    feats[match(state_A,orig.feats)] = "A"
    feats[match(state_0,orig.feats)] = "0"
    feats[match(state_0_mixed,orig.feats)] = "0"
    
    ## Gather all breakpoints
    breakpoints = sort(unique(c(breakpoints,d[,"start_snp_id"],d[,"stop_snp_id"])))
  }
  
  ## Mind the snps_annotation!
  all.bins = NULL
  for(chr in 1:24){
    bins.chr = NULL
    ## run p arm, should handle acrocentric chromosomes
    if(snps_annotation[chr,"start"]!=snps_annotation[chr,"p"]){
      bins.p.arm = cbind(c(snps_annotation[chr,"start"],
        breakpoints[breakpoints>snps_annotation[chr,"start"] & breakpoints<snps_annotation[chr,"p"]]),
        c(breakpoints[breakpoints>snps_annotation[chr,"start"] & breakpoints<snps_annotation[chr,"p"]],
          snps_annotation[chr,"p"]))
    } else {
      bins.p.arm = NULL
    }
    ## run q arm
    bins.q.arm = cbind(c(snps_annotation[chr,"q"],
      breakpoints[breakpoints>snps_annotation[chr,"q"] & breakpoints<snps_annotation[chr,"stop"]]),
      c(breakpoints[breakpoints>snps_annotation[chr,"q"] & breakpoints<snps_annotation[chr,"stop"]],
        snps_annotation[chr,"stop"]))
    bins.chr = rbind(bins.p.arm,bins.q.arm)
    if(nrow(bins.chr)>0){
      bins.chr = cbind(rep(chr, nrow(bins.chr)),bins.chr)
    }
    all.bins = rbind(all.bins, bins.chr)
  }
  colnames(all.bins) = c("chromosome","start_snp_id","stop_snp_id")
  all.bins = as.data.frame(all.bins)
  
  for(j in 1:nrow(samples)){
    sample_id = samples[j,"sample_id"]
    reference_id = samples[j,"reference_id"]
    print(paste("Creating features 2 sample ",j,"done"))
    d=dbGetQuery(con,paste("select * from feature_state_probabilities f inner join sample_features s
                          on s.feature_id=f.feature_id
                          where s.sample_id=",sample_id," and s.reference_id=",reference_id))
    ## state_AA = d[d$state_id==3 & d$het_logr_prob==1 & d$hom_logr_prob==1 & d$mbaf_prob==1 &
    ##              d$num_het_snps > MIN_REQUIRED_HET_SNPS &
    ##              d$num_het_snps/d$size > MIN_REQUIRED_HET_DENSITY, "feature_id"]
    ## state_AA = d[d$state_id==3 & d$het_logr_prob==1 & d$hom_logr_prob==1 & d$mbaf_prob==1 &
    ##   d$num_het_snps > MIN_REQUIRED_HET_SNPS, "feature_id"]

    ## state_gain = d[d$state_id%in%c(4,5,7,8) & d$het_logr_prob==1 & d$hom_logr_prob==1 & d$mbaf_prob==1,"feature_id"]
    ## state_A = d[d$state_id==2 & d$hom_logr_prob==1 & d$het_logr_prob==1 & d$mbaf_prob==1, "feature_id"]
    state_AA = d[d$state_id==3 & (d$het_logr_prob==1 | d$hom_logr_prob==1) & d$mbaf_prob==1 &
      d$num_het_snps > MIN_REQUIRED_HET_SNPS, "feature_id"]

    state_gain = d[d$state_id%in%c(4,5,7,8) & ( d$het_logr_prob==1 | d$hom_logr_prob==1 )& d$mbaf_prob==1,"feature_id"]
    state_A = d[d$state_id==2 & (d$hom_logr_prob==1 | d$het_logr_prob==1) & d$mbaf_prob==1, "feature_id"]

    state_0 = d[d$state_id==1 & d$het_logr_prob==1 & d$hom_logr_prob==1 & d$mbaf_prob==1, "feature_id"]
    state_0_mixed = d[d$state_id==2 & d$het_logr < -1 & d$hom_logr < -1, "feature_id"]

    orig.feats = sort(unique(d$feature_id))
    feats = rep("AB",length(orig.feats))
    feats[match(state_AA,orig.feats)] = "AA"
    feats[match(state_gain,orig.feats)] = "AAB"
    feats[match(state_A,orig.feats)] = "A"
    feats[match(state_0,orig.feats)] = "0"
    feats[match(state_0_mixed,orig.feats)] = "0"


    ## Assign states to all bins
    fff = cbind(d[match(orig.feats,d$feature_id),c(1,2,12:16)],feats)
    
    state = NULL
    for(k in 1:nrow(all.bins)){
      state[k] = as.character(fff[fff$start_snp_id<=all.bins[k,"stop_snp_id"] &
             fff$stop_snp_id>=all.bins[k,"start_snp_id"],"feats"])[1]
    }
    n = length(state)
    final = as.data.frame(cbind(rep(NA,n),rep(sample_id,n),rep(reference_id,n),rep(patient_id,n),
      all.bins$chromosome,
      snps[match(all.bins$start_snp_id,snps$snp_id),"position"],
      snps[match(all.bins$stop_snp_id,snps$snp_id),"position"],
      all.bins$start_snp_id,all.bins$stop_snp_id,state))
    colnames(final) = c("feature_id","sample_id","reference_id","patient_id","chromosome",
              "start_position","stop_position","start_snp_id","stop_snp_id","state")
    ## delete previous analysis
    rs = dbSendQuery(con,paste("delete from patient_feature_states where sample_id=",sample_id))
    dbWriteTable(con, "patient_feature_states", final, row.names=F, append=TRUE)
  }
  
  ## ###################################
  ## End of feature creation per patient
  ## ###################################

  ## #############################
  ## Compress features per patient
  ## #############################

  ## Assuming that all features have same start and stop locations across samples
  for(j in 1:nrow(samples)){
    sample_id = samples[j,"sample_id"]
    reference_id = samples[j,"reference_id"]
    ds=dbGetQuery(con,paste("select * from patient_feature_states",
      "where sample_id=",sample_id," and reference_id=",reference_id))
    if(j==1){
      d = ds
    } else {
      d = cbind(d,ds[,"state"])
    }
  }
  
  ## Compress/compact/fuse/merge neighbor features that have the same molecular state
  ## across all samples
  cur_j = 1
  cur_chromosome = d[1,"chromosome"]
  cur_start_pos  = d[1,"start_position"]
  cur_stop_pos   = d[1,"stop_position"]
  cur_start_sid  = d[1,"start_snp_id"]
  cur_stop_sid   = d[1,"stop_snp_id"]
  new_d = NULL
  for(j in 2:nrow(d)){
    if(j %% 100 == 0 ) {
      print(paste(j,"features merged out of ",nrow(d)))
    }
    ## Col 10 is the mol. state of the first sample
    ## if(d[j,"start_snp_id"]==cur_stop_sid &&
    ##    d[j,"chromosome"]==cur_chromosome &&
    ##    sum(!d[cur_j,10:ncol(d)]==rep("AB",ncol(d)-9))==0 &&
    ##    sum(d[j,10:ncol(d)]==rep("0",ncol(d)-9))>0 &&
    ##    d[j,"stop_snp_id"]-d[j,"start_snp_id"]==2){
    ##   cur_stop_sid = d[j,"stop_snp_id"]
    ##   cur_stop_pos = d[j,"stop_position"]
    ##   ## Replace the singleton 0 with AB
    ##   d[j,10:ncol(d)] = rep("AB",ncol(d)-9)
    ## } else
    if(d[j,"start_snp_id"]==cur_stop_sid &&
       d[j,"chromosome"]==cur_chromosome &&
       sum(!d[j,10:ncol(d)]==d[cur_j,10:ncol(d)])==0){
      cur_stop_sid = d[j,"stop_snp_id"]
      cur_stop_pos = d[j,"stop_position"]
    } else {
     new_d = rbind(new_d, cbind(cur_chromosome,cur_start_pos,cur_stop_pos,cur_start_sid,cur_stop_sid,
       d[cur_j,10:ncol(d)]))
     cur_j = j
     cur_chromosome = d[j,"chromosome"]
     cur_start_pos  = d[j,"start_position"]
     cur_stop_pos   = d[j,"stop_position"]
     cur_start_sid  = d[j,"start_snp_id"]
     cur_stop_sid   = d[j,"stop_snp_id"]     
   }
  }
  
  
  ## Transform the compressed features to mysql format
  final = NULL
  for(j in 1:nrow(samples)){
    sample_id = samples[j,"sample_id"]
    reference_id = samples[j,"reference_id"]
    n = nrow(new_d)
    final = cbind(rep(NA,n),rep(sample_id,n),rep(reference_id,n),rep(patient_id,n),
      new_d[,c(1:5,5+j)],rep(NA,n))
    colnames(final) = c("feature_id","sample_id","reference_id","patient_id","chromosome",
              "start_position","stop_position","start_snp_id","stop_snp_id","state","state_phased")
    ## delete previous analysis
    rs = dbSendQuery(con,paste("delete from compressed_patient_feature_states where sample_id=",sample_id))
    dbWriteTable(con, "compressed_patient_feature_states", final, row.names=F, append=TRUE)   
  }
 
##}
  ## #######################################
  ## End of compressing features per patient
  ## #######################################

}  ## End of looping over patients
 

  ## ##########################
  ## Start of haplotype phasing
  ## ##########################

for(pind in 1:nrow(patients)){
  
  patient_id = patients[pind,1]
  print(paste("Patient",patient_id))
  samples = dbGetQuery(con, paste("select * from acs_paired_samples",
    "where patient_id=",patient_id))

  features = dbGetQuery(con, paste("select * from compressed_patient_feature_states",
    "where patient_id=",patient_id,"order by sample_id"))
    
  ## Get the phased data for the patient
  phased.snps = dbGetQuery(con, paste("select snp_id, b_allele_haplotype from phased_data",
    "where sample_id=",reference_id))
  ab.snps = rep(NA,nrow(snps))
  ab.snps[phased.snps[,"snp_id"]] = phased.snps[,"b_allele_haplotype"]
  ## Assign sample 1 the "AA" haplotype
  ## Thus, all AI calls have the "A" haplotype - AA, AAB, A, AAAB, etc.
  ## There will be no "B" haplotypes - BB, BBA, BBBA, B, etc.
  ##sample.ids = sort(unique(features[,"sample_id"]))
  
  for(j in 1:nrow(samples)){
    sample_id = samples[j,"sample_id"]
    print(paste("Phasing sample ",j))
    ## Check if the sample is slotted into MEMORY
    check = dbGetQuery(con,paste("select * from data_memory where sample_id=",sample_id,"limit 1"))
    if(length(check)>0){
      sample.baf=dbGetQuery(con,paste("select snp_id, baf from data_memory where sample_id=",
        sample_id,"order by snp_id"))
    } else {
      sample.baf=dbGetQuery(con,paste("select snp_id, baf from data where sample_id=",
        sample_id,"order by snp_id"))
    }
    ## sample.baf = dbGetQuery(con, paste("select snp_id, baf from data",
    ##   "where sample_id=",sample.ids[j]))
    sample.sub.features = features[features[,"sample_id"]==sample_id,]
    sample.sub.features[,"state_phased"] = sample.sub.features[,"state"] 
    AI.features = sample.sub.features[sample.sub.features[,"state"]%in%c("A","AA","AAA","AAB"),]
    AI.features[,"state_phased"] = apply(AI.features,1,phase.region,vec1=sample.baf)
    sample.sub.features[sample.sub.features[,"state"]%in%c("A","AA","AAA","AAB"),] = AI.features
    ## delete previous analysis
    rs = dbSendQuery(con,paste("delete from compressed_patient_feature_states where sample_id=",sample_id))
    dbWriteTable(con, "compressed_patient_feature_states", sample.sub.features, row.names=F, append=TRUE)   
  }
}
  ## ########################
  ## End of haplotype phasing
  ## ########################





    ## sample_features = NULL
    ## for(j in 1:nrow(all.bins)){
    ##   if(j %% 100 == 0 ) {
    ##     print(paste(j,"features created"))
    ##   }

    ##   ## Fix the real start and stop of a segment
    ##   seg.start = which(final$snp_id==all.bins[j,"Start"])
    ##   seg.stop  = which(final$snp_id==all.bins[j,"Stop"])
    ##   if(length(seg.start)==0){
    ##     ## if SNP was blanked, insert it to "final" data frame
    ##     ins.snp = as.numeric(c(snps_all[snps_all$snp_id==all.bins[j,"Start"],]))
    ##     ins.row = cbind(ins.snp[3],ins.snp[1],ins.snp[2],NA,NA,NA,NA,NA,NA,NA,NA)
    ##     colnames(ins.row) = colnames(final)
    ##     final = rbind(final[1:max(which(final$snp_id<all.bins[j,"Start"])),],
    ##       ins.row,
    ##       final[min(which(final$snp_id>all.bins[j,"Start"])):nrow(final),])
    ##     seg.start = which(final$snp_id==all.bins[j,"Start"])
    ##   }  
    ##   if(length(seg.stop)==0){
    ##     ## if SNP was blanked
    ##     ins.snp = as.numeric(c(snps_all[snps_all$snp_id==all.bins[j,"Stop"],]))
    ##     ins.row = cbind(ins.snp[3],ins.snp[1],ins.snp[2],NA,NA,NA,NA,NA,NA,NA,NA)
    ##     colnames(ins.row) = colnames(final)
    ##     final = rbind(final[1:max(which(final$snp_id<all.bins[j,"Stop"])),],
    ##       ins.row,
    ##       final[min(which(final$snp_id>all.bins[j,"Stop"])):nrow(final),])
    ##     seg.stop  = which(final$snp_id==all.bins[j,"Stop"])
    ##   }
    ##   seg.cnai = median(final[(seg.start+1):seg.stop,"Smoothing"],na.rm=T)
    ##   seg.log2r = median(data[(seg.start+1):seg.stop,"LogRatio"],na.rm=T)
    ##   ## Calculate log2r for region, excluding heterozygous SNPs, and including homozygous SNPs
    ##   ## only, which helps to determine spurious single copy loss regions in areas of
    ##   ## copy neutral LOH, such as on chr. 9
    ##   sub.d = data[(seg.start+1):seg.stop,]
    ##   if(length(sub.d)>0){
    ##     seg.hom.log2r = median(sub.d[sub.d[,"DAF"]==0,"LogRatio"],na.rm=T)
    ##   } else {
    ##     seg.hom.log2r = NA
    ##   }
    ##   seg.daf   = median(final[(seg.start+1):seg.stop,"DAF"],na.rm=T)
      
    ##   ## Determine the molecular state of this feature  
    ##   sample_features = rbind(sample_features,cbind(
    ##     "NA",
    ##     sample_id,
    ##     reference_id,
    ##     all.bins[j,"Chromosome"],
    ##     final[seg.start+1,"PosBase"],
    ##     final[seg.stop,"PosBase"],        
    ##     all.bins[j,"Start"],
    ##     all.bins[j,"Stop"],                                       
    ##     sprintf("%1.4f",seg.cnai),
    ##     sprintf("%1.4f",seg.log2r),
    ##     sprintf("%1.4f",seg.daf),
    ##     sprintf("%1.4f",seg.hom.log2r)))
    ## }
    
    ## sample_features = as.data.frame(sample_features)
    ## colnames(sample_features)=c("feature_id","sample_id","reference_id","chromosome","start_position",
    ##           "stop_position", "start_snp_id", "stop_snp_id",
    ##           "mean_cnai", "mean_log2r", "mean_daf", "mean_hom_log2r")
    ## row.names(sample_features) = NULL
    ## sample_features[sample_features[,"mean_daf"]=="NaN","mean_daf"] = NA
    ## ## Save to sbepdb
    ## rs = dbGetQuery(con, paste("delete from sample_features where sample_id=",sample_id,"and reference_id=",reference_id))
    ## dbWriteTable(con, "sample_features", sample_features, row.names=F, append=TRUE)








## GLAD sge function
## segment.region = function(cg){
##   glad.cg = daglad(cg, mediancenter=F,  genomestep=F, smoothfunc="haarseg", 
##                  haarStartLevel=4, haarEndLevel=5, breaksFdrQ = 0.001,
##                  deltaN=0.15, forceGL=c(-0.2,0.2), amplicon=1, deletion=-2, DelBkpInAmp=T,
##                  verbose=F)
##   ##glad.cg$profileValues
##   smoothed = rbind(glad.cg$profileValues,cg$profileValuesNA)
##   smoothed = smoothed[order(c(glad.cg$profileValues[,"PosOrder"],glad.cg$profileValuesNA[,"PosOrder"])),]
## }





## sub.start = 59000000
## sub.stop  = 61000000
## sub.chr = 3
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")
## sub.start = 75000000
## sub.stop  = 80000000
## sub.chr = 16
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")
## sub.start = 1
## sub.stop  = 50000000
## sub.chr = 18
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")



## par(mfrow=c(4,1))
## sub.start = 1
## sub.stop  = 50000000
## sub.chr = 9
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")
## sub.start = 59000000
## sub.stop  = 61000000
## sub.chr = 3
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")
## sub.start = 75000000
## sub.stop  = 80000000
## sub.chr = 16
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")
## sub.start = 1
## sub.stop  = 50000000
## sub.chr = 18
## sub.final = final[final$Chromosome==sub.chr&final[,"PosBase"]>=sub.start&final[,"PosBase"]<=sub.stop,]
## plot(sub.final$PosBase,sub.final$LogRatio,cex=0.5,ylim=c(-2,1))
## lines(sub.final$PosBase,sub.final$Smoothing,col="red")
## abline(h=0,col="blue")



      ## cnai[which(cnai$Chromosome==chr & cnai$snp_id>=c.q & cnai$snp_id<=c.stop),"DAF"] = impute.daf.linear(sub.cnai)
      ## ##cnai[which(cnai$Chromosome==chr & cnai$snp_id>=c.q & cnai$snp_id<=c.stop),"DAF"] = impute.daf(sub.cnai, WINDOW_DAF_IMPUTE)
      ## cnai[which(cnai$Chromosome==chr & cnai$snp_id>=c.q & cnai$snp_id<=c.stop),"LogRatio"] =
      ##   cnai[loss.snps.q,"LogRatio"] - DAF_SCALING_FACTOR*cnai[loss.snps.q,"DAF"]

      ## ## Note that we are losing WINDOW/2 snps by the moving filter
      ## if(c.p-c.start+1>WINDOW){             # consider acrocentrics
      ##   ## filt.p = rollmean(cnai[cnai$snp_id>=c.start & cnai$snp_id<=c.p,"LogRatio"], WINDOW,na.pad=T)
      ##   ##filt.p = rollmedian(cnai[cnai$snp_id>=c.start & cnai$snp_id<=c.p,"LogRatio"], WINDOW,na.pad=T)
      ##   filt.p = rollmedian(cnai[cnai$snp_id>=c.start & cnai$snp_id<=c.p,"LogRatio"],WINDOW,na.pad=T)
      ##   ##filt.log2r = rollmedian(cnai[cnai$snp_id>=c.start & cnai$snp_id<=c.p,"LogRatio"],SMOOTHING_LOG2R_WINDOW,na.pad=T)
      ##   ##abs.cnai = abs(filt.p)+DAF_SCALING_FACTOR*cnai[cnai$snp_id>=c.start & cnai$snp_id<=c.p,"DAF"]
      ##   ## Get the index of the snps
      ##   ## gain.snps.index = which(abs.cnai>=NORMAL_THRESHOLD & filt.p>=GAIN_THRESHOLD)
      ##   ## loss.snps.index = which(abs.cnai>=NORMAL_THRESHOLD & filt.p<=LOSS_THRESHOLD)
      ##   ##gain.snps.index = which(abs.cnai>=NORMAL_THRESHOLD & filt.p>=GAIN_THRESHOLD)
      ##   loss.snps.index = which(filt.p<=LOSS_THRESHOLD)
      ##   ##gain.snps.p = gain.snps.index + (which(cnai$snp_id==c.start)-1)
      ##   loss.snps.p = loss.snps.index + (which(cnai$snp_id==c.start)-1)
      ##   ## if(length(gain.snps.p)>0){
      ##   ##   ##          cnai[gain.snps.p,"LogRatio"] = cnai[gain.snps.p,"LogRatio"] + DAF_SCALING_FACTOR*cnai[gain.snps.p,"DAF"]
      ##   ##   cnai[gain.snps.p,"LogRatio"] = filt.log2r[gain.snps.index]+DAF_SCALING_FACTOR*cnai[gain.snps.p,"DAF"]          
      ##   ## }
      ##   if(length(loss.snps.p)>0){    
      ##     cnai[loss.snps.p,"LogRatio"] = cnai[loss.snps.p,"LogRatio"] - DAF_SCALING_FACTOR*cnai[loss.snps.p,"DAF"]
      ##     ##cnai[loss.snps.p,"LogRatio"] = filt.log2r[loss.snps.index] - DAF_SCALING_FACTOR*cnai[loss.snps.p,"DAF"]
      ##   }
      ##   ## If the abs.cnai does not pass the 0.2 threshold, then the SNP does not get the DAF boost
      ## }
      ## ##  filt.q = rollmean(cnai[cnai$snp_id>=c.q & cnai$snp_id<=c.stop,"LogRatio"], WINDOW,na.pad=T)
      ## filt.q = rollmedian(cnai[cnai$snp_id>=c.q & cnai$snp_id<=c.stop,"LogRatio"], WINDOW,na.pad=T)
      ## ##filt.log2r = rollmedian(cnai[cnai$snp_id>=c.q & cnai$snp_id<=c.stop,"LogRatio"],SMOOTHING_LOG2R_WINDOW,na.pad=T)
      ## ##abs.cnai = abs(filt.q)+DAF_SCALING_FACTOR*cnai[cnai$snp_id>=c.q & cnai$snp_id<=c.stop,"DAF"]
      ## ##gain.snps.index = which(abs.cnai>=NORMAL_THRESHOLD & filt.q>=GAIN_THRESHOLD)
      ## ##loss.snps.index = which(abs.cnai>=NORMAL_THRESHOLD & filt.q<=LOSS_THRESHOLD)
      ## loss.snps.index = which(filt.q<=LOSS_THRESHOLD)      
      ## ##gain.snps.q = gain.snps.index + (which(cnai$snp_id==c.q)-1)
      ## loss.snps.q = loss.snps.index + (which(cnai$snp_id==c.q)-1)  
      ## ## if(length(gain.snps.q)>0){
      ## ##   ## cnai[gain.snps.q,"LogRatio"] = cnai[gain.snps.q,"LogRatio"] + DAF_SCALING_FACTOR*cnai[gain.snps.q,"DAF"]
      ## ##   cnai[gain.snps.q,"LogRatio"] = filt.log2r[gain.snps.index] + DAF_SCALING_FACTOR*cnai[gain.snps.q,"DAF"]
      ## ## }
      
      ## if(length(loss.snps.q)>0){
      ##   cnai[loss.snps.q,"LogRatio"] = cnai[loss.snps.q,"LogRatio"] - DAF_SCALING_FACTOR*cnai[loss.snps.q,"DAF"]
      ##   ##cnai[loss.snps.q,"LogRatio"] = filt.log2r[loss.snps.index] - DAF_SCALING_FACTOR*cnai[loss.snps.q,"DAF"]
      ## }


## interpolate.hets.daf = function(sub.data,WINDOW_DAF){
##   positions = which(sub.data!=0 & !is.na(sub.data))
##   if(length(positions)>WINDOW_DAF){
##     dd = rollmedian(sub.data[positions],WINDOW_DAF,na.pad=T)
##     dd[is.na(dd)] = 0
##     for(ind.dd in 1:(length(positions)-1)){
##       ## Imputing threshold, if two LOH informative SNPs (hets. in normal)
##       ## are more than 100 homozygous SNPs apart, do not assign a DAF for the homozygous snps
##       ## between them
##       ## if(ind.dd < length(positions)-IMPUTED_SEGMENT_SNPS_SPAN){
##       if(ind.dd < length(positions)-IMPUTED_SEGMENT_SNPS_SPAN &
##          ind.dd > IMPUTED_SEGMENT_SNPS_SPAN ){
##         if(positions[ind.dd+IMPUTED_SEGMENT_SNPS_SPAN] - positions[ind.dd] < IMPUTED_SEGMENT_THRES &
##            positions[ind.dd] - positions[ind.dd-IMPUTED_SEGMENT_SNPS_SPAN] < IMPUTED_SEGMENT_THRES ){        
##           sub.data[positions[ind.dd]:positions[ind.dd+1]] = dd[ind.dd]
##         } else {
##           ## Zap the DAF of all snps within a poor het-SNPs covered region to 0
##           sub.data[positions[ind.dd]:positions[ind.dd+1]] = 0
##         }
##       } else {        
##         sub.data[positions[ind.dd]:positions[ind.dd+1]] = dd[ind.dd]        
##       }
##     }
##     ## First interval
##     sub.data[1:positions[1]] = dd[1]
##     ## Last interval
##     sub.data[positions[length(positions)]:length(sub.data)] = dd[length(positions)]
##   }
##   sub.data
## }

## interpolate.homs.log2r = function(daf.data,log2r.data){
##   indxs = which(daf.data!=0)    
##   ## Check for chromosome Y, which has all its SNPs having DAF of 0.0
##   if(length(daf.data)>=length(indxs)+2){
##     log2r.data[indxs] = NA ## Blank all heterozygous SNPs
##     log2r.data = na.approx(log2r.data,rule=2) ## Replace het SNPs with the log2r of neighboring hom SNPs
##   }
##   log2r.data
## }

## interpolate.hets.log2r = function(daf.data,log2r.data){
##   indxs = which(daf.data==0)
##   ## Check for chromosome Y, which has all its SNPs having DAF of 0.0
##   if(length(daf.data)>=length(indxs)+2){
##     log2r.data[indxs] = NA ## Blank all homozygous SNPs
##     log2r.data = na.approx(log2r.data,rule=2) ## Replace hom SNPs with the log2r of neighboring het SNPs
##   }
##   log2r.data
## }





## Plotting functions

    ## plot.logr.baf.profile.chr = function(chrarr,hets.logr.glad,homs.logr.glad,baf.data,hets.daf.glad){
    ##     for(chri in 1:length(chrarr)){
    ##       chr = chrarr[chri]
    ##       xmax = snps[snps$snp_id==snps_annotation[snps_annotation$chromosome==chr,"stop"],"position"]
    ##       sub.final.hets = hets.logr.glad[hets.logr.glad$Chromosome==chr & !is.na(hets.logr.glad$LogRatio),]
    ##       sub.final.homs = homs.logr.glad[homs.logr.glad$Chromosome==chr & !is.na(homs.logr.glad$LogRatio),]
    ##       plot(sub.final.homs$PosBase,sub.final.homs$LogRatio,
    ##            ylim=c(DOUBLE_LOSS_THRESHOLD,GAIN_THRESHOLD),
    ##            xlim=c(0,xmax),
    ##            col="gray",ylab="Log2R",xaxs="i",xaxt="n",las=2,pch=".")
    ##       title(line=0,xlab=paste("Chromosome",chr),col.lab="black")
    ##       axis(1,at=seq(1,max(sub.final.homs$PosBase),by=10000000),
    ##            lab=round(seq(1,max(sub.final.homs$PosBase),by=10000000)/1000000))
    ##       points(sub.final.hets$PosBase,sub.final.hets$LogRatio,col="black",pch=".")
    ##       abline(h=0,col="green")
    ##       lines(sub.final.hets$PosBase,sub.final.hets$Smoothing,col="red")
    ##       lines(sub.final.homs$PosBase,sub.final.homs$Smoothing,col="blue")
    ##                 sub.baf = baf.data[baf.data$Chromosome==chr,]
    ##       ## plot BAF and DAF segmentation
    ##       sub.baf = baf.data[baf.data$Chromosome==chr,]
    ##       sub.daf = hets.daf.glad[hets.daf.glad$Chromosome==chr & !is.na(hets.daf.glad$LogRatio),]
    ##       plot(sub.baf$PosBase,sub.baf$baf,ylim=c(0,1),xlim=c(0,xmax),col="black",
    ##            ylab="BAF",yaxs="i",yaxt="n",xaxs="i",xaxt="n",pch=".")
    ##       title(line=0,xlab=paste("Chromosome",chr),col.lab="black")
    ##       axis(1,at=seq(1,max(sub.baf$PosBase),by=10000000),
    ##            lab=round(seq(1,max(sub.baf$PosBase),by=10000000)/1000000))
    ##       axis(2,at=c(0,0.5,1), labels=c(0,0.5,1), tick=T,las=2)
    ##       abline(h=0.5,col="green")
    ##       lines(sub.daf$PosBase,sub.daf$Smoothing,col="red")
    ##     }
    ## }


    ## plot.profile.chr = function(chrarr,hets.logr.glad,homs.logr.glad){
    ##   if(length(chrarr)>1){
    ##     for(chri in 1:length(chrarr)){
    ##       chr = chrarr[chri]
    ##       xmax = snps[snps$snp_id==snps_annotation[snps_annotation$chromosome==chr,"stop"],"position"]
    ##       sub.final.hets = hets.logr.glad[hets.logr.glad$Chromosome==chr,]
    ##       sub.final.homs = homs.logr.glad[homs.logr.glad$Chromosome==chr,]
    ##       plot(sub.final.homs$PosBase,sub.final.homs$LogRatio,
    ##            ylim=c(DOUBLE_LOSS_THRESHOLD,GAIN_THRESHOLD),
    ##            xlim=c(0,xmax),
    ##            col="gray",ylab="Log2R",xaxs="i",xaxt="n",las=2,pch=".")
    ##       title(line=0,xlab=paste("Chromosome",chr),col.lab="black")
    ##       axis(1,at=seq(1,max(sub.final.homs$PosBase),by=10000000),
    ##            lab=round(seq(1,max(sub.final.homs$PosBase),by=10000000)/1000000))
    ##       points(sub.final.hets$PosBase,sub.final.hets$LogRatio,col="black",pch=".")
    ##       abline(h=0,col="green")
    ##       lines(sub.final.hets$PosBase,sub.final.hets$Smoothing,col="red")
    ##       lines(sub.final.homs$PosBase,sub.final.homs$Smoothing,col="blue")
    ##     }
    ##   } else {
    ##     chr = chrarr
    ##     xmax = snps[snps$snp_id==snps_annotation[snps_annotation$chromosome==chr,"stop"],"position"]
    ##     sub.final.hets = hets.logr.glad[hets.logr.glad$Chromosome==chr,]
    ##     sub.final.homs = homs.logr.glad[homs.logr.glad$Chromosome==chr,]
    ##     plot(sub.final.homs$PosBase,sub.final.homs$LogRatio,
    ##          ylim=c(DOUBLE_LOSS_THRESHOLD,GAIN_THRESHOLD),
    ##          xlim=c(0,xmax),
    ##          col="gray",ylab="Log2R",xaxs="i",xaxt="n",las=2,pch=".")        
    ##     title(line=0,xlab=paste("Chromosome",chr),col.lab="black")
    ##     axis(1,at=seq(1,max(sub.final.homs$PosBase),by=10000000),
    ##          lab=round(seq(1,max(sub.final.homs$PosBase),by=10000000)/1000000))        
    ##     points(sub.final.hets$PosBase,sub.final.hets$LogRatio,col="black",pch=".")
    ##     abline(h=0,col="green")
    ##     lines(sub.final.hets$PosBase,sub.final.hets$Smoothing,col="red")
    ##     lines(sub.final.homs$PosBase,sub.final.homs$Smoothing,col="blue")
    ##   }
    ## }

    ## plot.profile.baf.chr = function(chrarr,baf.data,hets.daf.glad){
    ##   if(length(chrarr)>1){
    ##     for(chri in 1:length(chrarr)){
    ##       chr = chrarr[chri]
    ##       xmax = snps[snps$snp_id==snps_annotation[snps_annotation$chromosome==chr,"stop"],"position"]
    ##       sub.baf = baf.data[baf.data$Chromosome==chr,]
    ##       sub.daf = hets.daf.glad[hets.daf.glad$Chromosome==chr,]
    ##       plot(sub.baf$PosBase,sub.baf$baf,ylim=c(0,1),xlim=c(0,xmax),col="black",
    ##            ylab="BAF",yaxs="i",yaxt="n",xaxs="i",xaxt="n",pch=".")
    ##       title(line=0,xlab=paste("Chromosome",chr),col.lab="black")
    ##       axis(1,at=seq(1,max(sub.baf$PosBase),by=20000000),
    ##            lab=round(seq(1,max(sub.baf$PosBase),by=20000000)/1000000))
    ##       axis(2,at=c(0,0.5,1), labels=c(0,0.5,1), tick=T,las=2)
    ##       abline(h=0.5,col="gray")
    ##       lines(sub.daf$PosBase,sub.daf$Smoothing,col="red")
    ##     }
    ##   } else {
    ##     chr = chrarr
    ##     xmax = snps[snps$snp_id==snps_annotation[snps_annotation$chromosome==chr,"stop"],"position"]
    ##     sub.baf = baf.data[baf.data$Chromosome==chr,]
    ##     sub.daf = hets.daf.glad[hets.daf.glad$Chromosome==chr,]
    ##     plot(sub.baf$PosBase,sub.baf$baf,ylim=c(0,1),xlim=c(0,xmax),col="black",
    ##          ylab="BAF",yaxt="n",yaxs="i",xaxs="i",xaxt="n",pch=".")
    ##     title(line=0,xlab=paste("Chromosome",chr),col.lab="black")
    ##     axis(1,at=seq(1,max(sub.baf$PosBase),by=20000000),
    ##          lab=round(seq(1,max(sub.baf$PosBase),by=20000000)/1000000))
    ##     axis(2,at=c(0,0.5,1), labels=c(0,0.5,1), tick=T,las=2)
    ##     abline(h=0.5,col="gray")
    ##     lines(sub.daf$PosBase,sub.daf$Smoothing,col="red")
    ##   }
    ## }

    ## plot.profile = function(chr,start.pos,stop.pos,hets.logr.glad,homs.logr.glad){
    ##       sub.final.hets = hets.logr.glad[hets.logr.glad$Chromosome==chr &
    ##         hets.logr.glad$PosBase>=start.pos &
    ##         hets.logr.glad$PosBase<=stop.pos,]
    ##       sub.final.homs = homs.logr.glad[homs.logr.glad$Chromosome==chr &
    ##         homs.logr.glad$PosBase>=start.pos &
    ##         homs.logr.glad$PosBase<=stop.pos,]
    ##       plot(sub.final.homs$PosBase,sub.final.homs$LogRatio,cex=0.5,ylim=c(DOUBLE_LOSS_THRESHOLD,GAIN_THRESHOLD),col="gray",ylab="Log2R")
    ##       title(line=0,xlab=paste("Chromosome",chr),col.lab="black")          
    ##       points(sub.final.hets$PosBase,sub.final.hets$LogRatio,cex=0.5,col="black")
    ##       abline(h=0,col="gray")
    ##       lines(sub.final.hets$PosBase,sub.final.hets$Smoothing,col="red")
    ##       lines(sub.final.homs$PosBase,sub.final.homs$Smoothing,col="blue")
    ## }

    ## baf.data = cbind(snps$snp_id, snps, sample_data[,"baf"])
    ## colnames(baf.data) = c("snp_id", "Chromosome","PosBase","PosOrder","baf")
    ## print("Printing segmentation profile")
    ## bitmap(paste("logrbaf_p",patient_id,"_s",sample_id,".png",sep=""),
    ##        width=6400,height=6400,units="px",pointsize=16,type="png256",taa=4)
    ## par(mfrow=c(48,1),mar=c(2,5,1,1))
    ## plot.logr.baf.profile.chr(chr=c(1:24),hets.logr.glad,homs.logr.glad,baf.data,hets.daf.glad)
    ## dev.off()
