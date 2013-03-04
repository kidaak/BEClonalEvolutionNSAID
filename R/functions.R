## Connects to a MySQL server
connectDb <- function(dbname, user=FALSE, password=FALSE, host=FALSE) {
    require("RMySQL")
    ## If password is not provided, load the default stored at secured location in /etc/keys/.mysqlpass
    if(!password){
      password <- as.character(read.table("/etc/keys/.mysqlpass", as.is=T)[1,1])
    }
    ## If user is not provided, load "root" as default
    if(!user){
      user <- "root"
    }
    ## Use localhost if host is not provided
    if(!host){
      host <- "localhost"
    }
    db <- dbConnect(MySQL(), dbname=dbname, user=user, password=password, host=host)
    db
}

## Disconnects from MySQL server
closeDb <- function(db){
    dbDisconnect(db)
}

## Removes MT, UNKNOWN, etc. chromosomes and replaces X, Y with 23, 24
filterUnmappedSNPs <- function(snps){
  ## Replace chromosome X and Y by 23 and 24
  snps[snps[,"chromosome"]=="X","chromosome"] = 23
  snps[snps[,"chromosome"]=="Y","chromosome"] = 24
  ## Remove all extra snps, those unmapped on chromosomes (position is 0)
  ## and those mapped on MT, UNKNOWN, XY chromosomes.
  snps = snps[snps[,"chromosome"] %in% c(1:24),]
  snps = snps[snps[,"position"]>0,]
  snps
}

## Removes all values from vector x that do not fall within specified lower and upper boundaries
match.conditions <- function(x, low.bounds, up.bounds){
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


blank.homs.log2r <- function(daf.data,log2r.data){
  indxs = which(!is.na(daf.data))
  ## Check for chromosome Y, which has all its SNPs having DAF of 0.0
  if(length(daf.data)>=length(indxs)+2){
    log2r.data[indxs] = NA ## Blank all heterozygous SNPs
  }
  log2r.data
}

blank.hets.log2r <- function(daf.data,log2r.data){
  indxs = which(is.na(daf.data))    
  ## Check for chromosome Y, which has all its SNPs having DAF of 0.0
  if(length(daf.data)>=length(indxs)+2){
    log2r.data[indxs] = NA ## Blank all homozygous SNPs
  }
  log2r.data
}

## Haplotype phasing function

phase.region <- function(x,vec1){
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


## Removes low intensity SNPs in reference samples
## that may be consitutive CNVs that are not present in the germline genotype of
## the individual
## Based on visualization, we can discard SNPs with R lower than 0.15
## Typically that discards ~2,576 SNPs
removeLowIntensitySNP <- function(con, referenceID, threshold=0.15){
  ## Load data 
  d <- dbGetQuery(con, paste("select snp_id, r, baf from data",
    "where sample_id=",referenceID,"order by snp_id"))

  ## Filter missing values for SNPs: those SNPs that have NA for R or baf
  ind <- which(is.na(d[,"r"]))
  filt.1 <- d[ind,]
  if(length(ind)>0){
    d <- d[-ind,]
  }

  ## Visualize SNPs with lowest signal intensity, 0.5%
  ## lowest.r = as.numeric(quantile(d[,"r"], probs=0.005, na.rm=TRUE))
  ## hist(d[d[,"r"]<lowest.r,"r"],breaks=100)
  ind <- which(d[,"r"] <= threshold)
  filt.2 <- d[ind,]
  if(length(ind)>0){
    d <- d[-ind,]
  }

  ## Apply only filter 1 and filter 2 to SNPs
  all.filt <- rbind(filt.1,filt.2)
  all.filt <- cbind(rep(referenceID,nrow(all.filt)),all.filt[,"snp_id"])
  colnames(all.filt) <- c("sample_id","snp_id")
  res <- as.data.frame(all.filt)
  ## delete previous rows
  rs <- dbSendQuery(con,paste("delete from sample_filtered_snps where sample_id=",referenceID))
  dbWriteTable(con, "sample_filtered_snps", res, row.names=F, append=TRUE)
}


segmentSamplesOfAllPatients <- function(con){
  patients = dbGetQuery(con, "select distinct(patient_id) from acs_paired_samples")
  for(i in 1:nrow(patients)){   
    segmentSamplesOfPatient(con=con, patientID=patients[i,1])
  }
}

segmentSamplesOfPatient <- function(con, patientID, filterSNP=TRUE){
  
  samples = dbGetQuery(con, paste("select * from acs_paired_samples where patient_id=",patientID))
  
  for(i in 1:nrow(samples)){

    sample_id = samples[i,"sample_id"]
    reference_id = samples[i,"reference_id"]

    segment(con=con, sample_id, reference_id, filterSNP=filterSNP)
    
  } 
}




segmentTechnicalReplicates <- function(con){
  ## Patient 662, 2 blood samples
  bld.662.1 = 258
  bld.662.2 = 261
  
  ## Patient 437
  bld.437.1 = 168
  bld.437.2 = 332

  segment(con=con, bld.662.2, bld.662.1, filterSNP=TRUE)
  segment(con=con, bld.437.1, bld.437.2, filterSNP=TRUE)


}

callStatesTechnicalReplicates <- function(con){
  ## Patient 662, 2 blood samples
  bld.662.1 = 258
  bld.662.2 = 261
  
  ## Patient 437
  bld.437.1 = 168
  bld.437.2 = 332

  callSingleSampleFeatureStates(con=con, bld.662.2, bld.662.1, patient_id=662)
  callSingleSampleFeatureStates(con=con, bld.437.1, bld.437.2, patient_id=437)
}



segment <- function(con, sample_id, reference_id, filterSNP=TRUE){
  require("GLAD")
  
  ## ###############
  ## Constants
  ## ###############

  ## Compute log2 R ratio
  ## Add pseudo-signal to prevent -Inf and NaNs
  DOUBLE_LOSS_THRESHOLD = -3
  GAIN_THRESHOLD = 1
  NORMAL_AB_LOWER_BAF_THRESHOLD = 0.33
  NORMAL_AB_UPPER_BAF_THRESHOLD = 0.66
  PSEUDO = 0.000001
  CORFACTOR = 0.53
  
  ##MIN_REQUIRED_HET_DENSITY = 1/20000 ## One het SNP every 10 kb

  
  ## #################
  ## This has to be a pass-by-reference code 
  ## ####################

  ## Load chr, start SNP, stop SNP annotation
  snps_annotation = dbGetQuery(con, "select * from snps_annotation")
  ## Fix X=23 Y=24
  snps_annotation[23,1]=23
  snps_annotation[24,1]=24
  
  ## Get snps
  snps = dbGetQuery(con, "select chromosome, position, snp_id from snps order by chromosome, position")
  snps = filterUnmappedSNPs(snps)
  
  ## ###########################
  ## End of pass-by-reference code snippet
  ## Consider biobase asasydata storagemode lockedenvironment to provide pass by reference
  ## ##########################

  ## Check if the sample is slotted into MEMORY
  check = dbGetQuery(con,paste("select * from data_memory where sample_id=",sample_id,"limit 1"))
  if(length(check)>0){
    sample_data=dbGetQuery(con,paste("select r, baf, snp_id from data_memory where sample_id=",
      sample_id,"order by snp_id"))
  } else {
    sample_data=dbGetQuery(con,paste("select r, baf, snp_id from data where sample_id=",
      sample_id,"order by snp_id"))
  }
  
  check = dbGetQuery(con,paste("select * from data_memory where sample_id=",reference_id,"limit 1"))
  if(length(check)>0){
    reference_data=dbGetQuery(con,paste("select r, baf, snp_id from data_memory where sample_id=",
      reference_id,"order by snp_id"))
  } else {
    reference_data=dbGetQuery(con,paste("select r, baf, snp_id from data where sample_id=",
      reference_id,"order by snp_id"))
  }
  
  ## ###############
  ## Pre-processing
  ## ###############
  
  ## Filter STEP
  if(filterSNP){
    ## Get the bad probes from the reference and remove them from both reference and sample
    snp_ids.to.filter=dbGetQuery(con,paste("select snp_id from sample_filtered_snps where sample_id=",
      reference_id,"order by snp_id"))
    
    ## Get rid of some 
    sample_data = sample_data[which(!sample_data[,"snp_id"]%in%snp_ids.to.filter[,1]),]
    reference_data = reference_data[which(!reference_data[,"snp_id"]%in%snp_ids.to.filter[,1]),]

    ## 1. Remove all probes that have low signal in the reference sample for the patient
    ## These are CNV probes that the patient lacks, or genomic microregions that are missing
    snps = snps[which(!snps[,"snp_id"]%in%snp_ids.to.filter[,1]),]

  } else {
    ## No filtering    
  }
  
  ## 2. Take the Log2 (R sample / R reference) and
  ##    take the difference between B allele frequencies of sample and reference
  ##    or "DAF" - difference in allelic frequency
  ## Shift reference BAF data to get it centered at 0.5
  ## Normalize datasets, center at 0.55
  
  ## if(prev_reference_id!=reference_id){  ## This corrects BAF in reference, make sure this is done only once
  ##  reference_data[reference_data[,"baf"]<=CORFACTOR,"baf"]=reference_data[reference_data[,"baf"]<=CORFACTOR,"baf"]*0.5/CORFACTOR
  ##  reference_data[reference_data[,"baf"]>CORFACTOR,"baf"]=0.5+0.5*(reference_data[reference_data[,"baf"]>CORFACTOR,"baf"]-CORFACTOR)/(1-CORFACTOR)
  ## }  ## Correct BAF of reference

  ## Correct BAF in reference
  reference_data[reference_data[,"baf"]<=CORFACTOR,"baf"]=reference_data[reference_data[,"baf"]<=CORFACTOR,"baf"]*0.5/CORFACTOR
  reference_data[reference_data[,"baf"]>CORFACTOR,"baf"]=0.5+0.5*(reference_data[reference_data[,"baf"]>CORFACTOR,"baf"]-CORFACTOR)/(1-CORFACTOR)
  
  ## Correct BAF in sample
  sample_data[sample_data[,"baf"]<=CORFACTOR,"baf"]=sample_data[sample_data[,"baf"]<=CORFACTOR,"baf"]*0.5/CORFACTOR
  sample_data[sample_data[,"baf"]>CORFACTOR,"baf"]=0.5+0.5*(sample_data[sample_data[,"baf"]>CORFACTOR,"baf"]-CORFACTOR)/(1-CORFACTOR)
    
  ##
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
}



## #################################
## Single sample feature state calling
## #################################
callSingleSampleFeatureStates = function(con, sample_id, reference_id, patient_id){

  MIN_REQUIRED_HET_SNPS = 15 ## At least 5 SNPs confirm an LOH event

  ## Load chr, start SNP, stop SNP annotation
  snps_annotation = dbGetQuery(con, "select * from snps_annotation")
  ## Fix X=23 Y=24
  snps_annotation[23,1]=23
  snps_annotation[24,1]=24

  ## Get snps
  snps = dbGetQuery(con, "select chromosome, position, snp_id from snps order by chromosome, position")
  snps = filterUnmappedSNPs(snps)
  
  ## Get state calling thresholds
  sthres = dbGetQuery(con, "select * from state_calling_thresholds")
  
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
  

  breakpoints = NULL
  d=dbGetQuery(con,paste("select * from feature_state_probabilities f inner join sample_features s
                          on s.feature_id=f.feature_id
                          where s.sample_id=",sample_id," and s.reference_id=",reference_id))

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
  
  ## Gather all breakpoints
  breakpoints = sort(unique(c(breakpoints,d[,"start_snp_id"],d[,"stop_snp_id"])))

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



  d=dbGetQuery(con,paste("select * from feature_state_probabilities f inner join sample_features s
                          on s.feature_id=f.feature_id
                          where s.sample_id=",sample_id," and s.reference_id=",reference_id))
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
 
  ## #############################
  ## Skip compressing features 
  ## #############################

}



## Slots samples into MEMORY type database in the MySQL server
slotSamplesInMemory = function(con, sample_ids){

  ## Delete samples from memory
  bla = dbSendQuery(con, "delete from data_memory")

  for(i in 1:length(sample_ids)){
    print(paste("Loading sample",i))
    sample_id = sample_ids[i]
    dbSendQuery(con,paste("select data.sample_id,data.snp_id,data.r,data.baf from data where sample_id= ",
                          sample_id," INTO OUTFILE '/tmp/temptbl",i,".txt'",sep=""))
    dbSendQuery(con,paste("load data infile '/tmp/temptbl",i,".txt' into table data_memory",sep=""))
  }                                      

}



