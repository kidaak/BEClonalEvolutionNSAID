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


