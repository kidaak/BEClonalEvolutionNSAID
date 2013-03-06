## Connect to mySQL database
library(RMySQL)

passwd = as.character(read.table("/etc/keys/.mysqlpass",as.is=T)[1,1])
con = dbConnect(MySQL(), user="root", password=passwd, dbname="sbep", host="localhost")
## Get patient_id and the reference sample id
patients = t(dbGetQuery(con, "select distinct(patient_id) from acs_paired_samples"))

all.period.0 = NULL
all.period.1 = NULL
all.period.2 = NULL
period.0.scale = NULL
period.1.scale = NULL
period.2.scale = NULL
period.0.all.cnaloh = NULL
period.1.all.cnaloh = NULL
period.2.all.cnaloh = NULL
nc.users = c(1,3,4,5,6,7,8,9,10,12,13)
cf.users = c(2,11)

cnaloh.distrib.all = NULL
cnaloh.distrib.new = NULL
cnaloh.distrib.extinct = NULL
cnaloh.distrib.final = NULL

patient.on.nsaids.time = NULL
patient.off.nsaids.time = NULL


## Compute the average size of an alteration
for(pind in 1:length(patients)){
  patient_id = patients[pind]
  print(paste("Patient",patient_id))
  samples = dbGetQuery(con, paste("select * from acs_paired_samples where patient_id=",patient_id,
    "order by endoscopy_date, biopsy_level"))

  ## Calculate time
  if(pind %in% nc.users){  
    patient.on.nsaids.time = c(patient.on.nsaids.time,unique(samples[,"second_sampling_period_duration"]))
    patient.off.nsaids.time = c(patient.off.nsaids.time,unique(samples[,"first_sampling_period_duration"]))
  } else {
    patient.on.nsaids.time = c(patient.on.nsaids.time,unique(samples[,"first_sampling_period_duration"]))
    patient.off.nsaids.time = c(patient.off.nsaids.time,unique(samples[,"second_sampling_period_duration"]))
  }

  all.features = dbGetQuery(con, paste("select * from compressed_patient_feature_states",
    "where patient_id=",patient_id," order by start_snp_id, sample_id"))
  ## Re-order features by endoscopy date, biopsy level
  samples.order = rep(order(samples[,"sample_id"]),nrow(all.features)/nrow(samples))
  all.features.order.index = order(all.features[,"start_snp_id"],samples.order)
  all.features = all.features[all.features.order.index,]
  
  feat.states = matrix(all.features[,"state"],ncol=nrow(samples),byrow=T)
  feat.phased.states = matrix(all.features[,"state_phased"],ncol=nrow(samples),byrow=T)
  feat.binary.states = matrix(all.features[,"state"],ncol=nrow(samples),byrow=T)
  feat.binary.states[feat.binary.states!="AB"] = 1
  feat.binary.states[feat.binary.states!=1] = 0
  feat.size = matrix(all.features[,"stop_position"]-all.features[,"start_position"]+1,ncol=nrow(samples),byrow=T)
  feat.chr = matrix(all.features[,"chromosome"],ncol=nrow(samples),byrow=T)
  feat.start = matrix(all.features[,"start_position"],ncol=nrow(samples),byrow=T)
  feat.stop = matrix(all.features[,"stop_position"],ncol=nrow(samples),byrow=T)

  ## Sort according to the binary pattern
  patterns = apply(feat.binary.states,1,function(x) paste(x,collapse=""))
  order.index = order(patterns)

  sgas.summary = NULL
  sgas.summary = cbind(rep("",ncol(samples)),rep("",ncol(samples)),rep("",ncol(samples)),rep("",ncol(samples))) 
  for(i in 1:nrow(samples)){
    sgas.summary = cbind(sgas.summary,t(samples[i,]))
  }
  sgas.summary[1:nrow(sgas.summary),4] = names(samples)
  sgas.summary = rbind(sgas.summary,rep("",ncol(sgas.summary)))
  sgas.summary[nrow(sgas.summary),c(1:4)] = c("chromosome","start","stop","size")

  feat.binary.states = feat.binary.states[order.index,]
  feat.states = feat.states[order.index,]
  feat.phased.states = feat.phased.states[order.index,]
  feat.chr = feat.chr[order.index,]
  feat.start = feat.start[order.index,]
  feat.stop = feat.stop[order.index,]
  feat.size = feat.size[order.index,]

  colnames(sgas.summary) = c(1:ncol(sgas.summary))
  colnames(feat.binary.states) = c(1:ncol(feat.binary.states))
  colnames(feat.states) = c(1:ncol(feat.states))
  colnames(feat.phased.states) = c(1:ncol(feat.phased.states))
  colnames(feat.chr) = c(1:ncol(feat.chr))
  colnames(feat.start) = c(1:ncol(feat.start))
  colnames(feat.stop) = c(1:ncol(feat.stop))
  colnames(feat.size) = c(1:ncol(feat.size))
  
  binary.sga = rbind(sgas.summary,cbind(feat.chr[,1],feat.start[,1],feat.stop[,1],feat.size[,1],feat.binary.states))
  unphased.sga = rbind(sgas.summary,cbind(feat.chr[,1],feat.start[,1],feat.stop[,1],feat.size[,1],feat.states))
  phased.sga = rbind(sgas.summary,cbind(feat.chr[,1],feat.start[,1],feat.stop[,1],feat.size[,1],feat.phased.states))

  feat.binary.states = data.frame(as.numeric(feat.chr[,1]),feat.start[,1],feat.stop[,1],feat.binary.states)



  ## Baseline, off, on nsaids
  samples.period.0 = which(samples[,"endo_age_from_baseline"]==0)
  samples.period.1 = which(samples[,"endo_age_from_baseline"]<samples[,"nsaids_transition_age_from_baseline"] &
                           samples[,"endo_age_from_baseline"]>0)
  samples.period.2 = which(samples[,"endo_age_from_baseline"]>samples[,"nsaids_transition_age_from_baseline"])
  
  ## Plot total SGA
  feats.period.0 = feat.binary.states[,c(1:3,samples.period.0+3)]
  feats.period.1 = feat.binary.states[,c(1:3,samples.period.1+3)]
  feats.period.2 = feat.binary.states[,c(1:3,samples.period.2+3)]

  feat.freq.period.0 = apply(feats.period.0, 1, function(x) round(sum(as.numeric(x[4:length(x)]),na.rm=T)/(length(x)-3),digits=2))
  feat.freq.period.1 = apply(feats.period.1, 1, function(x) round(sum(as.numeric(x[4:length(x)]),na.rm=T)/(length(x)-3),digits=2))
  feat.freq.period.2 = apply(feats.period.2, 1, function(x) round(sum(as.numeric(x[4:length(x)]),na.rm=T)/(length(x)-3),digits=2))
  
  ## Select CNA/LOH events that are novel for the periods on/off nsaids 
  feats.period.1 = feats.period.1[feat.freq.period.0==0,]
  feats.period.2 = feats.period.2[feat.freq.period.0==0 & feat.freq.period.1==0,]
  
  new.endo.at.bx = c(TRUE,diff(samples[,"endo_age_from_baseline"])>0)
  
  for(i in 4:ncol(feat.binary.states)){
    feats = NULL
    ## Calculate the time elapsed since last endo
    endo.age.current.bx = samples[(i-3),"endo_age_from_baseline"]
    previous.endos.age = samples[which(samples[,"endo_age_from_baseline"]<endo.age.current.bx),"endo_age_from_baseline"]
    if(length(previous.endos.age)==0){
      endo.time = 1 ## Still baseline biopsies
    } else {
      endo.time = endo.age.current.bx - max(previous.endos.age)
    }

    if((i-3)%in%samples.period.0){
      sid=0
      feats = feat.binary.states[,c(1:3,i)]
      ## Save all SGAs present at baseline biopsies to "feat.freq.period.0"
      ## to mark them as observed
      feat.freq.period.0[feat.binary.states[,i]==1] = 1
    }
    if((i-3)%in%samples.period.1){
      sid=1
      feats = feat.binary.states[,c(1:3,i)]
      ## These are all the new SGAs      
      feats = feats[feat.freq.period.0==0,]

      ## Upon a new endoscopy, add all biopsies from the previous endoscopy to the list of existing SGAs
      if(new.endo.at.bx[i-3]){
        feats.2 = feat.binary.states[,which(samples[,"endo_age_from_baseline"]==max(previous.endos.age))+3]
        if(!is.null(ncol(feats.2))){
          feats.3 = apply(feats.2, 1, function(x) any(as.logical(as.numeric(as.character(x)))) )
        } else {
          feats.3 = as.logical(as.numeric(as.character(feats.2)))
        }
        feat.freq.period.0[feats.3] = 1
      }
    }
    if((i-3)%in%samples.period.2){
      sid=2
      feats = feat.binary.states[,c(1:3,i)]
      ## These are all the new SGAs      
      feats = feats[feat.freq.period.0==0,]
      
      if(new.endo.at.bx[i-3]){
        feats.2 = feat.binary.states[,which(samples[,"endo_age_from_baseline"]==max(previous.endos.age))+3]
        if(!is.null(ncol(feats.2))){
          feats.3 = apply(feats.2, 1, function(x) any(as.logical(as.numeric(as.character(x)))) )
        } else {
          feats.3 = as.logical(as.numeric(as.character(feats.2)))
        }
        feat.freq.period.0[feats.3] = 1
      }
    }    
    ## Divide by time elapsed since last endo
    cnaloh.distrib.new = rbind(cnaloh.distrib.new,cbind(pind,sid,t(table(cut(log10(
      feats[feats[,4]==1,3]-feats[feats[,4]==1,2]),
      breaks = c(-Inf,seq(2,7,by=1),Inf),right=T)))/endo.time
      ))
  }

  ## Calculate extinct events, events that are present at 0% frequency in LAST endo
  samples.period.0 = which(samples[,"endo_age_from_baseline"]<samples[,"nsaids_transition_age_from_baseline"])
  samples.period.1 = which(samples[,"endo_age_from_baseline"]>samples[,"nsaids_transition_age_from_baseline"] & 
                           samples[,"endo_coalescent_age"]>0)
  samples.period.2 = which(samples[,"endo_coalescent_age"]==0)

  feats.period.0 = feat.binary.states[,c(1:3,samples.period.0+3)]
  feats.period.1 = feat.binary.states[,c(1:3,samples.period.1+3)]
  feats.period.2 = feat.binary.states[,c(1:3,samples.period.2+3)]

  feat.freq.period.0 = apply(feats.period.0, 1, function(x) round(sum(as.numeric(x[4:length(x)]),na.rm=T)/(length(x)-3),digits=2))
  feat.freq.period.1 = apply(feats.period.1, 1, function(x) round(sum(as.numeric(x[4:length(x)]),na.rm=T)/(length(x)-3),digits=2))
  feat.freq.period.2 = apply(feats.period.2, 1, function(x) round(sum(as.numeric(x[4:length(x)]),na.rm=T)/(length(x)-3),digits=2))



  ##new.endo.at.bx = c(TRUE,diff(samples[,"endo_age_from_baseline"])>0)
  
  for(i in 4:ncol(feat.binary.states)){
    feats = NULL
    ## Calculate the time elapsed since last endo
    endo.age.current.bx = samples[(i-3),"endo_age_from_baseline"]
    ## Get the ages of biopsies in all future endoscopies
    future.bxs.age = samples[which(samples[,"endo_age_from_baseline"]>endo.age.current.bx),"endo_age_from_baseline"]
    if(length(future.bxs.age)==0){
      endo.time = 1 ## Still baseline biopsies
    } else {
      endo.time = min(future.bxs.age) - endo.age.current.bx
    }

    if((i-3)%in%samples.period.0){
      sid=0
      feats = feat.binary.states[,c(1:3,i)]
      ## Calculate all future events in all future endos
      feats.2 = feat.binary.states[,which(samples[,"endo_age_from_baseline"]%in%future.bxs.age)+3]
      if(!is.null(ncol(feats.2))){
        feats.3 = apply(feats.2, 1, function(x) any(as.logical(as.numeric(as.character(x)))) )
      } else {
        feats.3 = as.logical(as.numeric(as.character(feats.2)))
      }
      
      ## These are all the new SGAs that went extinct in future biopsies in future endos      
      feats = feats[feats.3==FALSE,]
    }
    if((i-3)%in%samples.period.1){
      sid=1
      feats = feat.binary.states[,c(1:3,i)]
      ## Calculate all future events in all future endos
      feats.2 = feat.binary.states[,which(samples[,"endo_age_from_baseline"]%in%future.bxs.age)+3]
      if(!is.null(ncol(feats.2))){
        feats.3 = apply(feats.2, 1, function(x) any(as.logical(as.numeric(as.character(x)))) )
      } else {
        feats.3 = as.logical(as.numeric(as.character(feats.2)))
      }      
      ## These are all the new SGAs that went extinct in future biopsies in future endos      
      feats = feats[feats.3==FALSE,]
    }
    if((i-3)%in%samples.period.2){
      sid=2
      feats = feat.binary.states[,c(1:3,i)]
      ## These are all the new SGAs      
      feats = feats[feat.freq.period.0==0,]
    }    
    ## Divide by time elapsed since last endo
    cnaloh.distrib.extinct = rbind(cnaloh.distrib.extinct,cbind(pind,sid,t(table(cut(log10(
      feats[feats[,4]==1,3]-feats[feats[,4]==1,2]),
      breaks = c(-Inf,seq(2,7,by=1),Inf),right=T)))/endo.time
      ))
  }


}

error.bar = function(x, y, upper, lower=upper, code,length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  ##arrows(x,y+upper, x, y-lower, angle=90, code=code, length=length, ...)
  arrows(x,y+upper, x, y, angle=90, code=code, length=length, ...)
}


library(Cairo)

## Calculate time



pdf("kostadinov_fig_new_sga_events_aug31.pdf",width=3,height=3,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),]
p1.mean = apply(p1[,3:ncol(p1)], 2, function(x) mean(x))
p1.sd = apply(p1[,3:ncol(p1)], 2, function(x) sd(x))
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),]
p2.mean = apply(p2[,3:ncol(p2)], 2, function(x) mean(x))
p2.sd = apply(p2[,3:ncol(p2)], 2, function(x) sd(x))
bp = barplot(rbind(p1.mean,p2.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGAs per biopsy per year",
        xlab = "Lesion size (bp)",names.arg=rep("",length(p1.mean)), ylim=c(0,70),las=1)
## Plot 95% CI
error.bar(bp,rbind(p1.mean,p2.mean),rbind(1.96*p1.sd/sqrt(nrow(p1)),1.96*p2.sd/sqrt(nrow(p2))),code=1,length=0.02)

axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))
legend("topright",legend=c("Off NSAID","On NSAID","95% CI"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,1),merge=T,box.col="black",cex=0.8)

for(i in 3:9){
  pval = t.test(p1[,i],p2[,i],alt="greater")$p.val
  text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+1,cex=0.8)
  print(pval)
}

dev.off()



## BOXPLOTS of categories
pdf("kostadinov_fig_new_sga_events_boxplot_aug31.pdf",width=3,height=3,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),3:ncol(cnaloh.distrib.new)]
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),3:ncol(cnaloh.distrib.new)]
boxplot.default(p1[,1],p2[,1],
        p1[,2],p2[,2],
        p1[,3],p2[,3],
        p1[,4],p2[,4],
        p1[,5],p2[,5],
        p1[,6],p2[,6],
        p1[,7],p2[,7],
        ylim=c(0,60),
        col=c("white","gray"),las=2,axes=F,outcex=0.5,whisklty="solid")
axis(1,at=seq(0.5,14.5,by=2),lty=1, c(0,expression(10^2),expression(10^3),expression(10^4),
                           expression(10^5),expression(10^6),expression(10^7),expression(10^8)) )
axis(2,at=seq(0,60,by=20))
box()
mtext(side=2,at=30,line=3, "Number of SGA per biopsy per year")
mtext(side=1,line=2, "Lesion size (bp)")


legend("topright", legend=c("Off NSAID","On NSAID","5th/95th percentiles","Median"),
       fill=c("white","darkgray",NA,NA),
       border=c("black","black",NA,NA),
       lty=c(0,0,1,1),lwd=c(0,0,1,3),merge=T,box.col="black",cex=0.8)

for(i in 1:7){
  pval = wilcox.test(p1[,i],p2[,i],alt="greater")$p.val
  ##text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+1,cex=0.8)
  print(pval)
}

dev.off()


## Boxplot correctly

matlab.boxplot <- function(m){
  m  <- as.matrix(m)
  bp <- boxplot(data.frame(m), plot=FALSE)
  bp$stats <- apply( m, 2, function(x) 
                   quantile(x, c(0.05,0.25, 0.5, 0.75, 0.95), na.rm=T) ) 
  tmp <- apply( m, 2, function(x){
    under <- x[ which( x < quantile(x, 0.05, na.rm=T) ) ]
    over  <- x[ which( x > quantile(x, 0.95, na.rm=T) ) ]
    return( c(under, over) )
  })   # always a matrix in this case
  bp$out   <- c(tmp)
  
  bp$group <- rep(1:ncol(tmp), each=nrow(tmp))
  
  bxp(bp)
}


##myboxplot.stats = function (x, coef = NULL, do.conf = TRUE, do.out = TRUE)
myboxplot.stats = function (x, coef = NULL, do.out = TRUE) 
{
  nna = !is.na(x)
  n = sum(nna)
  stats = quantile(x, c(.05,.25,.5,.75,.95), na.rm = TRUE)
  iqr = diff(stats[c(2, 4)])
  out = x < stats[1] | x > stats[5]
  ## conf = if (do.conf) 
  ##   stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n)
  ##list(stats = stats, n = n, conf = conf, out = x[out & nna])
  list(stats = stats, n = n, out = x[out & nna])
}

myboxplot.stats = function (x, coef = 1.5, do.conf = TRUE, do.out = TRUE) 
{
    if (coef < 0) 
        stop("'coef' must not be negative")
    nna <- !is.na(x)
    n <- sum(nna)
    stats <- stats::fivenum(x, na.rm = TRUE)
    iqr <- diff(stats[c(2, 4)])
    if (coef == 0) 
        do.out <- FALSE
    else {
        out <- if (!is.na(iqr)) {
            x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * 
                iqr)
        }
        else !is.finite(x)
        if (any(out[nna], na.rm = TRUE)) 
            stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- if (do.conf) 
        stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
    list(stats = stats, n = n, conf = conf, out = if (do.out) x[out & 
        nna] else numeric(0L))
}


boxplot.default = function (x, ..., range = 1.5, width = NULL, varwidth = FALSE, 
    notch = FALSE, outline = TRUE, names, plot = TRUE, border = par("fg"), 
    col = NULL, log = "", pars = list(boxwex = 0.8, staplewex = 0.5, 
        outwex = 0.5), horizontal = FALSE, add = FALSE, at = NULL) 
{
    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names)) 
        attributes(args)$names != ""
    else rep(FALSE, length.out = length(args))
    groups <- if (is.list(x)) 
        x
    else args[!namedargs]
    if (0L == (n <- length(groups))) 
        stop("invalid first argument")
    if (length(class(groups))) 
        groups <- unclass(groups)
    if (!missing(names)) 
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names"))) 
            attr(groups, "names") <- 1L:n
        names <- attr(groups, "names")
    }
    cls <- sapply(groups, function(x) class(x)[1L])
    cl <- if (all(cls == cls[1L])) 
        cls[1L]
    else NULL
    for (i in 1L:n) groups[i] <- list(myboxplot.stats(unclass(groups[[i]]), 
        range))
    stats <- matrix(0, nrow = 5L, ncol = n)
    conf <- matrix(0, nrow = 2L, ncol = n)
    ng <- out <- group <- numeric(0L)
    ct <- 1
    for (i in groups) {
        stats[, ct] <- i$stats
        conf[, ct] <- i$conf
        ng <- c(ng, i$n)
        if ((lo <- length(i$out))) {
            out <- c(out, i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct + 1
    }
    if (length(cl) && cl != "numeric") 
        oldClass(stats) <- cl
    z <- list(stats = stats, n = ng, conf = conf, out = out, 
        group = group, names = names)
    if (plot) {
        if (is.null(pars$boxfill) && is.null(args$boxfill)) 
            pars$boxfill <- col
        do.call("bxp", c(list(z, notch = notch, width = width, 
            varwidth = varwidth, log = log, border = border, 
            pars = pars, outline = outline, horizontal = horizontal, 
            add = add, at = at), args[namedargs]))
        invisible(z)
    }
    else z
}





pdf("kostadinov_fig_new_sga_events_boxplot_sep11.pdf",width=3,height=3,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),3:ncol(cnaloh.distrib.new)]
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),3:ncol(cnaloh.distrib.new)]
## Omit outliers
bp = boxplot(p1[,1],p2[,1],
        p1[,2],p2[,2],
        p1[,3],p2[,3],
        p1[,4],p2[,4],
        p1[,5],p2[,5],
        p1[,6],p2[,6],
        p1[,7],p2[,7],
        ylim=c(-2,100),
        col=c("white","gray"),las=2,axes=F,outline=T,outcex=0.7,whisklty="solid",medlty="solid",medlwd=1,
        at=sort(c(seq(1.2,13.2,by=2),seq(1.8,13.8,by=2))),
        boxwex=0.5, plot=F)
m = list(p1[,1],p2[,1],
        p1[,2],p2[,2],
        p1[,3],p2[,3],
        p1[,4],p2[,4],
        p1[,5],p2[,5],
        p1[,6],p2[,6],
        p1[,7],p2[,7])
bpq = lapply(m, function(x) quantile(x, c(0.05,0.25, 0.5, 0.75, 0.95), na.rm=T))
m2 = matrix(unlist(bpq),ncol=5,byrow=T)
bp$stats <- t(m2)
tmp <- lapply( m, function(x){
  under <- x[ which( x < quantile(x, 0.05, na.rm=T) ) ]
  over  <- x[ which( x > quantile(x, 0.95, na.rm=T) ) ]
  return( c(under, over) )
})   # always a matrix in this case

bp$out   <- c(unlist(tmp))
outlengths <- unlist(lapply(tmp, length))
bp$group <- rep(1:length(outlengths), times=outlengths)
bxp(bp,ylim=c(-2,100),
        boxfill=c("white","gray"),las=2,axes=F,staplewex=1,outline=T,outcex=0.5,whisklty="solid",medlty="solid",medlwd=1,
        at=sort(c(seq(1.2,13.2,by=2),seq(1.8,13.8,by=2))),
        boxwex=0.5)

axis(1,at=seq(0.5,14.5,by=2),lty=1, c(0,expression(10^2),expression(10^3),expression(10^4),
                           expression(10^5),expression(10^6),expression(10^7),expression(10^8)) )
axis(2,at=seq(0,100,by=20),lab=c(seq(0,80,by=20),">=100"),las=1)

for(i in 1:14){
  outliers = tmp[[i]][tmp[[i]]>100]
  if(length(outliers)>0){
    points(i-0.3+(1:length(outliers)/(length(outliers)+1)),rep(100,length(outliers)),cex=0.5)
  }
}

box()
mtext(side=2,line=3, "Number of SGA per biopsy per year")
mtext(side=1,line=2, "Lesion size (bp)")

legend("topright", legend=c("Off NSAID","On NSAID","Wilcox p<.05"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,0),lwd=c(0,0,1),pch=c(NA,NA,8),merge=T,box.col="black",cex=0.8,bg="white")

for(i in 1:7){
  pval = wilcox.test(p1[,i],p2[,i],alt="greater")$p.val
  if(pval<0.05){
    ## segments(i*2-1,-2,i*2,-2)
    ## segments(i*2-1,-2,i*2-1,-0.5)
    ## segments(i*2,-2,i*2,-0.5)
    points(i*2-0.5,-2,pch=8)
  }
  ##text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+1,cex=0.8)
  print(pval)
}

dev.off()


## R whiskers

pdf("kostadinov_fig_new_sga_events_boxplot_sep12.pdf",width=3,height=3,pointsize=10)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),3:ncol(cnaloh.distrib.new)]
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),3:ncol(cnaloh.distrib.new)]
## Omit outliers
bp = boxplot(p1[,1],p2[,1],
        p1[,2],p2[,2],
        p1[,3],p2[,3],
        p1[,4],p2[,4],
        p1[,5],p2[,5],
        p1[,6],p2[,6],
        p1[,7],p2[,7],
        ylim=c(-2,25),
        col=c("white","gray"),las=2,axes=F,outline=T,outcex=0.7,whisklty="solid",medlty="solid",medlwd=1,
        at=sort(c(seq(1.2,13.2,by=2),seq(1.8,13.8,by=2))),
        boxwex=0.5, plot=F)
bxp(bp,ylim=c(-2,25),
        boxfill=c("white","darkgray"),las=2,axes=F,staplewex=1,outline=F,outcex=0.5,whisklty="solid",medlty="solid",medlwd=1,
        at=sort(c(seq(1.1,13.1,by=2),seq(1.9,13.9,by=2))),
        boxwex=0.6)

axis(1,at=seq(0.5,14.5,by=2),lty=1, c(0,expression(10^2),expression(10^3),expression(10^4),
                           expression(10^5),expression(10^6),expression(10^7),expression(10^8)) )
axis(2,at=seq(0,30,by=5),las=1)

box()
mtext(side=2,line=3, "Number of SGA per biopsy per year")
mtext(side=1,line=2, "Lesion size (bp)")

legend("topright", legend=c("Off NSAID","On NSAID","Wilcox p<.05"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,0),lwd=c(0,0,1),pch=c(NA,NA,8),merge=T,box.col="black",cex=0.5,bg="white")

for(i in 1:7){
  pval = wilcox.test(p1[,i],p2[,i],alt="greater")$p.val
  if(pval<0.05){
    ## segments(i*2-1,-2,i*2,-2)
    ## segments(i*2-1,-2,i*2-1,-0.5)
    ## segments(i*2,-2,i*2,-0.5)
    points(i*2-0.5,-2,pch=8,cex=0.5)
  }
  ##text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+1,cex=0.8)
  print(pval)
}

dev.off()

## Boxplots use quantile(type=7 by default)












a = c(1:100)
boxplot(a,col=c("white"),las=2,axes=F,outline=F,whisklty="solid",medlty="solid",medlwd=1)






pdf("kostadinov_fig_extinct_sga_events_boxplot_sep12.pdf",width=3,height=3,pointsize=10)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),3:ncol(cnaloh.distrib.extinct)]
p2 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),3:ncol(cnaloh.distrib.extinct)]
boxplot(p1[,1],p2[,1],
        p1[,2],p2[,2],
        p1[,3],p2[,3],
        p1[,4],p2[,4],
        p1[,5],p2[,5],
        p1[,6],p2[,6],
        p1[,7],p2[,7],
        ylim=c(-2,25),
        col=c("white","darkgray"),las=2,axes=F,outline=F,outcex=0.5,whisklty="solid",staplewex=1,medlty="solid",medlwd=1,
        at=sort(c(seq(1.1,13.1,by=2),seq(1.9,13.9,by=2))),
        boxwex=0.6)
axis(1,at=seq(0.5,14.5,by=2),lty=1, c(0,expression(10^2),expression(10^3),expression(10^4),
                           expression(10^5),expression(10^6),expression(10^7),expression(10^8)) )
axis(2,at=seq(0,25,by=5),las=1)
box()
mtext(side=2,line=3, "Number of SGA per biopsy per year")
mtext(side=1,line=2, "Lesion size (bp)")


legend("topright", legend=c("Off NSAID","On NSAID","Wilcox p<.05"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,0),lwd=c(0,0,1),pch=c(NA,NA,8),merge=T,box.col="black",cex=0.5,bg="white")

for(i in 1:7){
  pval = wilcox.test(p1[,i],p2[,i],alt="less")$p.val
  if(pval<0.05){
    ## segments(i*2-1,-2,i*2,-2)
    ## segments(i*2-1,-2,i*2-1,-0.5)
    ## segments(i*2,-2,i*2,-0.5)
    points(i*2-0.5,-2,pch=8,cex=0.5)
  }
  ##text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+1,cex=0.8)
  print(pval)
}

dev.off()


## Sep 29

## Extinctions
off.nsaids = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),]
extinct.off.nsaids.total.rates = rowSums(off.nsaids[,3:ncol(off.nsaids)])

on.nsaids = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),]
extinct.on.nsaids.total.rates = rowSums(on.nsaids[,3:ncol(on.nsaids)])

## Novel
off.nsaids = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),]
new.off.nsaids.total.rates = rowSums(off.nsaids[,3:ncol(off.nsaids)])

on.nsaids = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),]
new.on.nsaids.total.rates = rowSums(on.nsaids[,3:ncol(on.nsaids)])



round(median(new.off.nsaids.total.rates),2)
round(quantile(new.off.nsaids.total.rates,probs=c(0.05,0.95)),2)

round(median(new.on.nsaids.total.rates),2)
round(quantile(new.on.nsaids.total.rates,probs=c(0.05,0.95)),2)

round(median(extinct.off.nsaids.total.rates),2)
round(quantile(extinct.on.nsaids.total.rates,probs=c(0.05,0.95)),2)

round(median(extinct.on.nsaids.total.rates),2)
round(quantile(extinct.off.nsaids.total.rates,probs=c(0.05,0.95)),2)

wilcox.test(new.off.nsaids.total.rates,new.on.nsaids.total.rates,alt="greater")
wilcox.test(extinct.off.nsaids.total.rates,extinct.on.nsaids.total.rates,alt="less")

length(new.off.nsaids.total.rates)
length(new.on.nsaids.total.rates)
length(extinct.off.nsaids.total.rates)
length(extinct.on.nsaids.total.rates)


##quantile(new.off.nsaids.total.rates,probs=c(0.05,0.95))












## Chisq test

chisq.test(rbind(p1.mean,p2.mean))

## For Xiaohong
## a = rbind(round(p1,digits=4),round(p2,digits=4))
## colnames(a) = c("Individual_id","On(1)_Off(0)_NSAID","Rate_0-10^2bp","Rate_10^2-10^3bp",
##           "Rate_10^3-10^4bp","Rate_10^4-10^5bp","Rate_10^5-10^6bp","Rate_10^6-10^7bp","Rate_10^7-10^8bp") 
## write.table(a,"p1.csv",row.names=F,quote=F,sep=",")










pdf("kostadinov_fig_extinct_sga_events_aug31.pdf",width=3,height=3,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p0 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),]

p0.mean = apply(p0[,3:ncol(p0)], 2, function(x) mean(x))
p0.sd = apply(p0[,3:ncol(p0)], 2, function(x) sd(x))
p1 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),]

p1.mean = apply(p1[,3:ncol(p1)], 2, function(x) mean(x))
p1.sd = apply(p1[,3:ncol(p1)], 2, function(x) sd(x))
bp = barplot(rbind(p0.mean,p1.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGA per biopsy per year",
        xlab = "Lesion size (bp)", ylim=c(0,70),names.arg=rep("",length(p0.mean)),las=1)

## Plot 95% CI
error.bar(bp,rbind(p0.mean,p1.mean),rbind(1.96*p0.sd/sqrt(nrow(p0)),1.96*p1.sd/sqrt(nrow(p1))),code=1,length=0.02)

## error.bar(bp,rbind(p0.mean,p1.mean),rbind(p0.sd,p1.sd),code=1,length=0.02)

axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))

legend("topright",legend=c("Off NSAID","On NSAID","95% CI"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,1),merge=T,box.col="black",cex=0.8)

for(i in 3:9){
  pval = t.test(p0[,i],p1[,i],alt="less")$p.val
  text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+3,cex=0.8)
  print(pval)
}

chisq.test(rbind(p0.mean,p1.mean))


dev.off()




## Plotting the total
pdf("kostadinov_fig_new_and_extinct_sga_events_aug31.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))

## Extinctions
off.nsaids = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),]
extinct.off.nsaids.total.rates = rowSums(off.nsaids[,3:ncol(off.nsaids)])

on.nsaids = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),]
extinct.on.nsaids.total.rates = rowSums(on.nsaids[,3:ncol(on.nsaids)])

t.test(off.nsaids.total.rates,on.nsaids.total.rates)

## Novel
off.nsaids = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),]
new.off.nsaids.total.rates = rowSums(off.nsaids[,3:ncol(off.nsaids)])

on.nsaids = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),]
new.on.nsaids.total.rates = rowSums(on.nsaids[,3:ncol(on.nsaids)])

t.test(off.nsaids.total.rates,on.nsaids.total.rates)
wilcox.test(off.nsaids.total.rates,on.nsaids.total.rates)


## boxplots
boxplot(new.off.nsaids.total.rates,new.on.nsaids.total.rates, extinct.off)

        


boxplot.stats <-  function (x, coef = NULL, do.conf = TRUE, do.out = TRUE) 
{ 
  nna <- !is.na(x) 
  n <- sum(nna) 
  stats <- quantile(x, c(.05,.25,.5,.75,.95), na.rm = TRUE) 
  iqr <- diff(stats[c(2, 4)]) 
  out <- x < stats[1] | x > stats[5] 
  conf <- if (do.conf) 
    stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n) 
  list(stats = stats, n = n, conf = conf, out = x[out & nna]) 
} 

pdf("kostadinov_fig_new_and_extinct_sga_events_boxplot_aug31.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))

boxplot(new.off.nsaids.total.rates,new.on.nsaids.total.rates,
        extinct.off.nsaids.total.rates, extinct.on.nsaids.total.rates,ylim=c(0,100),
        col=c("white","gray","white","gray"),names=c("New","New","Extinct","Extinct"),las=2,outpch=".")
mtext(side=2,at=50,line=3, "Number of SGA per biopsy per year")
legend("top",legend=c("Off NSAID","On NSAID","95%"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,1),merge=T,box.col=NA,cex=0.8)

pval = wilcox.test(new.on.nsaids.total.rates,new.off.nsaids.total.rates,alt="less")$p.val
text(paste("p=",round(pval,digits=3),sep=""),x=1.5,y=60,cex=0.8)

pval = wilcox.test(extinct.on.nsaids.total.rates,extinct.off.nsaids.total.rates,alt="greater")$p.val
text(paste("p=",round(pval,digits=3),sep=""),x=3.5,y=60,cex=0.8)

dev.off()





pdf("kostadinov_fig_new_and_extinct_sga_events_barplot_aug31.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))

new.mean = c(mean(new.off.nsaids.total.rates),mean(new.on.nsaids.total.rates))
new.sd = c(sd(new.off.nsaids.total.rates),sd(new.on.nsaids.total.rates))
extinct.mean = c(mean(extinct.off.nsaids.total.rates),mean(extinct.on.nsaids.total.rates))
extinct.sd = c(sd(extinct.off.nsaids.total.rates),sd(extinct.on.nsaids.total.rates))

bp = barplot(rbind(new.mean,extinct.mean),beside=T,col=c("white","darkgray"),
        space=c(0,1),
        ylab = "Number of SGA per biopsy per year",
        xlab = NA, ylim=c(0,200),names.arg=c("New SGA","Extinct SGA"),las=1)

## Plot 95% CI
new.off.nsaids.95.CI = 1.96*sd(new.off.nsaids.total.rates)/sqrt(length(new.off.nsaids.total.rates))
new.on.nsaids.95.CI = 1.96*sd(new.on.nsaids.total.rates)/sqrt(length(new.on.nsaids.total.rates))
extinct.off.nsaids.95.CI = 1.96*sd(extinct.off.nsaids.total.rates)/sqrt(length(extinct.off.nsaids.total.rates))
extinct.on.nsaids.95.CI = 1.96*sd(extinct.on.nsaids.total.rates)/sqrt(length(extinct.on.nsaids.total.rates))

error.bar(bp,rbind(new.mean,extinct.mean),rbind(c(new.off.nsaids.95.CI,new.on.nsaids.95.CI),
                                                c(extinct.off.nsaids.95.CI,extinct.on.nsaids.95.CI)),code=1,length=0.02)
legend("topleft",legend=c("Off NSAID","On NSAID","95% CI"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,1),merge=T,box.col=NA,cex=0.8)

pval = t.test(new.on.nsaids.total.rates,new.off.nsaids.total.rates,alt="less")$p.val
text(paste("p=",round(pval,digits=3),sep=""),x=2,y=100,cex=0.8)

pval = t.test(extinct.on.nsaids.total.rates,extinct.off.nsaids.total.rates,alt="greater")$p.val
text(paste("p=",round(pval,digits=3),sep=""),x=5,y=180,cex=0.8)

dev.off()
















pdf("kostadinov_fig_new_sga_events_aug24.pdf",width=3,height=3,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),]
p1.mean = apply(p1[,3:ncol(p1)], 2, function(x) mean(x))
p1.sd = apply(p1[,3:ncol(p1)], 2, function(x) sd(x))
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),]
p2.mean = apply(p2[,3:ncol(p2)], 2, function(x) mean(x))
p2.sd = apply(p2[,3:ncol(p2)], 2, function(x) sd(x))
bp = barplot(rbind(p1.mean,p2.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGAs per biopsy per year",
        xlab = "Lesion size (bp)",names.arg=rep("",length(p1.mean)), ylim=c(0,30),las=1)
## Plot 95% CI
error.bar(bp,rbind(p1.mean,p2.mean),rbind(1.96*p1.sd/sqrt(nrow(p1)),1.96*p2.sd/sqrt(nrow(p2))),code=1,length=0.02)

axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))
legend("topright",legend=c("Off NSAID","On NSAID","95% CI"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,1),merge=T,box.col="black",cex=0.8)

for(i in 3:9){
  pval = t.test(p1[,i],p2[,i])$p.val
  text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+1,cex=0.8)
  print(pval)
}

## Chisq test

chisq.test(rbind(p1.mean,p2.mean))

dev.off()







p1.summary = aggregate(p1, by=list(p1[,"pind"]), FUN=mean)
p2.summary = aggregate(p2, by=list(p2[,"pind"]), FUN=mean)
x2 = cbind(rowSums(p1.summary[,4:ncol(p1.summary)]),rowSums(p2.summary[,4:ncol(p2.summary)]))
pval = t.test(x2[,1],x2[,2])$p.val




CairoSVG(paste("",sep=""),width=2.2,height=2.2,pointsize=8,family="Helvetica")
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
bcols = c(rep("black",2),rep("darkgray",11))

my.barplot(rbind(rev(x2[,1]),rev(x2[,2])),beside=T,horiz=T,
           space=c(0,0),
           col=rbind(c(rep("white",2),rep("darkgray",11)),
             c(rep("darkgray",2),rep("white",11))),           
           ylab = NA, xaxt = "n",
           xlab = NA, xlim=c(0,150),names.arg=rep("",length(p0)),las=1)
axis(1,line=0,at=c(0,seq(5,20,by=5)),tick=T,lab=NA)
mtext(1,line=1,at=c(0,seq(5,20,by=5)), text=c(0,seq(5,20,by=5)))
mtext(2,line=1,at=c(seq(1,13))-0.5, text=c(seq(13,1)),las=1)
mtext(1,line=2,text="Sampling duration (Years)")
mtext(2,line=2,text="Participant ID")



        col=rbind(c(rep("darkgray",2),rep("white",11)),
                  c(rep("white",2),rep("darkgray",11))),


boxplot(x2[,1],horizontal=T,xlim=c(0,13),ylim=c(0,100),yaxs="i",medlty="blank",boxlty=1,staplelty=1,whisklty=1
           ,border=bcols,col=bcols,xlab=NA,ylab=NA,las=1,axes=F,outline=F)
axis(side=1,at=c(-3,-2,-1,0,1,2,3),lab=NA)
mtext(side=1,line=1,at=c(-3,-2,-1,0,1,2,3),
      text=c(expression(10^-3),"0.01","0.1","1","10","100",expression(10^3)),las=1)
axis(side=2,at=c(1:13),lab=NA)
bla = c("a","b","c","d","e","f","g","h","i","j","k","l","m")
bla = rev(bla)
mtext(side=2,line=1,at=c(1:13),text=bla,las=1)
mtext(side=2,line=2,at=7.5,text="Individuals")
mtext(side=1,line=2,at=0,text="SGA rate")
par(new=T)
bcols = c(rep("darkgray",2),rep("black",11))
boxplot(rev(x3),horizontal=T,ylim=c(-3,3),yaxs="i",medlty="blank",boxlty=1,staplelty=1,whisklty=1
           ,col=bcols,border=bcols,xlab=NA,ylab=NA,las=1,axes=F,outline=F)
legend("topright",legend=c("On NSAID","Off NSAID"),fill=c("darkgray","black"),box.col=NA,cex=0.7)
box()
dev.off()








pdf("kostadinov_fig_extinct_sga_events_aug16.pdf",width=3,height=3,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p0 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),]

p0.mean = apply(p0[,3:ncol(p0)], 2, function(x) mean(x))
p0.sd = apply(p0[,3:ncol(p0)], 2, function(x) sd(x))
p1 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),]

p1.mean = apply(p1[,3:ncol(p1)], 2, function(x) mean(x))
p1.sd = apply(p1[,3:ncol(p1)], 2, function(x) sd(x))
bp = barplot(rbind(p0.mean,p1.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGA per biopsy per year",
        xlab = "Lesion size (bp)", ylim=c(0,50),names.arg=rep("",length(p0.mean)),las=1)

## Plot 95% CI
error.bar(bp,rbind(p0.mean,p1.mean),rbind(1.96*p0.sd/sqrt(nrow(p0)),1.96*p1.sd/sqrt(nrow(p1))),code=1,length=0.02)

## error.bar(bp,rbind(p0.mean,p1.mean),rbind(p0.sd,p1.sd),code=1,length=0.02)

axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))

legend("topright",legend=c("Off NSAID","On NSAID","95% CI"),
       fill=c("white","darkgray",NA),
       border=c("black","black",NA),
       lty=c(0,0,1),merge=T,box.col="black",cex=0.8)

for(i in 3:9){
  pval = t.test(p0[,i],p1[,i])$p.val
  text(paste("p=",round(pval,digits=2),sep=""),x=(i-3)*4+2.5,y=p1.mean[i-2]+1.96*p1.sd[i-2]/sqrt(nrow(p1))+3,cex=0.8)
  print(pval)
}

chisq.test(rbind(p0.mean,p1.mean))


dev.off()
















## This is the old strange time normalization

## Plot new events

pdf("kostadinov_fig_new_sga_events_aug16.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),]

time.off.nsaids = patient.off.nsaids.time[p1[,"pind"]]
p1[,c(3:9)] = p1[,c(3:9)]/time.off.nsaids

p1.mean = apply(p1[,3:ncol(p1)], 2, function(x) mean(x))
p1.sd = apply(p1[,3:ncol(p1)], 2, function(x) sd(x))
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),]
time.on.nsaids = patient.on.nsaids.time[p2[,"pind"]]
p2[,c(3:9)] = p2[,c(3:9)]/time.on.nsaids

p2.mean = apply(p2[,3:ncol(p2)], 2, function(x) mean(x))
p2.sd = apply(p2[,3:ncol(p2)], 2, function(x) sd(x))
bp = barplot(rbind(p1.mean,p2.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGAs per biopsy per year",
        xlab = "Lesion size (bp)",names.arg=rep("",length(p1.mean)), ylim=c(0,20),las=1)
error.bar(bp,rbind(p1.mean,p2.mean),rbind(p1.sd,p2.sd),code=1,length=0.02)
axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))
legend("topright",legend=c("Off NSAID","On NSAID","p<0.05","p<0.01"),
       fill=c("white","darkgray",NA,NA),
       border=c("black","black",NA,NA),
       pch=c(NA,NA,4,8),box.col="black",cex=0.8)

for(i in 3:9){
  pval = wilcox.test(p1[,i],p2[,i])$p.val
  ## if(pval<0.01){
  ##   text((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,lab="***")
  ## } else
  if(pval<0.01){
    points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=8)
  } else if(pval<0.05){
    points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=4)
  }
  print(pval)
}
## [1] 0.595891
## [1] 0.1137371
## [1] 0.02077306
## [1] 0.5665799
## [1] 0.6886865
## [1] 0.06110919
## [1] 0.4049137

dev.off()

## Plot extinctions
pdf("kostadinov_fig_extinct_sga_events_aug16.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p0 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),]
time.off.nsaids = patient.off.nsaids.time[p0[,"pind"]]
p0[,c(3:9)] = p0[,c(3:9)]/time.off.nsaids

p0.mean = apply(p0[,3:ncol(p0)], 2, function(x) mean(x))
p0.sd = apply(p0[,3:ncol(p0)], 2, function(x) sd(x))
p1 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),]
time.on.nsaids = patient.on.nsaids.time[p1[,"pind"]]
p1[,c(3:9)] = p1[,c(3:9)]/time.on.nsaids

p1.mean = apply(p1[,3:ncol(p1)], 2, function(x) mean(x))
p1.sd = apply(p1[,3:ncol(p1)], 2, function(x) sd(x))
bp = barplot(rbind(p0.mean,p1.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGA per biopsy per year",
        xlab = "Lesion size (bp)", ylim=c(0,20),names.arg=rep("",length(p0.mean)),las=1)
error.bar(bp,rbind(p0.mean,p1.mean),rbind(p0.sd,p1.sd),code=1,length=0.02)
axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))

legend("topright",legend=c("Off NSAID","On NSAID","p<0.05","p<0.01"),
       fill=c("white","darkgray",NA,NA),
       border=c("black","black",NA,NA),
       pch=c(NA,NA,4,8),box.col="black",cex=0.8)

for(i in 3:9){
  pval = wilcox.test(p1[,i],p2[,i])$p.val
  ## if(pval<0.01){
  ##   text((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,lab="***")
  ## } else
  if(pval<0.01){
    points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=8)
  } else if(pval<0.05){
    points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=4)
  }
  print(pval)
}
## [1] 0.1283733
## [1] 0.4018579
## [1] 0.553781
## [1] 0.6609022
## [1] 0.9014483
## [1] 0.970036
## [1] 0.9929115

dev.off()







## Sum across all leasion sizes




pdf("kostadinov_fig_new_sga_events_any_lesion_size_aug16.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p1 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==1) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==2),]

p1 = cbind(p1[,c(1,2)],rowSums(p1[,3:9]))

time.off.nsaids = patient.off.nsaids.time[p1[,"pind"]]
p1[,3] = p1[,3]/time.off.nsaids

p1.mean = mean(p1[,3])
p1.sd = sd(p1[,3])
p2 = cnaloh.distrib.new[(cnaloh.distrib.new[,1]%in%nc.users & cnaloh.distrib.new[,2]==2) |
                        (cnaloh.distrib.new[,1]%in%cf.users & cnaloh.distrib.new[,2]==1),]

p2 = cbind(p2[,c(1,2)],rowSums(p2[,3:9]))
time.on.nsaids = patient.on.nsaids.time[p2[,"pind"]]
p2[,3] = p2[,3]/time.on.nsaids

p2.mean = mean(p2[,3])
p2.sd = sd(p2[,3])
bp = barplot(rbind(p1.mean,p2.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGAs per biopsy per year",
        xlab = "Lesion size (bp)",names.arg=rep("",length(p1.mean)), ylim=c(0,50),las=1)
error.bar(bp,rbind(p1.mean,p2.mean),rbind(p1.sd,p2.sd),code=1,length=0.02)
axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))
legend("topright",legend=c("Off NSAID","On NSAID","p<0.05","p<0.01"),
       fill=c("white","darkgray",NA,NA),
       border=c("black","black",NA,NA),
       pch=c(NA,NA,4,8),box.col="black",cex=0.8)

i=3
pval = wilcox.test(p1[,i],p2[,i])$p.val
if(pval<0.01){
  points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=8)
} else if(pval<0.05){
  points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=4)
}
print(pval)
## [1] 0.2601233

dev.off()

## Plot extinctions
pdf("kostadinov_fig_extinct_sga_events_any_lesion_size_aug16.pdf",width=2.5,height=2.5,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
p0 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==0) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==1),]

p0 = cbind(p0[,c(1,2)],rowSums(p0[,3:9]))

time.off.nsaids = patient.off.nsaids.time[p0[,"pind"]]
p0[,3] = p0[,3]/time.off.nsaids

p0.mean = mean(p0[,3])
p0.sd = sd(p0[,3])
p1 = cnaloh.distrib.extinct[(cnaloh.distrib.extinct[,1]%in%nc.users & cnaloh.distrib.extinct[,2]==1) |
                        (cnaloh.distrib.extinct[,1]%in%cf.users & cnaloh.distrib.extinct[,2]==0),]
p1 = cbind(p1[,c(1,2)],rowSums(p1[,3:9]))
time.on.nsaids = patient.on.nsaids.time[p1[,"pind"]]
p1[,3] = p1[,3]/time.on.nsaids

p1.mean = mean(p1[,3])
p1.sd = sd(p1[,3])
bp = barplot(rbind(p0.mean,p1.mean),beside=T,col=c("white","darkgray"),
        space=c(0,2),
        ylab = "Number of SGA per biopsy per year",
        xlab = "Lesion size (bp)", ylim=c(0,50),names.arg=rep("",length(p0.mean)),las=1)
error.bar(bp,rbind(p0.mean,p1.mean),rbind(p0.sd,p1.sd),code=1,length=0.02)
axis(1,at=seq(0,28,by=4)+1,lab=c(0,expression(10^2),expression(10^3),expression(10^4),
                          expression(10^5),expression(10^6),expression(10^7),expression(10^8)))

legend("topright",legend=c("Off NSAID","On NSAID","p<0.05","p<0.01"),
       fill=c("white","darkgray",NA,NA),
       border=c("black","black",NA,NA),
       pch=c(NA,NA,4,8),box.col="black",cex=0.8)

i=3
pval = wilcox.test(p1[,i],p2[,i])$p.val
if(pval<0.01){
  points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=8)
} else if(pval<0.05){
  points((i-3)*4+2.5,p1.mean[i-2]+p1.sd[i-2]+5,pch=4)
}
print(pval)
##[1] 0.7911821

dev.off()




























  my.barplot = function (height, width = 1, space = NULL, names.arg = NULL, 
    legend.text = NULL, beside = FALSE, horiz = FALSE, density = NULL, 
    angle = 45, col = NULL, border = par("fg"), main = NULL, 
    sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    xpd = TRUE, log = "", axes = TRUE, axisnames = TRUE, cex.axis = par("cex.axis"), 
    cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
    axis.lty = 0, offset = 0, add = FALSE, args.legend = NULL, 
    ...) 
{
    if (!missing(inside)) 
        .NotYetUsed("inside", error = FALSE)
    if (is.null(space)) 
        space <- if (is.matrix(height) && beside) 
            c(0, 1)
        else 0.2
    space <- space * mean(width)
    if (plot && axisnames && is.null(names.arg)) 
        names.arg <- if (is.matrix(height)) 
            colnames(height)
        else names(height)
    if (is.vector(height) || (is.array(height) && (length(dim(height)) == 
        1))) {
        height <- cbind(height)
        beside <- TRUE
        if (is.null(col)) 
            col <- "grey"
    }
    else if (is.matrix(height)) {
        if (is.null(col)) 
            col <- grey.colors(nrow(height))
    }
    else stop("'height' must be a vector or a matrix")
    if (is.logical(legend.text)) 
        legend.text <- if (legend.text && is.matrix(height)) 
            rownames(height)
    stopifnot(is.character(log))
    logx <- logy <- FALSE
    if (log != "") {
        logx <- length(grep("x", log)) > 0L
        logy <- length(grep("y", log)) > 0L
    }
    if ((logx || logy) && !is.null(density)) 
        stop("Cannot use shading lines in bars when log scale is used")
    NR <- nrow(height)
    NC <- ncol(height)
    if (beside) {
        if (length(space) == 2) 
            space <- rep.int(c(space[2L], rep.int(space[1L], 
                NR - 1)), NC)
        width <- rep(width, length.out = NR)
    }
    else {
        width <- rep(width, length.out = NC)
    }
    offset <- rep(as.vector(offset), length.out = length(width))
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    log.dat <- (logx && horiz) || (logy && !horiz)
    if (log.dat) {
        if (min(height + offset, na.rm = TRUE) <= 0) 
            stop("log scale error: at least one 'height + offset' value <= 0")
        if (logx && !is.null(xlim) && min(xlim) <= 0) 
            stop("log scale error: 'xlim' <= 0")
        if (logy && !is.null(ylim) && min(ylim) <= 0) 
            stop("log scale error: 'ylim' <= 0")
        rectbase <- if (logy && !horiz && !is.null(ylim)) 
            ylim[1L]
        else if (logx && horiz && !is.null(xlim)) 
            xlim[1L]
        else 0.9 * min(height, na.rm = TRUE)
    }
    else rectbase <- 0
    if (!beside) 
        height <- rbind(rectbase, apply(height, 2L, cumsum))
    rAdj <- offset + (if (log.dat) 
        0.9 * height
    else -0.01 * height)
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    if (horiz) {
        if (is.null(xlim)) 
            xlim <- range(rAdj, height + offset, na.rm = TRUE)
        if (is.null(ylim)) 
            ylim <- c(min(w.l), max(w.r))
    }
    else {
        if (is.null(xlim)) 
            xlim <- c(min(w.l), max(w.r))
        if (is.null(ylim)) 
            ylim <- range(rAdj, height + offset, na.rm = TRUE)
    }
    if (beside) 
        w.m <- matrix(w.m, ncol = NC)
    if (plot) {
        opar <- if (horiz) 
            par(xaxs = "i", xpd = xpd)
        else par(yaxs = "i", xpd = xpd)
        on.exit(par(opar))
        if (!add) {
            plot.new()
            plot.window(xlim, ylim, log = log, ...)
        }
        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
            ...) {
            if (horizontal) 
                rect(x1, y1, x2, y2, ...)
            else rect(y1, x1, y2, x2, ...)
        }
        if (beside) 
            xyrect(rectbase + offset, w.l, c(height) + offset, 
                w.r, horizontal = horiz, angle = angle, density = density, 
                col = col, border = border)
        else {
          for (i in 1:NC) {
            for (j in 1:NR) {
              mycol = col[,i]
              xyrect(height[1:NR, i], w.l[i], height[-1, i],
                     w.r[i], horizontal = horiz, col = mycol)
            }
          }
            ## for (i in 1L:NC) {
            ##     xyrect(height[1L:NR, i] + offset[i], w.l[i], 
            ##       height[-1, i] + offset[i], w.r[i], horizontal = horiz, 
            ##       angle = angle, density = density, col = col, 
            ##       border = border)
            ## }
        }
        if (axisnames && !is.null(names.arg)) {
            at.l <- if (length(names.arg) != length(w.m)) {
                if (length(names.arg) == NC) 
                  colMeans(w.m)
                else stop("incorrect number of names")
            }
            else w.m
            axis(if (horiz) 
                2
            else 1, at = at.l, labels = names.arg, lty = axis.lty, 
                cex.axis = cex.names, ...)
        }
        if (!is.null(legend.text)) {
            legend.col <- rep(col, length.out = length(legend.text))
            if ((horiz & beside) || (!horiz & !beside)) {
                legend.text <- rev(legend.text)
                legend.col <- rev(legend.col)
                density <- rev(density)
                angle <- rev(angle)
            }
            xy <- par("usr")
            if (is.null(args.legend)) {
                legend(xy[2L] - xinch(0.1), xy[4L] - yinch(0.1), 
                  legend = legend.text, angle = angle, density = density, 
                  fill = legend.col, xjust = 1, yjust = 1)
            }
            else {
                args.legend1 <- list(x = xy[2L] - xinch(0.1), 
                  y = xy[4L] - yinch(0.1), legend = legend.text, 
                  angle = angle, density = density, fill = legend.col, 
                  xjust = 1, yjust = 1)
                args.legend1[names(args.legend)] <- args.legend
                do.call("legend", args.legend1)
            }
        }
        title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
        if (axes) 
            axis(if (horiz) 
                1
            else 2, cex.axis = cex.axis, ...)
        invisible(w.m)
    }
    else w.m
}





  ## Separate baseline, 1st period, 2nd period
  samples.period.0 = which(samples[,"endo_age_from_baseline"]==0)
  samples.period.1 = which(samples[,"endo_age_from_baseline"]<samples[,"nsaids_transition_age_from_baseline"] &
                           samples[,"endo_age_from_baseline"]>0)
  samples.period.2 = which(samples[,"endo_age_from_baseline"]>samples[,"nsaids_transition_age_from_baseline"])
  ## Get summary statistics  
  for(i in 4:ncol(feat.binary.states)){
    if((i-3)%in%samples.period.0){sid=0}
    if((i-3)%in%samples.period.1){sid=1}
    if((i-3)%in%samples.period.2){sid=2}
    cnaloh.distrib.all = rbind(cnaloh.distrib.all,cbind(pind,sid,t(table(cut(log10(
      feat.binary.states[feat.binary.states[,i]==1,3]-feat.binary.states[feat.binary.states[,i]==1,2]),
      breaks = c(-Inf,seq(2,7,by=1),Inf),right=T)))))
  }


  ## Separate baseline, 1st period, 2nd period, final
  samples.period.0 = which(samples[,"endo_age_from_baseline"]==0)
  samples.period.1 = which(samples[,"endo_age_from_baseline"]<samples[,"nsaids_transition_age_from_baseline"] &
                           samples[,"endo_age_from_baseline"]>0)
  samples.period.2 = which(samples[,"endo_age_from_baseline"]>samples[,"nsaids_transition_age_from_baseline"] &
                           samples[,"endo_coalescent_age"]>0)
  samples.period.3 = which(samples[,"endo_coalescent_age"]==0)  
  ## Get summary statistics  
  for(i in 4:ncol(feat.binary.states)){
    if((i-3)%in%samples.period.0){sid=0}
    if((i-3)%in%samples.period.1){sid=1}
    if((i-3)%in%samples.period.2){sid=2}
    if((i-3)%in%samples.period.3){sid=3}
    cnaloh.distrib.final = rbind(cnaloh.distrib.final,cbind(pind,sid,t(table(cut(log10(
      feat.binary.states[feat.binary.states[,i]==1,3]-feat.binary.states[feat.binary.states[,i]==1,2]),
      breaks = c(-Inf,seq(2,7,by=1),Inf),right=T)))))
  }
