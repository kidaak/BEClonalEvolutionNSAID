
## Connect to mySQL database
library(RMySQL)
passwd = as.character(read.table("/etc/keys/.mysqlpass",as.is=T)[1,1])
con = dbConnect(MySQL(), user="root", password=passwd, dbname="sbep", host="localhost")

library(MASS)
library(RColorBrewer)
k = 9
my.cols <- rev(brewer.pal(k, "Reds"))

## Patient 662, 2 blood samples
bld.662.1 = 258
bld.662.2 = 261

## Patient 437
bld.437.1 = 168
bld.437.2 = 332


d1 = dbGetQuery(con, paste("select snp_id,r,baf from data",
 	   		   "where sample_id=",bld.662.1," order by snp_id"))

d2 = dbGetQuery(con, paste("select snp_id,r,baf from data",
 	   		   "where sample_id=",bld.662.2," order by snp_id"))

d3 = dbGetQuery(con, paste("select snp_id,r,baf from data",
 	   		   "where sample_id=",bld.437.1," order by snp_id"))

d4 = dbGetQuery(con, paste("select snp_id,r,baf from data",
 	   		   "where sample_id=",bld.437.2," order by snp_id"))


cor(d1[,"r"],d2[,"r"])
cor(d3[,"r"],d4[,"r"])

cor(d1[,"r"],d3[,"r"])
cor(d1[,"r"],d4[,"r"])
cor(d2[,"r"],d3[,"r"])
cor(d2[,"r"],d4[,"r"])


cor(d1[,"baf"],d2[,"baf"])
cor(d3[,"baf"],d4[,"baf"])

cor(d1[,"baf"],d3[,"baf"])
cor(d1[,"baf"],d4[,"baf"])
cor(d2[,"baf"],d3[,"baf"])
cor(d2[,"baf"],d4[,"baf"])

##
{
bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparisons_individual_m_baf.png",sep=""),
       width=250,height=250,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)
par(mfrow=c(1,1),mar=c(3,3,2,1))  
cor=cor(d1[,"baf"],d2[,"baf"])
plot(d1[,"baf"], d2[,"baf"],
     main=NA,xlab=NA,ylab=NA,pty="s",xlim=c(0,1),ylim=c(0,1),pch=".")
mtext(paste("Individual m blood #1 vs. blood #2",sep=""),side=3,line=1)
mtext(paste("Pearson r=",round(cor,digits=4),sep=""),side=3,line=0)
mtext(side=1,line=2,text="Blood #1 BAF")
mtext(side=2,line=2,text="Blood #2 BAF")
dev.off()

bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparisons_individual_f_baf.png",sep=""),
       width=250,height=250,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)
par(mfrow=c(1,1),mar=c(3,3,2,1))  
cor=cor(d3[,"baf"],d4[,"baf"])
plot(d3[,"baf"], d4[,"baf"],
     main=NA,xlab=NA,ylab=NA,pty="s",xlim=c(0,1),ylim=c(0,1),pch=".")
mtext(paste("Individual f blood #1 vs. blood #2",sep=""),side=3,line=1)
mtext(paste("Pearson r=",round(cor,digits=4),sep=""),side=3,line=0)
mtext(side=1,line=2,text="Blood #1 BAF")
mtext(side=2,line=2,text="Blood #2 BAF")
dev.off()

bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparisons_individual_mf_baf.png",sep=""),
       width=250,height=250,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)
par(mfrow=c(1,1),mar=c(3,3,2,1))  
cor=cor(d1[,"baf"],d3[,"baf"])
plot(d1[,"baf"], d3[,"baf"],
     main=NA,xlab=NA,ylab=NA,pty="s",xlim=c(0,1),ylim=c(0,1),pch=".")
mtext(paste("Individual m blood #1 vs. Individual f blood #1",sep=""),side=3,line=1)
mtext(paste("Pearson r=",round(cor,digits=4),sep=""),side=3,line=0)
mtext(side=1,line=2,text="Blood m #1 BAF")
mtext(side=2,line=2,text="Blood f #1 BAF")
dev.off()

}



## Compute r, baf correlations

pdf(paste("kostadinov_sga_events_distance_vs_spacetime_gej_normalized.pdf",sep=""),pointsize=10,width=8.5,height=11)
par(mfrow=c(1,2))

plot(d1[,"r"], d2[,"r"], title="Individual m",xlab=NA,ylab=NA,pty="s",xlim=c(0,5),ylim=c(0,5),pch=".")
z <- kde2d(d1[,"r"], d2[,"r"], n=20)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

plot(d1[,"r"], d3[,"r"], title="Individual m",xlab=NA,ylab=NA,pty="s",xlim=c(0,5),ylim=c(0,5),pch=".")
z <- kde2d(d1[,"r"], d3[,"r"], n=25)
contour(z, drawlabels=FALSE, nlevels=k-1, col=my.cols, add=TRUE)

patients = t(dbGetQuery(con, "select distinct(patient_id) from acs_paired_samples"))
pnum = c(1,12,2,3,4,5,6,7,8,9,13,10,11)
library(Cairo)




## Chromosome array
chrarr = c(1:24)
total.base.pairs = 0
for(chri in 1:length(chrarr)){
  chr = chrarr[chri]
  chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
  snp_id.start = chr.start.p.q.stop[,"start"] 
  snp_id.stop = chr.start.p.q.stop[,"stop"]
  snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
    "where snp_id>=",snp_id.start,"and",
    "snp_id <=",snp_id.stop,
    "order by position"))
  snp_pos.start = min(snp_pos[,"position"])
  snp_pos.stop = max(snp_pos[,"position"])
  total.base.pairs = total.base.pairs + snp_pos.stop - snp_pos.start +1 
}




##for (i in 1:length(sample_ids)){
##  sample_id    = sample_ids[i]
  ## reference_id = reference_ids[i]
  ## patient_id = patient_ids[i]
  ## date_level_id = date_level_ids[i]

img.width = 750
img.height = 180
##chr = -1 ## plot the entire genome
ylim.lower.log2r = -2.0
ylim.upper.log2r = 2.0
ylim.lower.baf = 0.0
ylim.upper.baf = 1.0
PSEUDO = 0.000001
het.upper = 0.66
het.lower = 0.33
y.offset.log2r = 0
y.offset.baf = 0

{
  bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparisons_individual_m_logr.png",sep=""),
         width=img.width,height=img.height,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)

  par(mfrow=c(1,1),mar=c(2,4,2,1))  

  xaxsticks = NULL
  current.position = 0
  ## Loop all chromosomes
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0

    data = d1
    data.ref = d2
    log2r = cbind(data[,"snp_id"],log2((data[,"r"]+PSEUDO)/(data.ref[,"r"]+PSEUDO)))
    colnames(log2r) = c("snp_id","log2r")
    ## Exclude NA SNPs
    log2r = log2r[which(!is.na(log2r[,"log2r"])),]
    ## Splush log2r values exceeding the y limits
    log2r[which(log2r[,"log2r"] < ylim.lower.log2r), "log2r"] = ylim.lower.log2r
    log2r[which(log2r[,"log2r"] > ylim.upper.log2r), "log2r"] = ylim.upper.log2r  

    if(chr==1){
      cor = cor(data[,"r"],data.ref[,"r"])
      plot(0,0,ylim=c(ylim.lower.log2r-y.offset.log2r,ylim.upper.log2r+y.offset.log2r),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i")
      axis(2,at=c(ylim.lower.log2r,0,-1,1,ylim.upper.log2r), labels=T, tick=T,las=2)
      title(ylab=bquote(Log[2] ~ R),line=2)
      mtext(paste("Individual m blood #1 vs. Individual m blood #2 comparison, Pearson r=",round(cor,digits=4),sep=""),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }

  cp = current.position - snp_pos.start

  points(cp+snp_pos[match(log2r[,"snp_id"],snp_pos[,"snp_id"]),"position"],
         log2r[,"log2r"], pch=".", col = "black")
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")
  dev.off()
}

{
  bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparisons_individual_f_logr.png",sep=""),
         width=img.width,height=img.height,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)

  par(mfrow=c(1,1),mar=c(2,4,2,1))  

  xaxsticks = NULL
  current.position = 0
  ## Loop all chromosomes
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0

    data = d3
    data.ref = d4
    log2r = cbind(data[,"snp_id"],log2((data[,"r"]+PSEUDO)/(data.ref[,"r"]+PSEUDO)))
    colnames(log2r) = c("snp_id","log2r")
    ## Exclude NA SNPs
    log2r = log2r[which(!is.na(log2r[,"log2r"])),]
    ## Splush log2r values exceeding the y limits
    log2r[which(log2r[,"log2r"] < ylim.lower.log2r), "log2r"] = ylim.lower.log2r
    log2r[which(log2r[,"log2r"] > ylim.upper.log2r), "log2r"] = ylim.upper.log2r  

    if(chr==1){
      cor = cor(data[,"r"],data.ref[,"r"])
      plot(0,0,ylim=c(ylim.lower.log2r-y.offset.log2r,ylim.upper.log2r+y.offset.log2r),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i")
      axis(2,at=c(ylim.lower.log2r,0,-1,1,ylim.upper.log2r), labels=T, tick=T,las=2)
      title(ylab=bquote(Log[2] ~ R),line=2)
      mtext(paste("Individual f blood #1 vs. Individual f blood #2 comparison, Pearson r=",round(cor,digits=4),sep=""),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }

  cp = current.position - snp_pos.start

  points(cp+snp_pos[match(log2r[,"snp_id"],snp_pos[,"snp_id"]),"position"],
         log2r[,"log2r"], pch=".", col = "black")
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")
  dev.off()
}

{
  bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparisons_individual_mf_logr.png",sep=""),
         width=img.width,height=img.height,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)

  par(mfrow=c(1,1),mar=c(2,4,2,1))  

  xaxsticks = NULL
  current.position = 0
  ## Loop all chromosomes
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0

    data = d1
    data.ref = d3
    log2r = cbind(data[,"snp_id"],log2((data[,"r"]+PSEUDO)/(data.ref[,"r"]+PSEUDO)))
    colnames(log2r) = c("snp_id","log2r")
    ## Exclude NA SNPs
    log2r = log2r[which(!is.na(log2r[,"log2r"])),]
    ## Splush log2r values exceeding the y limits
    log2r[which(log2r[,"log2r"] < ylim.lower.log2r), "log2r"] = ylim.lower.log2r
    log2r[which(log2r[,"log2r"] > ylim.upper.log2r), "log2r"] = ylim.upper.log2r  

    if(chr==1){
      cor = cor(data[,"r"],data.ref[,"r"])
      plot(0,0,ylim=c(ylim.lower.log2r-y.offset.log2r,ylim.upper.log2r+y.offset.log2r),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i")
      axis(2,at=c(ylim.lower.log2r,0,-1,1,ylim.upper.log2r), labels=T, tick=T,las=2)
      title(ylab=bquote(Log[2] ~ R),line=2)
      mtext(paste("Individual m blood #1 vs. Individual f blood #1 comparison, Pearson r=",round(cor,digits=4),sep=""),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }

  cp = current.position - snp_pos.start

  points(cp+snp_pos[match(log2r[,"snp_id"],snp_pos[,"snp_id"]),"position"],
         log2r[,"log2r"], pch=".", col = "black")
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")
  dev.off()
}






  bitmap(paste("supplemental_figure_blood_samples_technical_noise_comparison_f.png",sep=""),
         width=img.width,height=img.height,units="px",pointsize=10,type="pngalpha",bg="white",taa=4)

  par(mfrow=c(5,1),mar=c(2,4,2,1))  

  xaxsticks = NULL
  current.position = 0
  ## Loop all chromosomes
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0

    data = d1
    data.ref = d2
    log2r = cbind(data[,"snp_id"],log2((data[,"r"]+PSEUDO)/(data.ref[,"r"]+PSEUDO)))
    colnames(log2r) = c("snp_id","log2r")
    ## Exclude NA SNPs
    log2r = log2r[which(!is.na(log2r[,"log2r"])),]
    ## Splush log2r values exceeding the y limits
    log2r[which(log2r[,"log2r"] < ylim.lower.log2r), "log2r"] = ylim.lower.log2r
    log2r[which(log2r[,"log2r"] > ylim.upper.log2r), "log2r"] = ylim.upper.log2r  

    if(chr==1){
      plot(0,0,ylim=c(ylim.lower.log2r-y.offset.log2r,ylim.upper.log2r+y.offset.log2r),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i")
      axis(2,at=c(ylim.lower.log2r,0,-1,1,ylim.upper.log2r), labels=T, tick=T,las=2)
      title(ylab=bquote(Log[2] ~ R),line=2)
      mtext(paste("Individual m, Blood vs. blood comparison"),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }

  cp = current.position - snp_pos.start

  log2r.hets = log2r[match(hets,log2r[,"snp_id"]),]
  log2r.homs = log2r[-match(hets,log2r[,"snp_id"]),]    
  points(cp+snp_pos[match(log2r.homs[,"snp_id"],snp_pos[,"snp_id"]),"position"],
         log2r.homs[,"log2r"], pch=".", col = "black")
  points(cp+snp_pos[match(log2r.hets[,"snp_id"],snp_pos[,"snp_id"]),"position"],
         log2r.hets[,"log2r"], pch=".", col = "black")
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")



 ## Plot the profile only and annotations
    current.position = 0
  ## Loop all chromosomes
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0
          
    
    annot=dbGetQuery(con,paste("select * from sample_features s",
      "where s.sample_id=",sample_id," and s.reference_id=",reference_id,
      " and s.stop_snp_id >=",snp_id.start," and s.start_snp_id <=",snp_id.stop))
    
    annot[which(annot[,"hom_logr"] < ylim.lower.log2r), "hom_logr"] = ylim.lower.log2r
    annot[which(annot[,"hom_logr"] > ylim.upper.log2r), "hom_logr"] = ylim.upper.log2r
    annot[which(annot[,"het_logr"] < ylim.lower.log2r), "het_logr"] = ylim.lower.log2r
    annot[which(annot[,"het_logr"] > ylim.upper.log2r), "het_logr"] = ylim.upper.log2r

    if(chr==1){
      plot(0,0,ylim=c(ylim.lower.log2r-y.offset.log2r,ylim.upper.log2r+y.offset.log2r),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i")
      axis(2,at=c(ylim.lower.log2r,0,-1,1,ylim.upper.log2r), labels=T, tick=T,las=2)
      title(ylab=bquote(Log[2] ~ R),line=2)
      mtext(paste("Individual",patient_id,date_level_id),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }


  cp = current.position - snp_pos.start


  segments(cp+annot$start_position, annot$het_logr, cp+annot$stop_position, annot$het_logr, col="red")
  segments(cp+annot$start_position, annot$hom_logr, cp+annot$stop_position, annot$hom_logr, col="blue")

    ## profile = rbind(cbind(cp+annot$start_position, annot$hom_logr),
    ##   cbind(cp+annot$stop_position, annot$hom_logr))
    ## profile = profile[order(profile[,1]),]
    ## lines(profile[,1], profile[,2], col="blue")

    ## profile = rbind(cbind(cp+annot$start_position, annot$het_logr),
    ##   cbind(cp+annot$stop_position, annot$het_logr))
    ## profile = profile[order(profile[,1]),]
    ## lines(profile[,1], profile[,2], col="red")

    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")


  current.position = 0
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0
    
    
    ##Check if data is slotted into memory
    check = dbGetQuery(con,paste("select snp_id from data_memory 
                           where sample_id=",sample_id,
      " and snp_id =",snp_id.start))
    if(length(check)>1){
      data=dbGetQuery(con,paste("select snp_id, baf from data_memory 
                           where sample_id=",sample_id,
        " and snp_id >=",snp_id.start," and snp_id <=",snp_id.stop,
        " order by snp_id"))
      
      if(reference_id != prev_ref_id){
        
        data.ref=dbGetQuery(con,paste("select snp_id,  baf from data_memory
                           where sample_id=",reference_id,
          " and snp_id >=",snp_id.start," and snp_id <=",snp_id.stop,
          " order by snp_id"))
        hets = data.ref[data.ref["baf"]<=het.upper & data.ref["baf"]>=het.lower,"snp_id"]
      }
      prev_ref_id = reference_id
      
    } else {
      data=dbGetQuery(con,paste("select snp_id,  baf from data 
                           where sample_id=",sample_id,
        " and snp_id >=",snp_id.start," and snp_id <=",snp_id.stop,
        " order by snp_id"))
      
      if(reference_id != prev_ref_id){
        
        data.ref=dbGetQuery(con,paste("select snp_id,  baf from data
                           where sample_id=",reference_id,
          " and snp_id >=",snp_id.start," and snp_id <=",snp_id.stop,
          " order by snp_id"))
        hets = data.ref[data.ref["baf"]<=het.upper & data.ref["baf"]>=het.lower,"snp_id"]
      }
      prev_ref_id = reference_id
      
      
    }

    baf = cbind(data[,"snp_id"],data[,"baf"])
    colnames(baf) = c("snp_id","baf")
    # Exclude NA SNPs
    baf = baf[which(!is.na(baf[,"baf"])),]
    ## Normalize
    CORFACTOR = 0.53
    baf[baf[,"baf"]<=CORFACTOR,"baf"]=baf[baf[,"baf"]<=CORFACTOR,"baf"]*0.5/CORFACTOR
    baf[baf[,"baf"]>CORFACTOR,"baf"]=0.5+0.5*(baf[baf[,"baf"]>CORFACTOR,"baf"]-CORFACTOR)/(1-CORFACTOR)
    if(chr==1){
      plot(0,0,ylim=c(ylim.lower.baf-y.offset.baf,ylim.upper.baf+y.offset.baf),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i",pch=NA)
      axis(2,at=c(0,0.5,1), labels=c("0","0.5","1"), tick=T,las=2)
      title(ylab="mBAF",line=2)
      mtext(paste("Individual",patient_id,date_level_id),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }
  cp = current.position - snp_pos.start 
    baf.hets = baf[match(hets,baf[,"snp_id"]),]
    baf.homs = baf[-match(hets,baf[,"snp_id"]),]    
    points(cp+snp_pos[match(baf.homs[,"snp_id"],snp_pos[,"snp_id"]),"position"],
           baf.homs[,"baf"], pch=".", col = homs.col)
    points(cp+snp_pos[match(baf.hets[,"snp_id"],snp_pos[,"snp_id"]),"position"],
           baf.hets[,"baf"], pch=".", col = "black")
    ##segments(cp+annot$start_position, annot$mbaf, cp+annot$stop_position, annot$mbaf, col="red")

    ## profile = rbind(cbind(cp+annot$start_position, annot$mbaf),
    ##                 cbind(cp+annot$stop_position, annot$mbaf))
    ## profile = profile[order(profile[,1]),]
    ## lines(profile[,1], profile[,2], col="red")
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")

  ## Plot profile + annotations
  current.position = 0
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0
    
    
    annot=dbGetQuery(con,paste("select * from sample_features s",
      "where s.sample_id=",sample_id," and s.reference_id=",reference_id,
      " and s.stop_snp_id >=",snp_id.start," and s.start_snp_id <=",snp_id.stop))

    

    if(chr==1){
      plot(0,0,ylim=c(ylim.lower.baf-y.offset.baf,ylim.upper.baf+y.offset.baf),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i",pch=NA)
      axis(2,at=c(0,0.5,1), labels=c("0","0.5","1"), tick=T,las=2)
      title(ylab="mBAF",line=2)
      mtext(paste("Individual",patient_id,date_level_id),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }
    
  cp = current.position - snp_pos.start 

    
    segments(cp+annot$start_position, annot$mbaf, cp+annot$stop_position, annot$mbaf, col="red")
    ## profile = rbind(cbind(cp+annot$start_position, annot$mbaf),
    ##                 cbind(cp+annot$stop_position, annot$mbaf))
    ## profile = profile[order(profile[,1]),]
    ## lines(profile[,1], profile[,2], col="red")
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")

  ## Extra panel for annotations only.
  current.position = 0
  for(chri in 1:length(chrarr)){
  ##for(chri in 1:1){
    chr = chrarr[chri]
    ## Get snp_ids for 1-23 chromosomes
    chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
    
    snp_id.start = chr.start.p.q.stop[,"start"] 
    snp_id.stop = chr.start.p.q.stop[,"stop"]
    
    snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
      "where snp_id>=",snp_id.start,"and",
      "snp_id <=",snp_id.stop,
      "order by position"))
    snp_pos.start = min(snp_pos[,"position"])
    snp_pos.stop = max(snp_pos[,"position"])
    prev_ref_id = 0
    
    
      
    ## New annotations, based on compressed features
    annot2=dbGetQuery(con,paste("select * from compressed_patient_feature_states",
      "where sample_id=",sample_id," and reference_id=",reference_id,
      "and stop_snp_id >=",snp_id.start," and start_snp_id <=",snp_id.stop))

    if(chr==1){
      plot(0,0,ylim=c(ylim.lower.baf-y.offset.baf,ylim.upper.baf+y.offset.baf),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i",pch=NA)
      axis(2,at=c(0,0.25,0.5,0.75,1),
           labels=F,
           tick=T,las=2)
      axis(2,at=c(0.12,0.37,0.62,0.88),
           labels=c("HD","SD","CNLOH","GN"),
           tick=F,las=2)

      ##title(ylab="Annotations",line=2)
      mtext(paste("Individual",patient_id,date_level_id),side=3,line=0,at=total.base.pairs/2)
      current.position = snp_pos.start
    }

  state_0 = annot2[annot2[,"state_phased"]=="0",c("start_position","stop_position")]
  state_A = annot2[annot2[,"state_phased"]=="A",c("start_position","stop_position")]
  state_B = annot2[annot2[,"state_phased"]=="B",c("start_position","stop_position")]
  state_AA = annot2[annot2[,"state_phased"]=="AA",c("start_position","stop_position")]
  state_BB = annot2[annot2[,"state_phased"]=="BB",c("start_position","stop_position")]
  state_AAB = annot2[annot2[,"state_phased"]=="AAB",c("start_position","stop_position")]
  state_BBA = annot2[annot2[,"state_phased"]=="BBA",c("start_position","stop_position")]
  state_AAA = annot2[annot2[,"state_phased"]=="AAA",c("start_position","stop_position")]
  state_BBB = annot2[annot2[,"state_phased"]=="BBB",c("start_position","stop_position")]
  state_gain = annot2[annot2[,"state_phased"]=="AABB",c("start_position","stop_position")]

    
  cp = current.position - snp_pos.start 

  rect(cp+state_AA[,1],rep(0.51-y.offset.baf,nrow(state_AA)),
       cp+state_AA[,2],rep(0.74+y.offset.baf,nrow(state_AA)),
       col=sAA.col,border=sAA.col)

  rect(cp+state_BB[,1],rep(0.51-y.offset.baf,nrow(state_BB)),
       cp+state_BB[,2],rep(0.74+y.offset.baf,nrow(state_BB)),
       col=sAA.col,border=sAA.col)
    

  rect(cp+state_AAA[,1],rep(0.76-y.offset.baf,nrow(state_AAA)),
       cp+state_AAA[,2],rep(1+y.offset.baf,nrow(state_AAA)),
       col=sgain.col,border=sgain.col)
  
  rect(cp+state_BBB[,1],rep(0.76-y.offset.baf,nrow(state_BBB)),
       cp+state_BBB[,2],rep(1+y.offset.baf,nrow(state_BBB)),
       col=sgain.col,border=sgain.col)

  rect(cp+state_gain[,1],rep(0.76-y.offset.baf,nrow(state_gain)),
       cp+state_gain[,2],rep(1+y.offset.baf,nrow(state_gain)),
       col=sgain.col,border=sgain.col)

  rect(cp+state_AAB[,1],rep(0.76-y.offset.baf,nrow(state_AAB)),
       cp+state_AAB[,2],rep(1+y.offset.baf,nrow(state_AAB)),
       col=sgain.col,border=sgain.col)
  rect(cp+state_BBA[,1],rep(0.76-y.offset.baf,nrow(state_BBA)),
       cp+state_BBA[,2],rep(1+y.offset.baf,nrow(state_BBA)),
       col=sgain.col,border=sgain.col)

  rect(cp+state_A[,1],rep(0.26-y.offset.baf,nrow(state_A)),
       cp+state_A[,2],rep(0.49+y.offset.baf,nrow(state_A)),
       col=sA.col,border=sA.col)
  
  rect(cp+state_B[,1],rep(0.26-y.offset.baf,nrow(state_B)),
       cp+state_B[,2],rep(0.49+y.offset.baf,nrow(state_B)),
       col=sA.col,border=sA.col)

  rect(cp+state_0[,1],rep(0-y.offset.baf,nrow(state_0)),
       cp+state_0[,2],rep(0.24+y.offset.baf,nrow(state_0)),
       col=s0.col,border=s0.col)


    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    abline(v=current.position,col="white")
    abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    } else if(chr==24){
      mtext("Y", side = 1, line = 1, at = prev.position +  (current.position - prev.position)/2)
    } else {
      mtext(chr, side = 1, line = 1, at = prev.position + (current.position - prev.position)/2)
    }
  }
  axis(1,at=xaxsticks,labels=F,tick=T)
  box(col="black")

  
  dev.off()
  print(paste("Patient",patient_id,"sample",sample_id,"done."))
}



chrarr = c(1:24)
total.base.pairs = 0
chr.start = NULL
chr.stop = NULL
chr.snp_id.start = NULL
chr.snp_id.stop = NULL
for(chri in 1:length(chrarr)){
  chr = chrarr[chri]
  chr.start.p.q.stop = dbGetQuery(con,paste("select * from snps_annotation where chromosome=",chr)) 
  snp_id.start = chr.start.p.q.stop[,"start"] 
  snp_id.stop = chr.start.p.q.stop[,"stop"]
  chr.snp_id.start[chr] = snp_id.start
  chr.snp_id.stop[chr] = snp_id.stop  
  snp_pos = dbGetQuery(con,paste("select snp_id, position from snps",
    "where snp_id>=",snp_id.start,"and",
    "snp_id <=",snp_id.stop,
    "order by position"))
  snp_pos.start = min(snp_pos[,"position"])
  snp_pos.stop = max(snp_pos[,"position"])
  total.base.pairs = total.base.pairs + snp_pos.stop - snp_pos.start +1
  chr.start[chr] = snp_pos.start
  chr.stop[chr] = snp_pos.stop
}

CairoSVG("supplemental_figure_blood_samples_technical_replicates.svg",width=6.5,height=3,pointsize=8)
  par(mfrow=c(4,1),mar=c(2,3,1,1))  
  current.position = 0
  prev.position = 0
  xaxsticks = NULL

  for(chr in 1:24){
    snp_pos.start = chr.start[chr]
    snp_pos.stop = chr.stop[chr]
    if(chr%%2==0){      
      red = "gray"
      green = "gray"
    } else {
      red = "black"
      green = "black"
    }
    if(chr==1){
      plot(0,0,
           ylim=c(0,n),
           xlim=c(0,total.base.pairs),
           main=NA, xlab=NA, ylab=NA,
           xaxt="n",xaxs="i",yaxt="n",yaxs="i",pch=NA)
      ##abline(h=c(1:n),col="black") 
      axis(2,at=c(1:n)-0.5, labels=c(13:1),tick=F,las=1,line=0)
      axis(2,at=c(0:n), labels=NA,tick=T,las=1,line=0)
      ## axis(2,at=c(0.5:38.5), labels=rep(c(13:1),3), tick=F,las=2,line=0)
      ## axis(2,at=6.5, labels="2nd sampling interval", tick=F, las=0,line=1)
      ## axis(2,at=19.5, labels="1st sampling interval", tick=F, las=0,line=1)
      ## axis(2,at=32.5, labels="Baseline endo", tick=F, las=0,line=1)
      ## axis(2,at=19.5, labels="SGA", tick=F, las=0,line=2)
      current.position = snp_pos.start
    }
    
    cp = current.position - snp_pos.start 
    ##abline(v=cp,col="black") 
    if(chr%%2==0){
      rect(cp,0,cp+snp_pos.stop-snp_pos.start+1,n,col="grey90",border=NA)
    }
    ## Get annotation for current chromosome
    annot1 = all.sga[all.sga[,2]==chr,]
    ind = c(13:1)
    for(i in 1:n){
        if(length(annot1[annot1[,1]==ind[i] & annot1[,5]>0.0,1])>0){
          block = annot1[annot1[,1]==ind[i] & annot1[,5]>0.0,]
          rect(cp+block[,3],rep(i-1,nrow(block)),
               cp+block[,4],i-1+block[,5],
               col="black",border="black")
        }
    }
    
    ## Advance the current location
    prev.position = current.position
    current.position = current.position + snp_pos.stop - snp_pos.start +1
    ##abline(v=current.position,col="white")
    ##abline(v=prev.position,col="white")
    ## Plot chromosome character
    if(chr==23){
      mtext("X", side = 1, line = 0, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    } else if(chr==24){
      mtext("Y", side = 1, line = 0, at = prev.position +  (current.position - prev.position)/2)
    } else if(chr %in% c(1:22)){
      mtext(chr, side = 1, line = 0, at = prev.position + (current.position - prev.position)/2)
      xaxsticks = c(xaxsticks,current.position)
    }
  }
  mtext(side=1,line=1,at=total.base.pairs/2,text="Chromosome")
  mtext(side=2,line=2,at=n/2,text="Participant ID")
  box(col="black")
  dev.off()






mtext(side=1,line=1,"R of blood #1")
mtext(side=2,line=1,"R of blood #2")
abline(lm(y.var ~ x.var),col="red")




plot(x.var, y.var, title="Individual f",xlab=NA,ylab=NA)
mtext(side=1,line=1,"R of blood #1")
mtext(side=2,line=1,"R of blood #2")
abline(lm(y.var ~ x.var),col="red")

plot(x.var, y.var, title="Individual f vs. m",xlab=NA,ylab=NA)
mtext(side=1,line=1,"R of blood #1 (Individual f)")
mtext(side=2,line=1,"R of blood #1 (Individual m)")
abline(lm(y.var ~ x.var),col="red")

plot(x.var, y.var, title="Individual f vs. m",xlab=NA,ylab=NA)
mtext(side=1,line=1,"R of blood #2 (Individual f)")
mtext(side=2,line=1,"R of blood #2 (Individual m)")
abline(lm(y.var ~ x.var),col="red")



for(p in 1:length(patients)){
##  print(paste("Processing patient",p))
  patient_id = patients[p]

  ## Adjust y-axis
  if(pnames[p]%in%nsga500){
    sga_events_max=500
    sgayaxs = c("0","250", "500")
  } else if(pnames[p]%in%nsga1000){
    sga_events_max=1000
    sgayaxs = c("0","500", "1000")
  } else if(pnames[p]%in%nsga2200){
    sga_events_max=2200
    sgayaxs = c("0","1000", "2000")
  }


  features = dbGetQuery(con, paste("select * from compressed_patient_feature_states",
 	   		   "where patient_id=",patient_id," order by sample_id"))
  ## Convert features to 0/1
  features[features[,"state_phased"]!="AB","state_phased"] = 1
  features[features[,"state_phased"]=="AB","state_phased"] = 0
  samples = all.samples[all.samples[,"patient_id"]==patient_id,]

  comparisons = expand.grid(samples[samples[,"patient_id"]==patients[p],"sample_id"],
    samples[samples[,"patient_id"]==patients[p],"sample_id"])
  ## Exclude self comparisons
  comparisons =  comparisons[comparisons[,1]!=comparisons[,2],]
  distogram = as.data.frame(matrix(NA,nrow=nrow(comparisons),ncol=8))
  colnames(distogram) = c("sample_1","sample_2","spatial_distance","sga_distance","sga_events_distance","time_distance","date_1","date_2")
  distogram[,1] = comparisons[,1]
  distogram[,2] = comparisons[,2]
  distogram[,7] = samples[match(distogram[,"sample_1"],samples[,"sample_id"]),"endoscopy_date"]
  distogram[,8] = samples[match(distogram[,"sample_2"],samples[,"sample_id"]),"endoscopy_date"]
  

  for(i in 1:nrow(comparisons)){
    s1 = features[features[,"sample_id"]==comparisons[i,1],"state_phased"]
    s2 = features[features[,"sample_id"]==comparisons[i,2],"state_phased"]
    feat_sizes = features[features[,"sample_id"]==comparisons[i,1],"stop_position"] - features[features[,"sample_id"]==comparisons[i,1],"start_position"]+1
    
    ## Calculate amount (Mb) of the genome that is different, and the number of SGA events that are different
    distogram[i,"sga_distance"] = sum(feat_sizes[s1!=s2])/1000000
    distogram[i,"sga_events_distance"] = sum(s1!=s2)
    distogram[i,"spatial_distance"] = abs(
               abs(samples[samples[,"sample_id"]==comparisons[i,1],"endoscopy_les"] - samples[samples[,"sample_id"]==comparisons[i,1],"biopsy_level"]) -
               abs(samples[samples[,"sample_id"]==comparisons[i,2],"endoscopy_les"] - samples[samples[,"sample_id"]==comparisons[i,2],"biopsy_level"]))

    ## distogram[i,"spatial_distance"] = abs(samples[samples[,"sample_id"]==comparisons[i,1],"biopsy_level"] -
    ##                                samples[samples[,"sample_id"]==comparisons[i,2],"biopsy_level"])
    distogram[i,"time_distance"] = abs(round(as.numeric(difftime(
               samples[samples[,"sample_id"]==comparisons[i,1],"endoscopy_date"],
               samples[samples[,"sample_id"]==comparisons[i,2],"endoscopy_date"],units="days"))/365,digits=1))
  }

  timepoints = sort(unique(samples[,"endoscopy_date"]))
  plot(NA,NA,xlim=c(0,time_max),ylim=c(0,sga_events_max),xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  mtext(side=1,text="Time (years)",line=2)
  axis(1,at=c(0,5,10,15,20),lab=c(0,5,10,15,20),tick=T)  
  mtext(side=2,text="Number of SGA events distance",line=2)
  axis(2,at=sgayaxs,lab=sgayaxs,tick=T)
##  axis(2,at=c(0,500,1000,1500,2000),lab=c(0,500,1000,1500,2000),tick=T)
  mtext(side=3,text=paste("Individual",pnames[p]),line=1)
  x = abs(round(as.numeric(difftime(distogram[distogram[,"time_distance"]==0.0,"date_1"],timepoints[1],units="days"))/365,digits=1))
  y = distogram[distogram[,"time_distance"]==0.0,"sga_events_distance"]
  points(x,y,col="black")
  my.lm = lm(y ~ x)
  m = my.lm[[1]][2]
  x1 = min(x)
  x2 = max(x)
  y1 = m*(x1)+my.lm[[1]][1] 
  y2 = m*(x2)+my.lm[[1]][1]
    segments(x1,y1,x2,y2)
    legend("topright",border=NA,bty="n",
           legend=c(paste("R^2=",round(summary(my.lm)$r.squared,digits=2),sep=""),
             paste("p=",round(summary(my.lm)$coefficients[2,4],digits=2),sep="")))

##  abline(my.lm)

  ## plot(NA,NA,xlim=c(0,time_max),ylim=c(0,sga_events_max),xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  ## mtext(side=1,text="Time distance (years)",line=2)
  ## axis(1,at=c(0,5,10,15,20),lab=c(0,5,10,15,20),tick=T)  
  ## mtext(side=2,text="Number of SGA events distance",line=2)
  ## axis(2,at=c(0,500,1000,1500,2000),lab=c(0,500,1000,1500,2000),tick=T)
  ## mtext(side=3,text=paste("Individual",pnames[p]),line=1)
  ## points(distogram[,"time_distance"],distogram[,"sga_events_distance"],col="black")
  ## my.lm = lm(distogram[,"sga_events_distance"] ~ distogram[,"time_distance"])
  ## abline(my.lm)

  ## Time distances between biopsies separated by 1cm or less
  plot(NA,NA,xlim=c(0,time_max),ylim=c(0,sga_events_max),xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  mtext(side=1,text="Temporal distance (years)",line=2)
  axis(1,at=c(0,5,10,15,20),lab=c(0,5,10,15,20),tick=T)  
  mtext(side=2,text="Number of SGA events distance",line=2)
  axis(2,at=sgayaxs,lab=sgayaxs,tick=T)
##  axis(2,at=c(0,500,1000,1500,2000),lab=c(0,500,1000,1500,2000),tick=T)
  mtext(side=3,text=paste("Individual",pnames[p]),line=1)
  x = distogram[distogram[,"spatial_distance"]<=1,"time_distance"]
  y = distogram[distogram[,"spatial_distance"]<=1,"sga_events_distance"]
  points(x,y,col="black")
  my.lm = lm(y ~ x)
  m = my.lm[[1]][2]
  x1 = min(x)
  x2 = max(x)
  y1 = m*(x1)+my.lm[[1]][1] 
  y2 = m*(x2)+my.lm[[1]][1]
    segments(x1,y1,x2,y2)
    legend("topright",border=NA,bty="n",
           legend=c(paste("R^2=",round(summary(my.lm)$r.squared,digits=2),sep=""),
             paste("p=",round(summary(my.lm)$coefficients[2,4],digits=2),sep="")))

  ##abline(my.lm)

  
  ## plot(NA,NA,xlim=c(0,space_max),ylim=c(0,sga_events_max),xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  ## mtext(side=1,text="Space distance (cm)",line=2)
  ## axis(1,at=c(0,5,10,15),lab=c(0,5,10,15),tick=T)
  ## mtext(side=2,text="Number of SGA events distance",line=2)
  ## axis(2,at=c(0,500,1000,1500,2000),lab=c(0,500,1000,1500,2000),tick=T)
  ## mtext(side=3,text=paste("Individual",pnames[p]),line=1)
  ## points(distogram[,"spatial_distance"],distogram[,"sga_events_distance"],col="black")
  ## my.lm = lm(distogram[,"sga_events_distance"] ~ distogram[,"spatial_distance"])
  ## abline(my.lm)

  ## Only within the same time point
  plot(NA,NA,xlim=c(0,space_max),ylim=c(0,sga_events_max),xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  mtext(side=1,text="Spatial distance (cm)",line=2)
  axis(1,at=c(0,5,10,15),lab=c(0,5,10,15),tick=T)
  mtext(side=2,text="Number of SGA events distance",line=2)
  axis(2,at=sgayaxs,lab=sgayaxs,tick=T)
##  axis(2,at=c(0,500,1000,1500,2000),lab=c(0,500,1000,1500,2000),tick=T)
  mtext(side=3,text=paste("Individual",pnames[p]),line=1)
  x = distogram[distogram[,"time_distance"]==0.0,"spatial_distance"]
  y2 = distogram[distogram[,"time_distance"]==0.0,"sga_events_distance"]
  points(x,y2,col="black")
  my.lm = lm(y2 ~ x)
  if(is.na(my.lm$coefficients[2])){
    ##abline(v=x)
  } else {
    m = my.lm[[1]][2]
    x1 = min(x)
    x2 = max(x)
    y1 = m*(x1)+my.lm[[1]][1] 
    y2 = m*(x2)+my.lm[[1]][1]
    segments(x1,y1,x2,y2)
    legend("topright",border=NA,bty="n",
           legend=c(paste("R^2=",round(summary(my.lm)$r.squared,digits=2),sep=""),
             paste("p=",round(summary(my.lm)$coefficients[2,4],digits=2),sep="")))

##    abline(my.lm)
  }

  print(paste("Individual",pnames[p],wilcox.test(y,y2)$p.value))
  
}
dev.off()



