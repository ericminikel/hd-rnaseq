# Eric Minikel
# ki-coverage.r
# Perform characterization of HTT transcript in HD knock-in allelic series data

require(sqldf)
require(reshape)
options(stringsAsFactors=FALSE)
setwd('c:/sci/034ea/analysis/1/')

# matrix of base pair coverage for every base in Htt
coverage = read.table('allsamples.cov',header=TRUE)
names(coverage)

# list of filesizes of BAMs, as a proxy for library size
# if you want true library size in reads you could go back and re-calculate with samtools view -c
filesize = read.table('bam.filesize',header=FALSE)

# column indices of coverage: 1=contig, 2=location, 195=test, 3-194=data
covmat = as.matrix(coverage[,3:194]) # matrix of just the real data
samplenames = colnames(covmat)

# allow Pythonic negative indexing
negdex = function(v,idx) {
  return (v[length(v) + 1 + idx]  )
}

# extract metadata from sample names
lane = sapply(strsplit(samplenames,"_"),negdex,idx=-1)
sex  = sapply(strsplit(samplenames,"_"),negdex,idx=-2)
zyg  = sapply(strsplit(samplenames,"_"),negdex,idx=-3)
qlen = sapply(strsplit(samplenames,"_"),negdex,idx=-4)
tiss = sapply(strsplit(samplenames,"_"),negdex,idx=-5)
id   = sapply(strsplit(samplenames,"_"),"[[",1)
sid  = samplenames
metadata=data.frame(id,lane,sex,zyg,qlen,tiss,sid)

# parse sample number such as "845" out of the full sample name
metadata$sampleno = sapply(sapply(sapply(metadata$sid,strsplit,"_"),"[[",1),substr,2,8)
# do same for the full path of BAM files
filesize$sampleno = sapply(sapply(sapply(filesize$V2,strsplit,"[/|_]"),"[[",6),substr,1,8)
# grab the BAM library size in bytes and add to metadata
metadata$bambytes = filesize$V1[match(metadata$sampleno,filesize$sampleno)]

# convert Q length to number
metadata$q = as.numeric(sapply(metadata$qlen,substr,start=2,stop=4))

# add colors for tissues for plotting
coltiss = data.frame(tiss=c("Liver","striatum","cortex"),k=c("brown","orange","blue"))
metadata$tk = coltiss$k[match(metadata$tiss,coltiss$tiss)]

# add colors for q lengths for plotting
qs = sort(unique(metadata$q))
qs_inverse = 16**2 - qs # have high Q be black, low Q be gray
qh = as.hexmode(qs_inverse)
qk = paste('#',qh,qh,qh,sep='')
colq = data.frame(q=qs,k=qk)
metadata$qk = colq$k[match(metadata$q,colq$q)]

# check that the metadata table turned out how I wanted it
head(metadata)

# write out the metadata table for other uses
write.table(metadata,'metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE)

# data frame of the genomic sites to which the rows in the coverage matrix (covmat) refer
sites = data.frame(coverage[,c(1,2)])
Mmu5 = sites$contig == '5'
ki = sites$contig != '5'

# now I have three things:
# 1. covmat contains the data
# 2. sites contains the genomic location information for the rows of covmat
# 3. metadata contains the sample information for the columns of covmat


# plot total depth over all mouse Htt for each tissue
for (tissue in unique(metadata$tiss)) {
  png(paste(tissue,'.totalcov.png',sep=''),width=600,height=400)
  depth = rowSums(covmat[,metadata$tiss==tissue]) # sum depth
  plot(sites$base[sites$contig=='5'],depth[sites$contig=='5'],type='h',lwd=2,
       main=paste(tissue,' total depth',sep=''),xlab='Mmu5 location',
       ylab='base pair depth')
  dev.off()
}


# need to get exon structure
# Ensembl tracking_id for proper Htt transcript is ENSMUST00000080036
# grab exon structure from the gtf file using these bash commands:
# cd /data/HD/analysis/evm/038ea/analysis/1/cufflinks/1/847_Liver_Q92_HET_F_L5
# cat isoforms.fpkm_tracking | grep Htt | less
# cd /data/HD/public/GRCm38.73/
# cat *.gtf | grep ENSMUST00000080036 | grep -w exon | wc -l # 67 - correct
# cat *.gtf | grep ENSMUST00000080036 | grep -w exon | cut -f4,5 > Mmu5.Htt.exons

exons = read.table('Mmu5.Htt.exons',header=FALSE,sep='\t') 
colnames(exons) = c('start','stop')
sites$exon = findInterval(sites$base,exons$start) # extract exon numbers corresponding to each base in Htt
sites$exon[sites$contig != '5'] = NA # assign no exon number for the knock-in construct
plot(sites$base[sites$contig=='5'],sites$exon[sites$contig=='5'],type='h') # verify that exon numbers were assigned correctly

# this function is like findInterval but it makes sure the basepair is both
# greater than a start index AND lower than a stop index
# like BEDtools, this function uses 0-based starts and 1-based stops
findBoundedInterval = function(basepair,maptable) {
  temp1 = (which(maptable$start <= basepair & maptable$stop > basepair))
  if (length(temp1) == 0) {
    return (0)
  } else {
    return (temp1)
  }
}

# use this function to calculate the exon number for every base
sites$exon = sapply(sites$base,findBoundedInterval,maptable=exons)

# infer intron boundaries from the exons table
introns = data.frame(exons$stop[1:66],exons$start[2:67])
colnames(introns) = c('start','stop')
# calculate intron number for every base
sites$intron = sapply(sites$base,findBoundedInterval,maptable=introns)

# confirm that exons and introns were assigned correctly
plot(sites$base[Mmu5],sites$exon[Mmu5],type='h',col='blue')
points(sites$base[Mmu5],sites$intron[Mmu5],type='h',col='gray')
points(sites$base[Mmu5],sites$exon[Mmu5],type='h',col='blue')

# calculate total depth across all 192 samples at each base
sites$sumdepth = rowSums(covmat)
# plot total depth across endogenous Mmu5 Htt transcript colored by intron/exon status
exonic = sites$exon > 0
plot(sites$exon[Mmu5 & exonic], sites$sumdepth[Mmu5 & exonic], pch='.')
plot(sites$intron[Mmu5 & !exonic], sites$sumdepth[Mmu5 & !exonic], pch='.')

# calculate average depth of each exon
sql_query = "
select   avg(s.sumdepth) avdepth
from     sites s
where    s.exon > 0
group by s.exon
order by s.exon
;
"
# note that in the sql query I average across all the _base pairs_ in the exon
# while in the line below I average across all the 192 _samples_ in the dataset
exons$avdepth = unlist(sqldf(sql_query))/length(samplenames)
# plot average depth by exon number
png('av.depth.by.exon.png',width=600,height=400)
barplot(exons$avdepth,col='black',main='Average depth by exon',xlab='Exon',
        names.arg=1:dim(exons)[1],ylab='Average total depth',
        cex.names=.6,las=3)
dev.off()

# intron average depth, exactly like above
sql_query = "
select   avg(s.sumdepth) avdepth
from     sites s
where    s.intron > 0
group by s.intron
order by s.intron
;
"
introns$avdepth = unlist(sqldf(sql_query))/length(samplenames)
png('av.depth.by.intron.png',width=600,height=400)
barplot(introns$avdepth,col='black',main='Average depth by intron',xlab='Intron',
        names.arg=1:dim(introns)[1],ylab='Average total depth',
        cex.names=.6,las=3)
dev.off()

# intron exon ratio
# potentially a measure of intron inclusion normalized for 3' enrichment
# ratio of intron av depth to that of the exon to its right
# right as opposed to left or both because this avoids exon 1
# which has the complicated knock in situation
ieratio = introns$avdepth / exons$avdepth[2:dim(exons)[1]]
png('ieratio.png',width=600,height=400)
barplot(ieratio,col='black',main='Intron/exon ratio',
        xlab='Intron (n) / exon (n+1)',ylab='Depth ratio',
        names.arg=1:dim(introns)[1],cex.names=.6,las=3)
dev.off()

# make one general plot of total Htt coverage colored by exon vs. intron for ALL samples
depth = rowSums(covmat[sites$contig=='5',])
png('total.htt.coverage.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron > 0], depth[sites$intron > 0], col='gray', type='h',
     main='Total coverage of Htt, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$contig=='5']))
points(sites$base[sites$exon > 0], depth[sites$exon > 0], col='blue', type='h')
legend('topleft',c('intron','exon'),col=c('gray','blue'),pch=NA,lwd=3)
dev.off()

# highlight the few areas my presentation will focus on - around introns 1, 49, 59 and 3'UTR
png('total.htt.coverage.areas.of.interest.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron > 0], depth[sites$intron > 0], col='gray', type='h',
     main='Total coverage of Htt, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$contig=='5']))
points(sites$base[sites$exon > 0], depth[sites$exon > 0], col='lightblue', type='h')
points(sites$base[sites$exon %in% c(1,2,49,50,59,60,67)], depth[sites$exon %in% c(1,2,49,50,59,60,67)], col='blue', type='h')
points(sites$base[sites$intron %in% c(1,49,59)], depth[sites$intron %in% c(1,49,59)], col='#777777', type='h')
dev.off()

# focus on intron 1
png('total.intron1.coverage.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==1], depth[sites$intron ==1], col='gray', type='h',
     main='Total coverage of intron 1, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth[sites$intron ==1]), xlim=range(sites$base[sites$intron ==1 | sites$exon %in% c(1,2)]))
points(sites$base[sites$exon %in% c(1,2)], depth[sites$exon %in% c(1,2)], col='blue', type='h')
legend('topleft',c('intron','exon'),col=c('gray','blue'),pch=NA,lwd=3)
dev.off()

# focus on intron 49
png('total.intron49.coverage.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==49], depth[sites$intron ==49], col='gray', type='h',
     main='Total coverage of intron 49, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth[sites$intron ==49]), xlim=range(sites$base[sites$intron ==49 | sites$exon %in% c(49,50)]))
points(sites$base[sites$exon %in% c(49,50)], depth[sites$exon %in% c(49,50)], col='blue', type='h')
legend('topleft',c('intron','exon'),col=c('gray','blue'),pch=NA,lwd=3)
dev.off()

# focus on intron 59
png('total.intron59.coverage.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==59], depth[sites$intron ==59], col='gray', type='h',
     main='Total coverage of intron 59, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth[sites$intron ==59]), xlim=range(sites$base[sites$intron ==59 | sites$exon %in% c(59,60)]))
points(sites$base[sites$exon %in% c(59,60)], depth[sites$exon %in% c(59,60)], col='blue', type='h')
legend('topleft',c('intron','exon'),col=c('gray','blue'),pch=NA,lwd=3)
dev.off()

# focus on 3'UTR
png('total.3utr.coverage.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron > 66], depth[sites$intron > 66], col='gray', type='h',
     main="Total coverage of Htt 3'UTR (exon 67), all samples", xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$intron > 66 | sites$exon > 66]))
points(sites$base[sites$exon > 66], depth[sites$exon > 66], col='blue', type='h')
legend('topleft',c('intron','exon'),col=c('gray','blue'),pch=NA,lwd=3)
dev.off()


# does relationship between exon and depth depend on:
# 1. Q length
# 2. Zygosity
# 3. Sex
# 4. Tissue

# create a table of exon data by sample to run models on
exon_avgs = aggregate(covmat, list(sites$exon), mean)# get average depth for each exon for each sample
exon_avgs_rel = melt(exon_avgs,id=c("Group.1")) # convert to a relational database 3-column format
colnames(exon_avgs_rel) = c("exon","sid","avdepth")
exon_avgs_rel = exon_avgs_rel[exon_avgs_rel$exon > 0,] # remove data for all introns combined
exon_avgs_data = merge(exon_avgs_rel, metadata, by="sid") # merge in the metadata for each sample

# see which variables modulate the relationship between exon number and depth
# complicated to know how to model zygosity because
# WT are available only for striatum
m = lm(avdepth ~ exon + sex + zyg * q * tiss, data=exon_avgs_data)
summary(m)

# note that lm is essentially ancova, where the significance
# of interaction term * tells you if the slope is different
# based on the other variable

m = lm(avdepth ~ exon * tiss, data=exon_avgs_data)
summary(m)
# yes: the higher the exon, the more coverage (.93)
# Liver much less than cortex (-18) and striatum slightly more (1.12)
# exon:tissLiver - .74, i.e. the rise is less rapid in Liver
# exon:tissstriatum, not significantly diff in striatum vs. cortex

m = lm(avdepth ~ exon * tiss * q, data=subset(exon_avgs_data, zyg!='WT'))
summary(m)

# plot exon vs. depth for diff tissues
png('exon.vs.depth.by.tissue.png',width=600,height=400)
plot(exon_avgs_data$exon, exon_avgs_data$avdepth, pch=19, col=exon_avgs_data$tk,
     main = "Exon number vs. depth by tissue", xlab='Htt exon number',
     ylab = "Average depth")
legend('topleft',coltiss$tiss,col=coltiss$k,pch=19)
dev.off()

# plot exon vs. depth for diff q lengths
png('exon.vs.depth.by.q.png',width=600,height=400)
het = exon_avgs_data$zyg == 'HET'
stri = exon_avgs_data$tiss == 'striatum'
# plot exon vs. depth for diff Q lengths, within striatum
plot(exon_avgs_data$exon[het&stri], exon_avgs_data$avdepth[het&stri], pch=19, 
     col=exon_avgs_data$qk[het&stri],
     main = "Exon number vs. depth by Q length", xlab='Htt exon number',
     ylab = "Average depth", cex=.8)
legend('topleft',legend=colq$q,col=colq$k,pch=19)
for (q in qs) {
  abline(lm(exon_avgs_data$avdepth[het&stri&exon_avgs_data$q==q] ~ exon_avgs_data$exon[het&stri&exon_avgs_data$q==q]),col=colq$k[colq$q==q])
}
dev.off()

# print out the slopes for each q length in the model
for (q in qs) {
  m = lm(exon_avgs_data$avdepth[het&stri&exon_avgs_data$q==q] ~ exon_avgs_data$exon[het&stri&exon_avgs_data$q==q])
  print(q)
  print(summary(m)$coefficients[2,1])
}
# results:
# [1] 20
# [1] 1.03311
# [1] 80
# [1] 1.053802
# [1] 92
# [1] 0.8473141
# [1] 111
# [1] 1.025569
# [1] 140
# [1] 0.8107097
# [1] 175
# [1] 0.8132042
# so the slopes are ~1.00 for Q20, Q80 and Q111, and ~.80 for Q111, Q140z and Q175z

m = lm(avdepth ~ exon * qlen, data = exon_avgs_data[het&stri,])
summary(m)
# the intercepts for diff q lengths are not significantly different from Q111
# except for Q92 and Q140z which are borderline at p = .04 and p = .06 respectively
# on the other hand, the slopes for Q92, Q140z and Q175z are all highly significantly
# different, p < 1e-8
# tl;dr: same coverage but less 5' degradation in those three lines.

# create same sort of table for introns by sample
intron_avgs = aggregate(covmat, list(sites$intron), mean)
intron_avgs_rel = melt(intron_avgs,id=c("Group.1"))
colnames(intron_avgs_rel) = c("intron","sid","avdepth")
intron_avgs_rel = intron_avgs_rel[intron_avgs_rel$intron > 0,]
intron_avgs_data = merge(intron_avgs_rel, metadata, by="sid")


# plot intron vs. depth for diff tissues
png('intron.vs.depth.by.tissue.png',width=600,height=400)
plot(intron_avgs_data$intron, intron_avgs_data$avdepth, pch=19, col=intron_avgs_data$tk,
     main = "intron number vs. depth by tissue", xlab='Htt intron number',
     ylab = "Average depth")
legend('topleft',coltiss$tiss,col=coltiss$k,pch=19)
dev.off()

# plot intron vs. depth for diff q lengths
png('intron.vs.depth.by.q.png',width=600,height=400)
het = intron_avgs_data$zyg == 'HET'
stri = intron_avgs_data$tiss == 'striatum'
# plot intron vs. depth for diff Q lengths, within striatum
plot(intron_avgs_data$intron[het&stri], intron_avgs_data$avdepth[het&stri], pch=19, 
     col=intron_avgs_data$qk[het&stri],
     main = "intron number vs. depth by Q length", xlab='Htt intron number',
     ylab = "Average depth", cex=.8)
legend('topleft',legend=colq$q,col=colq$k,pch=19)
for (q in qs) {
  abline(lm(intron_avgs_data$avdepth[het&stri&intron_avgs_data$q==q] ~ intron_avgs_data$intron[het&stri&intron_avgs_data$q==q]),col=colq$k[colq$q==q])
}
dev.off()

# double check format of data
dim(covmat)
head(exon_avgs_data)


# re-plot detail for intron 1 and 49
png('intron.49.detail.png',width=600,height=400)
plot(sites$base[sites$intron==49], sites$sumdepth[sites$intron==49], type='h', col='gray',
     xlab = 'Mmu5 base pair position', ylab = 'Total depth', main = 'Coverage of intron 49')
points(sites$base[sites$exon==49], sites$sumdepth[sites$exon==49], type='h', col='blue')
points(sites$base[sites$exon==50], sites$sumdepth[sites$exon==50], type='h', col='blue')
dev.off()

png('intron.1.detail.png',width=600,height=400)
plot(sites$base[sites$intron==1], sites$sumdepth[sites$intron==1], type='h', col='gray',
     xlab = 'Mmu5 base pair position', ylab = 'Total depth', main = 'Coverage of intron 1')
points(sites$base[sites$exon==1], sites$sumdepth[sites$exon==1], type='h', col='blue')
points(sites$base[sites$exon==2], sites$sumdepth[sites$exon==2], type='h', col='blue')
dev.off()

# examine: do the peaks in intron 1 coverage match predicted polyA sites??
introns[1,]
# chr5:34762109-34782742
# Retrieved from UCSC genome browser
# fed into polyah, accessed 2013-11-22
# http://linux1.softberry.com/berry.phtml?topic=polyah&group=programs&subgroup=promoter
# results:
# >mm10_dna range=chr5:34762109-34782742 5'pad=0 3'pad=0 stran
# Length of sequence-     20634
# 13 potential polyA sites were predicted
# Pos.:    678 LDF-  1.44
# Pos.:   1234 LDF-  4.09
# Pos.:   1285 LDF-  2.66
# Pos.:   2863 LDF-  5.24
# Pos.:   3535 LDF-  0.44
# Pos.:   4411 LDF-  1.57
# Pos.:   4659 LDF-  3.35
# Pos.:   8234 LDF-  6.39
# Pos.:  13574 LDF-  4.00
# Pos.:  13579 LDF-  3.79
# Pos.:  13583 LDF-  3.45
# Pos.:  18569 LDF-  4.31
# Pos.:  19641 LDF-  5.44
polyah = read.table('mmu5.htt.intron1.polyah.txt',skip=3) # grab the above data into a data frame
polyah = polyah[,c(2,4)] # only the relative base pair position and prediction strength columns
colnames(polyah) = c("relpos","ldf")
polyah$abspos = polyah$relpos + introns$start[1] # convert to absolute base pair position

# plot intron 1 with polyAH prediction sites superimposed
png('intron.1.detail.polyah.png',width=600,height=400)
plot(sites$base[sites$intron==1], sites$sumdepth[sites$intron==1], type='h', col='gray',
     xlab = 'Mmu5 base pair position', ylab = 'Total depth', main = 'Coverage of intron 1')
points(sites$base[sites$exon==1], sites$sumdepth[sites$exon==1], type='h', col='blue')
points(sites$base[sites$exon==2], sites$sumdepth[sites$exon==2], type='h', col='blue')
points(polyah$abspos, polyah$ldf*(1500/max(polyah$ldf)), type='h', col='red')
legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
dev.off()


# this giant loop is going to create one plot for every region of interest X every tissue X every genotype
ranges = list()
ranges[['KI_CONSTRUCT']] = sites$contig != '5'
ranges[['EXON1_INTRON1_EXON2']] = sites$exon %in% c(1,2) | sites$intron == 1
ranges[['EXON49_INTRON49_EXON50']] = sites$exon %in% c(49,50) | sites$intron == 49
ranges[['EXON59_INTRON59_EXON60']] = sites$exon %in% c(59,60) | sites$intron == 59
for (rangename in names(ranges)) {
  range = ranges[[rangename]]
  use_sites = sites[range,] # x-axis range of plot - where in Htt are we plotting
  ytop = max(rowMeans(covmat[range,])) # determine y axis maximum to be used for all plots of this region for consistency
  for (nowtiss in c('striatum','cortex','Liver')) { # loop over tissues
    genotypes = list()
    genotypes[['ALL_WT']] = metadata$zyg == 'WT' & metadata$tiss == nowtiss  # treat all WT mice together regardless of what Q length their littermates are
    genotypes[['Q20_HET']]  = metadata$q==20  & metadata$zyg == 'HET' & metadata$tiss == nowtiss
    genotypes[['Q80_HET']]  = metadata$q==80  & metadata$zyg == 'HET' & metadata$tiss == nowtiss
    genotypes[['Q92_HET']]  = metadata$q==92  & metadata$zyg == 'HET' & metadata$tiss == nowtiss
    genotypes[['Q111_HET']] = metadata$q==111 & metadata$zyg == 'HET' & metadata$tiss == nowtiss
    genotypes[['Q140_HET']] = metadata$q==140 & metadata$zyg == 'HET' & metadata$tiss == nowtiss
    genotypes[['Q175_HET']] = metadata$q==175 & metadata$zyg == 'HET' & metadata$tiss == nowtiss
    for (genotype in names(genotypes)) { # loop over genotypes
      png(paste(rangename,nowtiss,genotype,'cov.png',sep='.'),width=600,height=400)
      gt_depth = rowMeans(covmat[range, genotypes[[genotype]]]) # calculate avg depth for each base for this genotype
      plot(use_sites$base, gt_depth, type='h', col='gray', ylim=c(0,ytop), main=paste(rangename,genotype,nowtiss),
           xlab = 'Mmu5 base pair position', ylab='Average depth')
      points(use_sites$base[use_sites$exon > 0],  gt_depth[use_sites$exon > 0], type='h', col='blue')
      points(polyah$abspos, polyah$ldf*(ytop/max(polyah$ldf)), type='h', col='red') # superimpose polyAH precition sites again
      legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
      dev.off()
    }
  }
}


# now try to match up ki construct with mouse genome.
# mouse CAG tract in mm10 is Mmu5:34761957-34761978,
# CAG tract in ki construct is base 191-251
# human intron 1 in ki construct starts at base 400, exon 1 starts at base 51


sites$kiplotpos = NA # initialize a new column to hold X axis plotting position for combined plots
sites$kiplotpos[sites$contig != '5'] = sites$base[sites$contig != '5'] # for the ki construct, it is itself
sites$kiplotpos[sites$contig == '5' & sites$base %in% 34761957:34761978] = 191:(191+21) # map CAG tract to first 21 bases of ki CAG tract
sites$kiplotpos[sites$contig == '5' & sites$base %in% (34761957-191):34761957] = 0:191 # align pre-CAG pieces
sites$kiplotpos[sites$contig == '5' & sites$base %in% 34761978:(34761978+1700)] = 251:(251+1700) # align post-CAG pieces
# now the sections I've chosen to plot of the endogenous gene are 5:34761766-34763678
# and the ki construct corresponds to just part of that

# important note about above:
# if you actually smith-waterman align the ki construct and endogenous gene sequences, you'll find that the sequence similarity between
# human and mouse is not all in one contiguous block as I've assumed above. this ensuing analysis should be re-done more carefully lining
# up the bases between the ki construct and Htt gene. we discussed this at the 2013-12-12 meeting in Jim's office but it was decided to hold
# off on re-doing that until we get the actual sequences of the Wheeler and Zeitlin knock-in alleles from Tammy

# create boolean vectors for indexing
ki = sites$contig != '5'
mmu5 = sites$contig == '5'
cag = sites$kiplotpos %in% 191:251
exon1 = sites$kiplotpos %in% 51:400

# check depth of the two
par(mfrow=c(2,1))
plot(sites$kiplotpos[mmu5], sites$sumdepth[mmu5], type='h', col='violet')
points(sites$kiplotpos[ki], sites$sumdepth[ki], type='h', col='orange')
dev.off()

# now plot the ki construct and endogenous gene depth for each tissue-genotype combination
range = !is.na(sites$kiplotpos) #  indices to plot on x axis
ytop = 2*max(rowMeans(covmat[range,])) # maximum of ylim for plotting
for (nowtiss in c('striatum','cortex','Liver')) {
  genotypes = list()
  genotypes[['ALL_WT']] = metadata$zyg == 'WT' & metadata$tiss == nowtiss
  genotypes[['Q20_HET']]  = metadata$q==20  & metadata$zyg == 'HET' & metadata$tiss == nowtiss
  genotypes[['Q80_HET']]  = metadata$q==80  & metadata$zyg == 'HET' & metadata$tiss == nowtiss
  genotypes[['Q92_HET']]  = metadata$q==92  & metadata$zyg == 'HET' & metadata$tiss == nowtiss
  genotypes[['Q111_HET']] = metadata$q==111 & metadata$zyg == 'HET' & metadata$tiss == nowtiss
  genotypes[['Q140_HET']] = metadata$q==140 & metadata$zyg == 'HET' & metadata$tiss == nowtiss
  genotypes[['Q175_HET']] = metadata$q==175 & metadata$zyg == 'HET' & metadata$tiss == nowtiss
  
  for (genotype in names(genotypes)) {
    png(paste('ki-mmu5.overlay',nowtiss,genotype,'png',sep='.'),width=600,height=400)
    gt_depth = rowMeans(covmat[range, genotypes[[genotype]]])
    plot(NA,NA, xaxt='n', xlab = 'position within knock-in construct', ylab='average base pair depth',
         xlim=c(0,max(sites$kiplotpos,na.rm=TRUE)), ylim=c(0,ytop), main=paste('ki-mmu5 overlay',genotype,nowtiss))
    points(sites$kiplotpos[mmu5], rowMeans(covmat[mmu5,genotypes[[genotype]]]), type='l', lwd=4, col='gray')
    points(sites$kiplotpos[ki], rowMeans(covmat[ki,genotypes[[genotype]]]), type='l', lwd=4, col='orange')
    axis(side=1, at=c(0,51,191,251,400,951), labels = rep('',6), lwd.ticks=4, pos=0)
    mtext(side=1, at=c(0,(191+51)/2,(251+191)/2,(400+251)/2,(800+400)/2), 
          text=c("5'UTR","Exon 1","CAG","Exon 1","Intron 1") )
    legend('topright',c('Mmu5 coverage','KI construct coverage'),col=c('gray','orange'),lwd=4,pch=NA)
    dev.off()
  }
}




## POLYADENYLATION ANALYSIS
# this part is still a work in progress
# See polya-search.bash - I used that script (after trying many things that didn't work) to pull out the
# 5' alignment locations of reads that appeared to be polyadenylated. here I read that data into R and play
# with it:

# read polyA read 5' locations in
polya = read.table('polya.try4.txt',header=FALSE)$V1
length(polya)

# read a list of sites in the Htt gene where the reference genome contains >= 7 consecutive A
# these are likely spots where the polyA reads might be false positives, although my requirement that the
# polyA reads also have a soft-clipped CIGAR string should already have reduced the false positive rate
# considerably
genomic_7a_sites = read.table('mmu5.htt.7a.sites.txt',header=FALSE)
# define a region around the 7 A site to exclude from the polyA read analysis
# max read length is 75 so ignore 75 bp on either side
genomic_7a_sites$start = genomic_7a_sites$V1 - 75 
genomic_7a_sites$stop = genomic_7a_sites$V1 + 75
genomic_7a_sites = genomic_7a_sites[,c("start","stop")]
# now find all the indices in the polya table that might be false positives as they are near a genomic site with 7 As
# the inner sapply will return a list where some entries have > 1 index b/c they are near multiple sites, so wrap outer sapply to take first item from each
false_positive_matches = sapply(sapply(polya,findBoundedInterval,maptable=genomic_7a_sites),"[[",1)
polya = polya[false_positive_matches==0] # remove the potential false positives
length(polya)

# plot Htt coverage with introns of interest highlighted and polyA reads as red below the plot
png('total.htt.coverage.areas.of.interest.polya.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron > 0], depth[sites$intron > 0], col='gray', type='h',
     main='Total coverage of Htt, all samples, with polyA', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$contig=='5']))
points(sites$base[sites$exon > 0], depth[sites$exon > 0], col='lightblue', type='h')
points(sites$base[sites$exon %in% c(1,2,49,50,59,60,67)], depth[sites$exon %in% c(1,2,49,50,59,60,67)], col='blue', type='h')
points(sites$base[sites$intron %in% c(1,49,59)], depth[sites$intron %in% c(1,49,59)], col='#777777', type='h')
points(polya, rep(-1000, length(polya)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

# now go through and zoom in on each intron of interest
png('total.intron1.coverage.polya.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==1], depth[sites$intron ==1], col='gray', type='h',
     main='Total coverage of intron 1, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth[sites$intron ==1]), xlim=range(sites$base[sites$intron ==1 | sites$exon %in% c(1,2)]))
points(sites$base[sites$exon %in% c(1,2)], depth[sites$exon %in% c(1,2)], col='blue', type='h')
points(polya, rep(-1000, length(polya)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

png('total.intron49.coverage.polya.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==49], depth[sites$intron ==49], col='gray', type='h',
     main='Total coverage of intron 49, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth[sites$intron ==49]), xlim=range(sites$base[sites$intron ==49 | sites$exon %in% c(49,50)]))
points(sites$base[sites$exon %in% c(49,50)], depth[sites$exon %in% c(49,50)], col='blue', type='h')
points(polya, rep(-1000, length(polya)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

png('total.intron59.coverage.polya.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==59], depth[sites$intron ==59], col='gray', type='h',
     main='Total coverage of intron 59, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth[sites$intron ==59]), xlim=range(sites$base[sites$intron ==59 | sites$exon %in% c(59,60)]))
points(sites$base[sites$exon %in% c(59,60)], depth[sites$exon %in% c(59,60)], col='blue', type='h')
points(polya, rep(-1000, length(polya)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

png('total.intron61.coverage.polya.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron ==61], depth[sites$intron ==61], col='gray', type='h',
     main='Total coverage of intron 61, all samples', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=c(0,2000), xlim=range(sites$base[sites$intron ==61 | sites$exon %in% c(61,62)]))
points(sites$base[sites$exon %in% c(61,62)], depth[sites$exon %in% c(61,62)], col='blue', type='h')
points(polya, rep(-1000, length(polya)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

png('total.3utr.coverage.polya.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron > 66], depth[sites$intron > 66], col='gray', type='h',
     main="Total coverage of Htt 3'UTR (exon 67), all samples", xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$intron > 66 | sites$exon > 66]))
points(sites$base[sites$exon > 66], depth[sites$exon > 66], col='blue', type='h')
points(polya, rep(-1000, length(polya)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

# cluster polya sites to remove singletons - throw out ones that aren't within 100 of another.
polya = sort(polya)
keep = rep(TRUE,length(polya))
for (i in 2:(length(polya)-1)) {
  if (polya[i] > polya[i-1]+100 && polya[i] < polya[i+1]-100) {
    keep[i] = FALSE
  }
}
# manually inspect first and last
head(polya)
tail(polya)
polya_clustered = polya[keep] # choose the clustered reads
length(polya_clustered) # 308 raw, 276 after masking, 257 after clustering

# re-plot with only the clustered reads
png('total.htt.coverage.areas.of.interest.polya.clustered.all.samples.png',width=600,height=400)
plot(sites$base[sites$intron > 0], depth[sites$intron > 0], col='gray', type='h',
     main='Total coverage of Htt, all samples, with polyA clusters', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$contig=='5']))
points(sites$base[sites$exon > 0], depth[sites$exon > 0], col='lightblue', type='h')
points(sites$base[sites$exon %in% c(1,2,49,50,59,60,67)], depth[sites$exon %in% c(1,2,49,50,59,60,67)], col='blue', type='h')
points(sites$base[sites$intron %in% c(1,49,59)], depth[sites$intron %in% c(1,49,59)], col='#777777', type='h')
points(polya_clustered, rep(-1000, length(polya_clustered)), type='h', col='red')
legend('topleft',c('intron','exon','polyA read'),col=c('gray','blue','red'),pch=NA,lwd=3)
dev.off()

# for my own reference - plot with the exon and intron numbers like e1, i1, etc.
png('total.htt.coverage.numbering.guide.all.samples.png',width=1200,height=800)
plot(sites$base[sites$intron > 0], depth[sites$intron > 0], col='gray', type='h',
     main='Total coverage of Htt, all samples, with polyA clusters', xlab='Mmu5 base pair position',
     ylab='Total depth', ylim=range(depth), xlim=range(sites$base[sites$contig=='5']))
points(sites$base[sites$exon > 0], depth[sites$exon > 0], col='lightblue', type='h')
points(sites$base[sites$exon %in% c(1,2,49,50,59,60,67)], depth[sites$exon %in% c(1,2,49,50,59,60,67)], col='blue', type='h')
points(sites$base[sites$intron %in% c(1,49,59)], depth[sites$intron %in% c(1,49,59)], col='#777777', type='h')
points(polya_clustered, rep(-1000, length(polya_clustered)), type='h', col='red')
text(rowMeans(introns[,c("start","stop")]),rep(10000,dim(introns)[1]),labels=paste("i",rownames(introns),sep=''))
text(rowMeans(exons[,c("start","stop")]),rep(12000,dim(exons)[1]),labels=paste("e",rownames(exons),sep=''))
dev.off()

# now try to group the clusters so as to make one plot of each cluster
polya_group = polya_clustered %/% 10000
polya[duplicated(polya_group)]
polya_zoom_centers = sort(polya_clustered[!duplicated(polya_group)])
zoomctr=polya_zoom_centers[1] # sites on the x axis to zoom in on

for (zoomctr in polya_zoom_centers) {
  xlim = c(zoomctr-10000,zoomctr+10000) # plot +- 10kb around the zoom center
  exon = findBoundedInterval(zoomctr,maptable=exons)
  intron = findBoundedInterval(zoomctr,maptable=introns)
  zonetitle=''
  if (exon > 0) {
    zonetitle = paste("exon",exon,sep='')
  } else if (intron > 0) {
    zonetitle = paste("intron",intron,sep='')
  }
  range = sites$base > min(xlim) & sites$base < max(xlim)
  ylim = c(0, max(depth[range],na.rm=TRUE))
  polyacount = sum(polya_clustered %/% 10000 == zoomctr %/% 10000)
  png(paste('polya.zoom.',zoomctr,'.png',sep=''),width=1200,height=800)
  plot(sites$base[sites$intron > 0], depth[sites$intron > 0], col='gray', type='h',
       main=paste('Total coverage of Htt near ',polyacount,'polyA read cluster at',zonetitle), xlab='Mmu5 base pair position',
       ylab='Total depth', ylim=ylim, xlim=xlim)
  points(sites$base[sites$exon > 0], depth[sites$exon > 0], col='blue', type='h')
  points(polya_clustered[polya_clustered %/% 10000 == zoomctr %/% 10000], rep(-1000, length(polya_clustered[polya_clustered %/% 10000 == zoomctr %/% 10000])), type='h', col='red')
  text(rowMeans(introns[,c("start","stop")]),rep(max(ylim),dim(introns)[1]),labels=paste("i",rownames(introns),sep=''))
  text(rowMeans(exons[,c("start","stop")]),rep(max(ylim),dim(exons)[1]),labels=paste("e",rownames(exons),sep=''))
  dev.off()
}

