require(sqldf)
require(reshape)

options(stringsAsFactors=FALSE)
setwd('c:/sci/034ea/analysis/1/')
coverage = read.table('allsamples.cov',header=TRUE)
names(coverage)

# column indices of coverage: 1=contig, 2=location, 195=test, 3-194=data
covmat = as.matrix(coverage[,3:194])
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

# convert Q length to number
metadata$q = as.numeric(sapply(metadata$qlen,substr,start=2,stop=4))

# add colors for tissues
coltiss = data.frame(tiss=c("Liver","striatum","cortex"),k=c("brown","orange","blue"))
metadata$tk = coltiss$k[match(metadata$tiss,coltiss$tiss)]

# add colors for q lengths
qs = sort(unique(metadata$q))
qs_inverse = 16**2 - qs # have high Q be black, low Q be gray
qh = as.hexmode(qs_inverse)
qk = paste('#',qh,qh,qh,sep='')

colq = data.frame(q=qs,k=qk)
metadata$qk = colq$k[match(metadata$q,colq$q)]

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
# cd /data/HD/analysis/evm/038ea/analysis/1/cufflinks/1/847_Liver_Q92_HET_F_L5
# cat isoforms.fpkm_tracking | grep Htt | less
# cd /data/HD/public/GRCm38.73/
# cat *.gtf | grep ENSMUST00000080036 | grep -w exon | wc -l # 67 - correct
# cat *.gtf | grep ENSMUST00000080036 | grep -w exon | cut -f4,5 > Mmu5.Htt.exons

exons = read.table('Mmu5.Htt.exons',header=FALSE,sep='\t')
colnames(exons) = c('start','stop')
sites$exon = findInterval(sites$base,exons$start)
sites$exon[sites$contig != '5'] = NA
plot(sites$base[sites$contig=='5'],sites$exon[sites$contig=='5'],type='h')

# like findInterval but it makes sure the basepair is both
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

sites$exon = sapply(sites$base,findBoundedInterval,maptable=exons)

introns = data.frame(exons$stop[1:66],exons$start[2:67])
colnames(introns) = c('start','stop')
sites$intron = sapply(sites$base,findBoundedInterval,maptable=introns)

plot(sites$base[Mmu5],sites$exon[Mmu5],type='h',col='blue')
points(sites$base[Mmu5],sites$intron[Mmu5],type='h',col='gray')
points(sites$base[Mmu5],sites$exon[Mmu5],type='h',col='blue')

sites$sumdepth = rowSums(covmat)
exonic = sites$exon > 0
plot(sites$exon[Mmu5 & exonic], sites$sumdepth[Mmu5 & exonic], pch='.')
plot(sites$intron[Mmu5 & !exonic], sites$sumdepth[Mmu5 & !exonic], pch='.')

sql_query = "
select   avg(s.sumdepth) avdepth
from     sites s
where    s.exon > 0
group by s.exon
order by s.exon
;
"
exons$avdepth = unlist(sqldf(sql_query))/length(samplenames)
png('av.depth.by.exon.png',width=600,height=400)
barplot(exons$avdepth,col='black',main='Average depth by exon',xlab='Exon',
        names.arg=1:dim(exons)[1],ylab='Average total depth',
        cex.names=.6,las=3)
dev.off()

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
# ratio of intron av depth to that of the exon to its right
# right as opposed to left or both because this avoids exon 1
# which has the complicated knock in situation
ieratio = introns$avdepth / exons$avdepth[2:dim(exons)[1]]
png('ieratio.png',width=600,height=400)
barplot(ieratio,col='black',main='Intron/exon ratio',
        xlab='Intron (n) / exon (n+1)',ylab='Depth ratio',
        names.arg=1:dim(introns)[1],cex.names=.6,las=3)
dev.off()

# does relationship between exon and depth depend on:
# 1. Q length
# 2. Zygosity
# 3. Sex
# 4. Tissue

exon_avgs = aggregate(covmat, list(sites$exon), mean)
exon_avgs_rel = melt(exon_avgs,id=c("Group.1"))
colnames(exon_avgs_rel) = c("exon","sid","avdepth")
exon_avgs_rel = exon_avgs_rel[exon_avgs_rel$exon > 0,]
exon_avgs_data = merge(exon_avgs_rel, metadata, by="sid")

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

for (q in qs) {
  m = lm(exon_avgs_data$avdepth[het&stri&exon_avgs_data$q==q] ~ exon_avgs_data$exon[het&stri&exon_avgs_data$q==q])
  print(q)
  print(summary(m)$coefficients[2,1])
}

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

dim(covmat)
head(exon_avgs_data)

# plot the exon 1 / intron 1 area in particular

# what's up with intron 49?
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

# do these peaks match predicted polyA sites??
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
polyah = read.table('mmu5.htt.intron1.polyah.txt',skip=3)
polyah = polyah[,c(2,4)]
colnames(polyah) = c("relpos","ldf")
polyah$abspos = polyah$relpos + introns$start[1]

png('intron.1.detail.polyah.png',width=600,height=400)
plot(sites$base[sites$intron==1], sites$sumdepth[sites$intron==1], type='h', col='gray',
     xlab = 'Mmu5 base pair position', ylab = 'Total depth', main = 'Coverage of intron 1')
points(sites$base[sites$exon==1], sites$sumdepth[sites$exon==1], type='h', col='blue')
points(sites$base[sites$exon==2], sites$sumdepth[sites$exon==2], type='h', col='blue')
points(polyah$abspos, polyah$ldf*(1500/max(polyah$ldf)), type='h', col='red')
legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
dev.off()

sum(metadata$zyg == 'WT' & metadata$tiss == nowtiss) # 48

max(ex1to2[,6:12]) # 23
ytop = 25


ytop = 25
ex1to2 = subset(sites, exon %in% c(1,2) | intron == 1)
incl = sites$exon %in% c(1,2) | sites$intron == 1
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
    png(paste('ex1to2.cov',nowtiss,genotype,'png',sep='.'),width=600,height=400)
    gt_depth = rowMeans(covmat[incl, genotypes[[genotype]]])
    plot(ex1to2$base, gt_depth, type='h', col='gray', ylim=c(0,ytop), main=paste('Exon 1 - 2,',genotype,nowtiss),
         xlab = 'Mmu5 base pair position', ylab='Average depth')
    points(ex1to2$base[ex1to2$exon > 0],  gt_depth[ex1to2$exon > 0], type='h', col='blue')
    points(polyah$abspos, polyah$ldf*(ytop/max(polyah$ldf)), type='h', col='red')
    legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
    dev.off()
  }
}


ranges = list()
ranges[['KI_CONSTRUCT']] = sites$contig != '5'
ranges[['EXON1_INTRON1_EXON2']] = sites$exon %in% c(1,2) | sites$intron == 1
ranges[['EXON49_INTRON49_EXON50']] = sites$exon %in% c(49,50) | sites$intron == 49
ranges[['EXON59_INTRON59_EXON60']] = sites$exon %in% c(59,60) | sites$intron == 59

for (rangename in names(ranges)) {
  range = ranges[[rangename]]
  use_sites = sites[range,] 
  ytop = max(rowMeans(covmat[range,]))
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
      png(paste(rangename,nowtiss,genotype,'cov.png',sep='.'),width=600,height=400)
      gt_depth = rowMeans(covmat[range, genotypes[[genotype]]])
      plot(use_sites$base, gt_depth, type='h', col='gray', ylim=c(0,ytop), main=paste(rangename,genotype,nowtiss),
           xlab = 'Mmu5 base pair position', ylab='Average depth')
      points(use_sites$base[use_sites$exon > 0],  gt_depth[use_sites$exon > 0], type='h', col='blue')
      points(polyah$abspos, polyah$ldf*(ytop/max(polyah$ldf)), type='h', col='red')
      legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
      dev.off()
    }
  }
}

# next: look at reads and see if any 


# create a sample-wise slope for 5' degradation of exons
# and a factor for 
for (sid in samplenames) {
  # model of 5' deg based on exons 2:67
  m = lm(avdepth ~ exon, data=subset(exon_avgs_data, sid==sid & exon > 1))
  
}

