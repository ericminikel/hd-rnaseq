options(stringsAsFactors=FALSE)
require(sqldf)

setwd('c:/sci/026rplcl/analysis/tx/')
files = list.files(pattern="\\.cov$")
files

metadata = read.table('metadata.txt',sep='\t',header=TRUE)

metadata
coverage = read.table(files[1])[,1:5]
colnames(coverage) = c("contig","start","stop","transcript","relpos")
coverage$base = coverage$start + coverage$relpos

for (file in files) {
  sid = substr(file,1,5) # get the sample id, e.g. "MGH01"
  coverage$tempcol = read.table(file)[,6] # get just the depth column. other columns are same in all files.
  colnames(coverage)[dim(coverage)[2]] = sid # assign sample id as name of new col
}
head(coverage)

# look only at HTT transcript
coverage = subset(coverage, transcript == 'NM_002111')

covmat = as.matrix(coverage[,7:(24+7-1)])
rownames(covmat) = coverage$base
samplenames = colnames(covmat)
head(covmat)

sites = coverage[,1:6]

sites$sumdepth = rowSums(covmat)

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


# cat genes_ERCC.gtf.bed | grep NM_002111 | grep -w exon > htt.exons
exons = read.table('htt.exons')[,2:3]
colnames(exons) = c('start','stop')
sites$exon = sapply(sites$base,findBoundedInterval,maptable=exons)
introns = data.frame(exons$stop[1:66],exons$start[2:67])
colnames(introns) = c('start','stop')
sites$intron = sapply(sites$base,findBoundedInterval,maptable=introns)

head(sites)

png('htt.coverage.total.png',width=600,height=400)
plot(sites$base, sites$sumdepth, type='h', col='gray', xlab='base pair position', ylab='total depth, 24 LCL samples',
     main = 'HTT transcript coverage in 24 LCLs')
points(sites$base[sites$exon > 0], sites$sumdepth[sites$exon > 0], type='h', col='blue')
dev.off()


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

introns[1,]
# chr4:3076815-3088665
# retrieved DNA from UCSC Genome Browser > PolyAH
# http://linux1.softberry.com/berry.phtml?topic=polyah&group=programs&subgroup=promoter
# >hg19_dna range=chr4:3076815-3088665 5'pad=0 3'pad=0 strand=
#   Length of sequence-     11851
# 9 potential polyA sites were predicted
# Pos.:    867 LDF-  3.54
# Pos.:   1396 LDF-  3.09
# Pos.:   1409 LDF-  2.68
# Pos.:   2711 LDF-  5.30
# Pos.:   4569 LDF-  1.31
# Pos.:   7328 LDF-  6.73
# Pos.:   8838 LDF-  1.15
# Pos.:   9941 LDF-  0.76
# Pos.:  11597 LDF-  1.47
polyah = read.table('polyah.txt',skip=3)
polyah = polyah[,c(2,4)]
colnames(polyah) = c("relpos","ldf")
polyah$abspos = polyah$relpos + introns$start[1]

ytop = max(sites$sumdepth[sites$intron==1])
png('intron.1.detail.polyah.png',width=600,height=400)
plot(sites$base[sites$intron==1], sites$sumdepth[sites$intron==1], type='h', col='gray',
     xlab = 'chr4 base pair position', ylab = 'Total depth', main = 'Coverage of intron 1')
points(sites$base[sites$exon==1], sites$sumdepth[sites$exon==1], type='h', col='blue')
points(sites$base[sites$exon==2], sites$sumdepth[sites$exon==2], type='h', col='blue')
points(polyah$abspos, polyah$ldf*(ytop/max(polyah$ldf)), type='h', col='red')
legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
dev.off()

# plot hd and control separately
hd = metadata$cag1 > 35
ctrl = metadata$cag1 <= 35
samplesets = list()
samplesets[['hd']] = metadata$cag1 > 35
samplesets[['control']] = metadata$cag1 <= 35
for (sampleset in names(samplesets)) {
  use = samplesets[[sampleset]]
  # intron 1 detail
  depth = rowMeans(covmat[,use])
  ytop = max(depth[sites$intron==1])
  png(paste('intron.1.detail.polyah.',sampleset,'.png',sep=''),width=600,height=400)
  plot(sites$base[sites$intron==1], depth[sites$intron==1], type='h', col='gray',
       xlab = 'chr4 base pair position', ylab = 'mean depth', main = paste('Coverage of intron 1 in ',sampleset,' samples',sep=''))
  points(sites$base[sites$exon==1], depth[sites$exon==1], type='h', col='blue')
  points(sites$base[sites$exon==2], depth[sites$exon==2], type='h', col='blue')
  points(polyah$abspos, polyah$ldf*(ytop/max(polyah$ldf)), type='h', col='red')
  legend('topright',c('intron coverage','exon coverage','polyA site prediction strength'),col=c('gray','blue','red'),lwd=2,pch=NA)
  dev.off()
  # overall HTT
  png(paste('htt.coverage.',sampleset,'.png',sep=''),width=600,height=400)
  plot(sites$base, depth, type='h', col='gray', xlab='base pair position', ylab='mean depth',
       main = paste('HTT transcript coverage in ',sampleset,' samples',sep=''))
  points(sites$base[sites$exon > 0], depth[sites$exon > 0], type='h', col='blue')
  dev.off()
  # 3'UTR detail
  png(paste('3utr.detail.polyah.',sampleset,'.png',sep=''),width=600,height=400)
  plot(sites$base[sites$exon==67], depth[sites$exon==67], type='h', col='blue',
       xlab = 'chr4 base pair position', ylab = 'mean depth', main = paste("Coverage of 3UTR in ",sampleset,' samples',sep=''))
  dev.off()  
}
