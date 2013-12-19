mkdir /data/HD/analysis/038ea
cd /data/HD/dataset/mouse_ki_allele_series_rnaseq

find /data/HD/dataset/mouse_ki_allele_series_rnaseq -type "f" -name "*.clipped.fastq.gz" | sort > /data/HD/analysis/038ea/fastq.list

cd /data/HD/analysis/038ea
mkdir analysis
mkdir analysis/1
mkdir analysis/1/fastqc

cat fastq.list | awk '{print "fastqc "$0" -f fastq --nogroup -o /data/HD/analysis/038ea/analysis/1/fastqc/"}' > fastqc.bash
cat fastqc.bash | awk '{print "bsub -q medium -W 12:00 \""$0"\""}' | bash

cd src
bsub.py medium 23:00 "python aggregate-fastqc.py"

ls /data/HD/analysis/038ea/analysis/1/fastqc
# downloaded and messed with files in PostgreSQL

cd /data/HD/dataset/GRCm38.73
bsub.py medium 05:00 "cat Mus_musculus.GRCm38.73.dna.toplevel.fa hdh_q_ki_construct.fa > GRCm38.73.plus.ki.fa"

bsub.py medium 01:00 "samtools faidx GRCm38.73.plus.ki.fa"
bsub.py medium 05:00 "bwa index GRCm38.73.plus.ki.fa"
bsub.py medium 11:00 "gmap_build -d grcm3873ki -k 15 -s none GRCm38.73.plus.ki.fa"

cat fastq.list | awk -F"[./]" '{print $0" "$8}' > fastq.2col
mkdir analysis/1/gsnap/
mkdir analysis/1/gsnap/1

cat fastq.list | awk -F"[_./]" '{print $0" "$8}' | head
# ideally we want filename, name, RGID=five digits, RGNAME=numberL_striatum_Q___HET_F_L8, RG-lib=LB8

cat fastq.2col | awk ' NR%2 == 1 {print "gsnap --format=sam --read-group-id=RGUNKNOWN --read-group-name=$i --read-group-library=LBUNKNOWN --read-group-platform=illumina -N 1 -m 10 -d mm10_e $1 $fastqdir/HD_ESC_pool_${i}_2_fastq.txt | samtools view -Sbh - > $workdir/$i.bam ; samtools sort $i.bam $i.srtd ; 
    samtools index $i.srtd.bam

# http://stackoverflow.com/questions/10220348/putting-two-consecutive-lines-into-one-line-with-perl-awk
cat fastq.list | sed 'N;s/\n/ /' | awk -F"[./]" '{print $0" "$8}' > fastq.3col

cat fastq.list | sed 'N;s/\n/ /' | awk -F"[./]" '{print $0" "$8}' | awk -F"[./_]" '{print $0" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$19" "$20}' > fastq.12col

workdir=/data/HD/analysis/038ea/analysis/1/gsnap/1/
cat fastq.12col | awk -v workdir=$workdir '{print "bsub -q long -m general \"gsnap --gunzip --format\=sam --read-group-id="$4" --read-group-library="$12" --read-group-platform=illumina -N 1 -m 10 -d grcm3873ki "$1" "$2" | samtools view -Sbh - > "workdir $3".bam; samtools sort "workdir $3".bam "workdir $3".srtd; samtools index "workdir $3".srtd.bam\""}' | bash
# awk: warning: escape sequence `\=' treated as plain `='
# probably could have just said --format=sam instead of --format\=sam

#######################################
# 2013-10-22 post-alignment recap

# there are 192 samples
wc -l fastq.12col
# 192 fastq.12col
# alignments are in
cd /data/HD/analysis/038ea/analysis/1/gsnap/1
ls *.bam | wc -l
# 379
# ls *srtd.bam | wc -l
187
# 379 - 187 = 192 non-srtd bams. thus the alignments do appear to have run for all samples.

# now just need to sort them.
# assume that none of the sorting and indexing succeeded.
mkdir junk
mv *.bai junk
mv *.srtd.bam junk

ll | sort -n -k5,5 | head
# most suspicious are these below 5 GB:
-rw-rw----. 1 em476 miket   361337372 Oct  3 17:26 30727_1_712L_cortex_Q80_HET_M_L2.bam
-rw-rw----. 1 em476 miket  1309935161 Oct  2 12:50 21025_1_458L_cortex_Q175_HET_F_L6.bam
-rw-rw----. 1 em476 miket  2428764160 Oct  2 19:36 28269_1_646L_cortex_Q20_HET_M_L2.bam
-rw-rw----. 1 em476 miket  2550628352 Oct  2 19:38 23457_1_522L_cortex_Q111_HET_F_L2.bam
-rw-rw----. 1 em476 miket  2909732864 Oct  2 19:38 514_Liver_Q111_HET_M_L8.bam
-rw-rw----. 1 em476 miket  3203760128 Oct  2 19:35 20927_1_450L_striatum_Q175_HET_M_L1.bam
-rw-rw----. 1 em476 miket  3593273344 Oct  2 19:38 23528_1_527L_striatum_Q111_WT_F_L7.bam

ll | sort -n -k5,5 | tail
# all look good




##########################
# now do sorting and indexing
# test run:
samtools view -h 20914_1_449L_striatum_Q175_WT_M_L1.bam | head -10000 | samtools view -Sbh > test.bam
for bamfile in test.bam
do
  bsub -q short -W 00:20 "workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/1; tmpdir=\$(mktemp -d); trap 'rm -rfv \$tmpdir' EXIT; cp $workdir/$bamfile \$tmpdir; samtools sort \$tmpdir/$bamfile \$tmpdir/$bamfile.srtd; samtools index \$tmpdir/$bamfile.srtd.bam; mv \$tmpdir/$bamfile.srtd.bam $workdir/srtd; mv \$tmpdir/$bamfile.srtd.bam.bai $workdir/srtd; rm \$tmpdir/$bamfile"
done

# tips from Ashok
# tmp_dir=$(env TMPDIR=${outputdir} mktemp -d ) # didn't use this one
# trap 'rm -rfv $tmp_dir' EXIT # did use this


cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1
tmpdir=/scratch/miket/tmpEric
cd $workdir
mkdir srtd
for bamfile in *.bam
do
  bsub -q medium -W 20:00 -m general "workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/1; tmpdir=\$(mktemp -d); trap 'rm -rfv \$tmpdir' EXIT;  cp $workdir/$bamfile \$tmpdir; samtools sort \$tmpdir/$bamfile \$tmpdir/$bamfile.srtd; samtools index \$tmpdir/$bamfile.srtd.bam; mv \$tmpdir/$bamfile.srtd.bam $workdir/srtd; mv \$tmpdir/$bamfile.srtd.bam.bai $workdir/srtd; rm \$tmpdir/$bamfile"
done

########################
# Check that FASTQ and BAM sizes match

# check that FASTQ file sizes are concordant
cd /data/HD/analysis/evm/038ea
cat fastq.12col | awk '{print "du -sb "$1}' | bash | awk '{print $1}' > fastq.1.filesize
cat fastq.12col | awk '{print "du -sb "$2}' | bash | awk '{print $1}' > fastq.2.filesize
paste fastq.1.filesize fastq.2.filesize | awk '{print $1 - $2}' | sort -n | head
paste fastq.1.filesize fastq.2.filesize | awk '{print $1 - $2}' | sort -n | tail
# range -100MB to +80MB


cd /data/HD/analysis/evm/038ea
cat fastq.12col | awk '{print "du -sb "$1}' | bash > fastq.1.filesize
cat fastq.12col | awk '{print "du -sb "$2}' | bash > fastq.2.filesize
paste fastq.1.filesize fastq.2.filesize | awk '{print $1 - $3" "$2" "$4}' | sort -n | tail
# confirmed

# given that the FASTQ themselves are like 3 GB, this is just a few percent.
# can plausibly be compression differences and not differences in
# number of lines or something like that.

# given that the two FASTQs are close enough,
# we'll look at FASTQ size vs. BAM size.

cat fastq.12col | awk '{print "du -sb "$1}' | bash > fastq.1.filesize

ls analysis/1/gsnap/1/srtd/*.bam | awk '{print "du -sb "$1}' | bash > bam.filesize

# switch to R
R
fastq = read.table('fastq.1.filesize',header=FALSE)
bam = read.table('bam.filesize',header=FALSE)
bam$name = substr(bam$V2,25,55)
fastq$name = substr(fastq$V2,71,101)
bam$bamsize=bam$V1
fastq$fastqsize=fastq$V1
merged = merge(fastq,bam,"name")[,c("name","bamsize","fastqsize")]
png('bam.vs.fastq.size.png',width=500,height=500)
plot(merged$bamsize,merged$fastqsize,pch=19)
dev.off()
png('bam.vs.fastq.size.with.labels.png',width=500,height=500)
plot(merged$bamsize,merged$fastqsize,pch=19)
text(merged$bamsize,merged$fastqsize,labels=merged$name,pos=4,cex=.8)
dev.off()
# examining these plots there are 7 clear outliers where alignment
# must just not have finished.
head(sort(merged$fastqsize/merged$bamsize))
#[1] 0.3514891 0.3537249 0.3537428 0.3538352 0.3545493 0.3546114
tail(sort(merged$fastqsize/merged$bamsize),n=10)
# [1] 0.3820659 0.3841930 0.3924874 0.4298696 0.7758920 0.8830711 0.9025634
# [8] 1.1106061 2.1686593 6.2669840
badbams = merged$name[merged$fastqsize/merged$bamsize > .40] # 7 entries
write.table(badbams,'badbams.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
q()
n


##########################
# Re-align the 7 failed BAMs

cat badbams.txt
grep -f badbams.txt fastq.12col | wc -l 
# 7
grep -f badbams.txt fastq.12col | head -1

# realign just those 7
# first need to get base filename out of full path, for temporary directory access. add 2 more columns.
cat fastq.12col | awk -F'[/ ]' '{print $7" "$14" "$0}' > fastq.14col
# not worth my time bc the first line isn't one of the 7 that need realignment anyway.
mkdir analysis/1/gsnap/2/
workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/2/
grep -f badbams.txt fastq.14col | awk -v workdir=$workdir '{print "bsub -q big -m general \"workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/2; tmpdir=\\\$(mktemp -d); trap '\''rm -rfv \$tmpdir'\'' EXIT; cp "$3" \\\$tmpdir; cp "$4" \\\$tmpdir; gsnap --gunzip --format\=sam --read-group-id="$6" --read-group-library="$14" --read-group-platform=illumina -N 1 -m 10 -d grcm3873ki \\\$tmpdir/"$1" \\\$tmpdir/"$2" | samtools view -Sbh - > \\\$tmpdir/"$5".bam; cp \\\$tmpdir/"$5".bam "workdir";\""}' | bash

#awk: warning: escape sequence `\$' treated as plain `$'
#awk: warning: escape sequence `\=' treated as plain `='


###########################
# Check that the 7 failed BAMs are fixed

cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/2

[em476@eris1n2 2]$ ll
total 45045384
-rw-rw----. 1 em476 miket 8054421376 Nov  9 22:59 20927_1_450L_striatum_Q175_HET_M_L1.bam
-rw-rw----. 1 em476 miket 8194416037 Nov  9 09:29 21025_1_458L_cortex_Q175_HET_F_L6.bam
-rw-rw----. 1 em476 miket 6395104031 Nov 10 06:15 23457_1_522L_cortex_Q111_HET_F_L2.bam
-rw-rw----. 1 em476 miket 7731821382 Nov  8 19:30 23528_1_527L_striatum_Q111_WT_F_L7.bam
-rw-rw----. 1 em476 miket 7535565600 Nov  9 16:41 28269_1_646L_cortex_Q20_HET_M_L2.bam
-rw-rw----. 1 em476 miket 8167264705 Nov  8 21:42 28301_1_648L_striatum_Q20_HET_M_L6.bam
[em476@eris1n2 2]$ for file in *
> do
> ll ../1/$file
> done
-rwxr-x---. 1 em476 hdmiket 3203760128 Oct  2 19:35 ../1/20927_1_450L_striatum_Q175_HET_M_L1.bam
-rwxr-x---. 1 em476 hdmiket 1309935161 Oct  2 12:50 ../1/21025_1_458L_cortex_Q175_HET_F_L6.bam
-rwxr-x---. 1 em476 hdmiket 2550628352 Oct  2 19:38 ../1/23457_1_522L_cortex_Q111_HET_F_L2.bam
-rwxr-x---. 1 em476 hdmiket 3593273344 Oct  2 19:38 ../1/23528_1_527L_striatum_Q111_WT_F_L7.bam
-rwxr-x---. 1 em476 hdmiket 2428764160 Oct  2 19:36 ../1/28269_1_646L_cortex_Q20_HET_M_L2.bam
-rwxr-x---. 1 em476 hdmiket 6971037326 Oct  3 20:55 ../1/28301_1_648L_striatum_Q20_HET_M_L6.bam



cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1
tmpdir=/scratch/miket/tmpEric
cd $workdir
mkdir srtd
for bamfile in *.bam
do
  bsub -q medium -W 20:00 -m general "workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/1; tmpdir=\$(mktemp -d); trap 'rm -rfv \$tmpdir' EXIT;  cp $workdir/$bamfile \$tmpdir; samtools sort \$tmpdir/$bamfile \$tmpdir/$bamfile.srtd; samtools index \$tmpdir/$bamfile.srtd.bam; mv \$tmpdir/$bamfile.srtd.bam $workdir/srtd; mv \$tmpdir/$bamfile.srtd.bam.bai $workdir/srtd; rm \$tmpdir/$bamfile"
done





# check now with the re-aligned BAMs in there

cd /data/HD/analysis/evm/038ea/
cat fastq.12col | awk '{print "du -sb "$1}' | bash > fastq.1.filesize
ls *.bam | awk '{print "du -sb "$1}' | bash > bam.filesize

# switch to R
R
fastq = read.table('/data/HD/analysis/evm/038ea/fastq.1.filesize',header=FALSE)
bam = read.table('bam.filesize',header=FALSE)
bam$name = substr(bam$V2,1,lapply(bam$V2,length)-4)
fastq$name = substr(fastq$V2,71,101)
bam$bamsize=bam$V1
fastq$fastqsize=fastq$V1
merged = merge(fastq,bam,"name")[,c("name","bamsize","fastqsize")]
png('bam.vs.fastq.size.png',width=500,height=500)
plot(merged$bamsize,merged$fastqsize,pch=19)
dev.off()
png('bam.vs.fastq.size.with.labels.png',width=500,height=500)
plot(merged$bamsize,merged$fastqsize,pch=19)
text(merged$bamsize,merged$fastqsize,labels=merged$name,pos=4,cex=.8)
dev.off()
# examining these plots there are 7 clear outliers where alignment
# must just not have finished.
head(sort(merged$fastqsize/merged$bamsize))
#[1] 0.3514891 0.3537249 0.3537428 0.3538352 0.3545493 0.3546114
tail(sort(merged$fastqsize/merged$bamsize),n=10)
# [1] 0.3820659 0.3841930 0.3924874 0.4298696 0.7758920 0.8830711 0.9025634
# [8] 1.1106061 2.1686593 6.2669840
badbams = merged$name[merged$fastqsize/merged$bamsize > .40] # 7 entries
write.table(badbams,'badbams.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
q()
n

cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/2
workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/2
tmpdir=/scratch/miket/tmpEric
cd $workdir
mkdir srtd
for bamfile in *.bam
do
  bsub -q medium -W 20:00 -m general "workdir=/data/HD/analysis/evm/038ea/analysis/1/gsnap/2; tmpdir=\$(mktemp -d); trap 'rm -rfv \$tmpdir' EXIT;  cp $workdir/$bamfile \$tmpdir; samtools sort \$tmpdir/$bamfile \$tmpdir/$bamfile.srtd; samtools index \$tmpdir/$bamfile.srtd.bam; mv \$tmpdir/$bamfile.srtd.bam $workdir/srtd; mv \$tmpdir/$bamfile.srtd.bam.bai $workdir/srtd; rm \$tmpdir/$bamfile"
done


cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1
mkdir bad
for file in ../2/srtd/*.bam; do ll $file; done
mv srtd/20927_1_450L_striatum_Q175_HET_M_L1.bam.srtd.bam      bad/
mv srtd/21025_1_458L_cortex_Q175_HET_F_L6.bam.srtd.bam        bad/
mv srtd/23457_1_522L_cortex_Q111_HET_F_L2.bam.srtd.bam        bad/
mv srtd/23528_1_527L_striatum_Q111_WT_F_L7.bam.srtd.bam       bad/
mv srtd/28269_1_646L_cortex_Q20_HET_M_L2.bam.srtd.bam         bad/
mv srtd/28301_1_648L_striatum_Q20_HET_M_L6.bam.srtd.bam       bad/
mv srtd/30727_1_712L_cortex_Q80_HET_M_L2.bam.srtd.bam         bad/


mv ../2/srtd/20927_1_450L_striatum_Q175_HET_M_L1.bam.srtd.bam  ./srtd/
mv ../2/srtd/21025_1_458L_cortex_Q175_HET_F_L6.bam.srtd.bam    ./srtd/
mv ../2/srtd/23457_1_522L_cortex_Q111_HET_F_L2.bam.srtd.bam    ./srtd/
mv ../2/srtd/23528_1_527L_striatum_Q111_WT_F_L7.bam.srtd.bam   ./srtd/
mv ../2/srtd/28269_1_646L_cortex_Q20_HET_M_L2.bam.srtd.bam     ./srtd/
mv ../2/srtd/28301_1_648L_striatum_Q20_HET_M_L6.bam.srtd.bam   ./srtd/
mv ../2/srtd/30727_1_712L_cortex_Q80_HET_M_L2.bam.srtd.bam     ./srtd/


mv ../2/srtd/20927_1_450L_striatum_Q175_HET_M_L1.bam.srtd.bam.bai ./srtd/
mv ../2/srtd/21025_1_458L_cortex_Q175_HET_F_L6.bam.srtd.bam.bai ./srtd/
mv ../2/srtd/23457_1_522L_cortex_Q111_HET_F_L2.bam.srtd.bam.bai ./srtd/
mv ../2/srtd/23528_1_527L_striatum_Q111_WT_F_L7.bam.srtd.bam.bai ./srtd/
mv ../2/srtd/28269_1_646L_cortex_Q20_HET_M_L2.bam.srtd.bam.bai ./srtd/
mv ../2/srtd/28301_1_648L_striatum_Q20_HET_M_L6.bam.srtd.bam.bai ./srtd/
mv ../2/srtd/30727_1_712L_cortex_Q80_HET_M_L2.bam.srtd.bam.bai ./srtd/



cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1/srtd
du -sb *.bam > bam.filesize
cat /data/HD/analysis/evm/038ea/fastq.12col | awk '{print "du -sb "$1}' | bash > fastq.filesize
# switch to R
R
fastq = read.table('fastq.filesize',header=FALSE)
bam = read.table('bam.filesize',header=FALSE)
bam$name = substr(bam$V2,1,21)
fastq$name = substr(fastq$V2,71,91)
bam$bamsize=bam$V1
fastq$fastqsize=fastq$V1
merged = merge(fastq,bam,"name")[,c("name","bamsize","fastqsize")]
dim(merged) # 192 3

png('bam.vs.fastq.size.png',width=500,height=500)
plot(merged$bamsize,merged$fastqsize,pch=19)
dev.off()
png('bam.vs.fastq.size.with.labels.png',width=500,height=500)
plot(merged$bamsize,merged$fastqsize,pch=19)
text(merged$bamsize,merged$fastqsize,labels=merged$name,pos=4,cex=.8)
dev.off()

write.table(merged,"merged.size.txt",sep="\t",col.names=TRUE,row.names=FALSE)

badbams = merged$name[merged$fastqsize/merged$bamsize > .41]
badbams
# [1] "514_Liver_Q111_HET_M_" # .68 ratio

q()
n

# there are a few other marginal cases based on the plot
# also this time there are several BAMs much larger than the FASTQ, all from liver
# I didn't see them the first time I did this when I had accidentally left out liver
# perhaps those just have high information entropy???? due to contamination or something?
# 

# now need to quick align #514

cd /data/HD/analysis/evm/038ea/tempdata
rm *
bsub.py medium 08:00 "gunzip -c /data/HD/dataset/mouse_ki_allele_series_rnaseq/20130424-H22819.FASTQs/514_Liver_Q111_HET_M_L8.LB9_1.clipped.fastq.gz | split -d -a 3 -l 1000000 - tmp1.list"
bsub.py medium 08:00 "gunzip -c /data/HD/dataset/mouse_ki_allele_series_rnaseq/20130424-H22819.FASTQs/514_Liver_Q111_HET_M_L8.LB9_2.clipped.fastq.gz | split -d -a 3 -l 1000000 - tmp2.list"

for i in {000..177}
do
  bsub.py short 00:50 "gsnap --format=sam --read-group-id=514 --read-group-library=LB9 --read-group-platform=illumina -N 1 -m 10 -d grcm3873ki tmp1.list$i tmp2.list$i | samtools view -Sbh - > tmp.list$i.bam"
done

# now to re-merge
cd /data/HD/analysis/evm/038ea/tempdata
bsub.py medium 23:00 "samtools merge 514_Liver_Q111_HET_M_L8.LB9.bam `ls tmp.list*.bam | tr '\n' ' '`"

# and then sort
bsub.py medium 08:00 "rm /data/HD/analysis/evm/038ea/analysis/1/srtd/514_Liver_Q111_HET_M_L8.LB9.bam; samtools sort /data/HD/analysis/evm/038ea/tempdata/514_Liver_Q111_HET_M_L8.LB9.bam /data/HD/analysis/evm/038ea/analysis/1/srtd/514_Liver_Q111_HET_M_L8.LB9.bam.srtd; samtools index /data/HD/analysis/evm/038ea/analysis/1/srtd/514_Liver_Q111_HET_M_L8.LB9.bam.srtd.bam"

# re-ran the R script above. the file is STILL too small, no change. must be corrupted in some other way and not just
# that the job didn't finish.

# see http://unix.stackexchange.com/questions/10421/output-from-ls-has-newlines-but-displays-on-a-single-line-why
# ls -C forces ls to use columns, but it still has SOME newlines
# difficult to use sed for this. see http://stackoverflow.com/questions/1251999/sed-how-can-i-replace-a-newline-n

# now get ready for cufflinks
# cufflinks runs separately on each sample
# the GTF file is:
/data/HD/public/GRCm38.73/Mus_musculus.GRCm38.73.gtf

cd /data/HD/analysis/evm/038ea/analysis/1
mkdir cufflinks
mkdir cufflinks/1
cat ../../fastq.12col | awk '{print "cufflinks -o cufflinks/1/"$3" -G /data/HD/public/GRCm38.73/Mus_musculus.GRCm38.73.gtf gsnap/1/srtd/"$3".bam.srtd.bam"}' > cufflinks.bash
# examine bash script
cat cufflinks.bash | awk '{print "bsub.py medium 23:00 \""$0"\""}' | bash

# adapted from https://gist.github.com/ericminikel/7533289
cd /data/HD/analysis/evm/038ea/analysis/1/cufflinks/1
mkdir ../combined
# prep header rows
cat 847_Liver_Q92_HET_F_L5/isoforms.fpkm_tracking | head -1 | awk '{print "SAMPLENAME\t"$0}' > ../combined/all.isoforms.fpkm_tracking
cat 847_Liver_Q92_HET_F_L5/genes.fpkm_tracking    | head -1 | awk '{print "SAMPLENAME\t"$0}' > ../combined/all.genes.fpkm_tracking
# combine files
for subdir in *
do
    # combine isoforms data into one file, adding a column for sample name
    cat $subdir/isoforms.fpkm_tracking | tail -n +2 | awk -v subdir=$subdir '{print subdir"\t"$0}' >> ../combined/all.isoforms.fpkm_tracking
    # combine gene symbol data into one file, adding a column for sample name
    cat $subdir/genes.fpkm_tracking | tail -n +2 | awk -v subdir=$subdir '{print subdir"\t"$0}' >> ../combined/all.genes.fpkm_tracking
done

for subdir in *
do
 echo $subdir
done
# comparison is accomplished with cuffdiff
# the way to parallelize is prob with -p, NOT by splitting the GTF into separate files
# see: http://www.biostars.org/p/11001/ 

# now for the bedtools part

cd /data/HD/analysis/evm/038ea/analysis/1/cufflinks/1/23541_1_528L_striatum_Q111_HET_F_L5
cat isoforms.fpkm_tracking | grep Htt
# this one is longest: 5:34761743-34912534
# and it looks good in UCSC genome browser too
cd /data/HD/public/GRCm38.73
cat Hdh_Q_ki_construct      799     3846369812      50      51

cat > htt.bed
5	34761743	34912534
Hdh_Q_ki_construct	0	799
^D

# bedtools calculations
cd /data/HD/analysis/evm/038ea/analysis/1/
mkdir gsnap/1/httcov/
mkdir gsnap/1/httcov/1
cat htt.bed
bamfile=gsnap/1/srtd/test.bam.srtd.bam # first try on my small "test" file
coverageBed -abam $bamfile -b htt.bed -split -d > $bamfile.htt.cov
for bamfile in gsnap/1/srtd/*.bam
do
    bsub -m general -q medium -W 20:00 "cd /data/HD/analysis/evm/038ea/analysis/1; coverageBed -abam $bamfile -b htt.bed -split -d > $bamfile.htt.cov"
done

mv gsnap/1/srtd/*.cov gsnap/1/httcov/1/

cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1/httcov/1
mkdir ../combined
mkdir ../1col
# create a file to hold just the contig and base columns
cat > ../combined/locations.txt 
contig base
^D
cat test.bam.srtd.bam.htt.cov | awk '{print $1" "$2+$4}' >> ../combined/locations.txt
# now create a 1-column version of all the coverage files
for covfile in *
do
  #echo -n $covfile" " >> ../combined/header_row.txt
  echo $covfile | awk -F"." '{print $1}' > ../1col/$covfile.1col # add a header row with the sample name (and some extraneous stuff)
  cat $covfile | awk '{print $5}' >> ../1col/$covfile.1col # add just the depth column ($5)
done
# paste them all together
paste ../combined/locations.txt ../1col/* > ../combined/allsamples.cov # paste all the samples together into a matrix
# remove clutterfiles
rm ../1col/*

# final file:
/data/HD/analysis/evm/038ea/analysis/1/gsnap/1/httcov/combined/allsamples.cov
# Use RSEM to get TPM values

rsem-prepare-reference --help
rsem-calculate-expression --help

cd ~/bin
wget http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.7.tar.gz
gunzip rsem-1.2.7.tar.gz
tar -xvf rsem-1.2.7.tar
cd rsem-1.2.7
make
# then add ~/bin/rsem-1.2.7 to PATH


