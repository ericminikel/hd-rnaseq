# polya-search.bash
# try to pull out reads that appear to be polyadenylated from the knock-in allelic series BAMs for Htt transcript

###########
# WHAT WORKED THE BEST:
#
cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1/srtd
for bamfile in *.bam
do
  samtools view $bamfile 5:34761744-34912534 | awk '$10 ~ /(AAAAAAA$)/ && $6 ~ /S$/ {print $4}' >> polya.try4.txt
done
############


# everything else below is stuff I tried that worked less well:




# Htt: 5:34761744-34912534

$ samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam 5:34761744-34912534 | grep AAAAAAAAAA | less
HWI-ST1131:130420:D24VUACXX:4:1106:3136:100195  99      5       34884962        40      33M     =       34885093        181     AACACATTCTGTTTCTTTAAAAAAAAAAAAAAA       ?@@FDADDH
HWI-ST1131:130420:D24VUACXX:4:2106:20882:81538  163     5       34884972        40      23M     =       34885063        141     GTTTCTTTAAAAAAAAAAAAAAA @BBFFFFFHHHHHJJJJIJJJJJ R
HWI-ST1131:130420:D24VUACXX:4:2304:6210:2773    99      5       34884972        40      24M1S   =       34885063        141     GTTTCTTTAAAAAAAAAAAAAAAAA       CBCFFFFFHGHHHJJJJ
HWI-ST1131:130420:D24VUACXX:4:1214:14321:37721  83      5       34884978        40      50M     =       34884824        -204    TTAAAAAAAAAAAAAAAAGTCCAATGTTACAGGAGCCCTCAAGACCCTG
HWI-ST1131:130420:D24VUACXX:4:1216:10627:17487  147     5       34884982        40      38M     =       34884801        -219    AAAAAAAAAAAAAAGTCCAATGTTACAGGAGCCCTCAA  DFHJGIJJI
HWI-ST1131:130420:D24VUACXX:4:1212:5876:21691   147     5       34912487        40      33M9S   =       34912393        -136    TTGTACAGTTGAATTAAAAAAGCGATTAAACCAAAAAAAAAA      ?


# HWI-ST1131:130420:D24VUACXX:4:1303:13728:29704  339     1       4748220 0       22M     =       4721991 -26251  AAAAAAAAAAAAAAAAAAAAAC  B=8?DEFEDF?8FFDD=B??=8  RG:Z:28288      M

zcat /data/HD/dataset/mouse_ki_allele_series_rnaseq/20130420-H22816.FASTQs/28288_1_647L_striatum_Q20_WT_M_L4.LB15_1.clipped.fastq.gz | grep AAAAAAAAAAAAAAAAAAAA | less
# sure enough, they're all of different length
/data/HD/dataset/mouse_ki_allele_series_rnaseq/20130420-H22816.FASTQs/28288_1_647L_striatum_Q20_WT_M_L4.LB15_2.clipped.fastq.gz 

cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1/srtd
bsub.py medium 08:00 "samtools view samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam | grep AAAAAAAAAAAAAAA > /data/HD/analysis/evm/038ea/analysis/1/28288.A15.sam"

samtools view samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam | grep AAAAAAAAAAAAAAA | 



# http://eqtl.uchicago.edu/RNA_Seq_data/Software/polyA/readme

# Htt + 1000 bases each side: 5:34760744-34913534
$ samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam 5:34760744-34913534 | egrep "(AAAAAAAAAA|TTTTTTTTTT)" | less
HWI-ST1131:130420:D24VUACXX:4:1106:3136:100195  99      5       34884962        40      33M     =       34885093        181     AACACATTCTGTTTCTTTAAAAAAAAAAAAAAA       ?@@FDADDH
HWI-ST1131:130420:D24VUACXX:4:2106:20882:81538  163     5       34884972        40      23M     =       34885063        141     GTTTCTTTAAAAAAAAAAAAAAA @BBFFFFFHHHHHJJJJIJJJJJ R
HWI-ST1131:130420:D24VUACXX:4:2304:6210:2773    99      5       34884972        40      24M1S   =       34885063        141     GTTTCTTTAAAAAAAAAAAAAAAAA       CBCFFFFFHGHHHJJJJ
HWI-ST1131:130420:D24VUACXX:4:1214:14321:37721  83      5       34884978        40      50M     =       34884824        -204    TTAAAAAAAAAAAAAAAAGTCCAATGTTACAGGAGCCCTCAAGACCCTG
HWI-ST1131:130420:D24VUACXX:4:1216:10627:17487  147     5       34884982        40      38M     =       34884801        -219    AAAAAAAAAAAAAAGTCCAATGTTACAGGAGCCCTCAA  DFHJGIJJI
HWI-ST1131:130420:D24VUACXX:4:2211:11329:80224  99      5       34893728        40      38M     =       34893840        162     CTAGCCTGATACCATGCTTTTTTTTTTTTTTTTTTTTT  CCCFFFFFG
HWI-ST1131:130420:D24VUACXX:4:2110:18589:44330  147     5       34893745        40      36M     =       34893444        -337    TTTTTTTTTTTTTTTTTTTTTAAAGTCTGTCTCCGG    DDDDDDBFG
HWI-ST1131:130420:D24VUACXX:4:1308:2014:100092  83      5       34898684        40      50M     =       34896933        -1801   TTTTTTTTTTAATAAGAAATTTGAGATTTCTTCTCAGTCACTTATTTGG
HWI-ST1131:130420:D24VUACXX:4:1212:5876:21691   147     5       34912487        40      33M9S   =       34912393        -136    TTGTACAGTTGAATTAAAAAAGCGATTAAACCAAAAAAAAAA      ?


# go through all reads mapped to Htt+1000 and grep for > 9 A or > 9 T, then require that CIGAR string (col 6) have an S, and if so, grab the base where it aligned (col 4)
# http://www.pement.org/awk/awk1line.txt

samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam 5:34760744-34913534 | egrep "(AAAAAAAAAA|TTTTTTTTTT)" | awk '$6 ~ /S/ {print $4}'

for bamfile in *.bam
do
  samtools view $bamfile 5:34760744-34913534 | egrep "(AAAAAAAAAA|TTTTTTTTTT)" | awk '$6 ~ /S/ {print $4}' >> polya.try1.txt
done
wc -l polya.try1.txt # 539

for bamfile in *.bam
do
  samtools view $bamfile 5:34760744-34913534 | egrep "AAAAAAAAAAAAAAA" | awk '$6 ~ /S/ {print $4}' >> polya.try2.txt
done
wc -l polya.try2.txt # 186

for bamfile in *.bam
do
  samtools view $bamfile 5:34760744-34913534 | awk '$6 ~ /S/ && $10 ~ /AAAAAAAAAA/ {print $4}' >> polya.try3.txt
done
wc -l polya.try3.txt # 


cat ../../../28288.A15.sam | awk '{ if ($7 == "=") print $3" "$8; else print $7" "$8}' > ../../../28288.polya.mates.txt

# > 19 A in trimmed fastq 
zcat /data/HD/dataset/mouse_ki_allele_series_rnaseq/20130420-H22816.FASTQs/28288_1_647L_striatum_Q20_WT_M_L4.LB15_1.clipped.fastq.gz | grep AAAAAAAAAAAAAAAAAAAA | wc -l
# 12331
# equiv in untrimmed
cd /data/HD/dataset/mouse_ki_allele_series_rnaseq/20130420-H22816.FASTQs
zcat 28288_1_647L_striatum_Q20_WT_M_L4.LB15_1.fastq.gz | grep AAAAAAAAAAAAAAAAAAAA | wc -l
# 15162

# and mates?
zcat /data/HD/dataset/mouse_ki_allele_series_rnaseq/20130420-H22816.FASTQs/28288_1_647L_striatum_Q20_WT_M_L4.LB15_2.clipped.fastq.gz | grep AAAAAAAAAAAAAAAAAAAA | wc -l
# 15205

# therefore based on trimmed reads there ought to be 15205 + 12331 = ~27.5K reads with 20+ A
# and indeed we have that many in BAMs...
samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam | grep AAAAAAAAAAAAAAAAAAAA | wc -l
# 28120

samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam 5:34760744-34913534 | awk '$10 ~ /(^AAAAA|AAAAA$)/ && $6 ~ /S$/ {print $0}' | less

# http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr5:34761740-34912521&hgsid=356085703&knownGene=pack&hgFind.matches=uc008xdc.2,



# are there chimeric reads between construct and elsewhere??
$ samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam Hdh_Q_ki_construct:1-799 | awk '$7 != "=" {print $0}' | less # nothing
$ samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam Hdh_Q_ki_construct:1-799 | less
HWI-ST1131:130420:D24VUACXX:4:2214:14055:92976  419     Hdh_Q_ki_construct      1       3       50M     =       118     167     GCCCAAGATGGCTGAGCGCC
HWI-ST1131:130420:D24VUACXX:4:1107:5240:69691   355     Hdh_Q_ki_construct      3       3       50M     =       113     160     CCAAGATGGCTGAGCGCCTT
HWI-ST1131:130420:D24VUACXX:4:1107:19035:87068  355     Hdh_Q_ki_construct      3       3       50M     =       72      119     CCAAGATGGCTGAGCGCCTT
HWI-ST1131:130420:D24VUACXX:4:1313:6559:91568   419     Hdh_Q_ki_construct      3       3       50M     =       81      128     CCAAGATGGCTGAGCGCCTT
HWI-ST1131:130420:D24VUACXX:4:2103:2709:6007    419     Hdh_Q_ki_construct      4       3       50M     =       65      111     CAAGATGGCTGAGCGCCTTG
HWI-ST1131:130420:D24VUACXX:4:2112:15910:74749  419     Hdh_Q_ki_construct      20      3       50M     =       108     138     CTTGGTTCCGCTTCTGCCTG
HWI-ST1131:130420:D24VUACXX:4:2212:13297:60796  419     Hdh_Q_ki_construct      53      3       50M     =       119     116     CATTCATTGCCTTGCTGCTA
HWI-ST1131:130420:D24VUACXX:4:2103:2709:6007    339     Hdh_Q_ki_construct      65      3       49M1S   =       4       -111    TGCTGCTAAGTGGCGCCGCG
HWI-ST1131:130420:D24VUACXX:4:1107:19035:87068  403     Hdh_Q_ki_construct      72      3       50M     =       3       -119    AAGTGGCGCCGCGTAGTGCC
HWI-ST1131:130420:D24VUACXX:4:1313:6559:91568   339     Hdh_Q_ki_construct      81      3       50M     =       3       -128    CGCGTAGTGCCAGTAGGCTC
HWI-ST1131:130420:D24VUACXX:4:2112:15910:74749  339     Hdh_Q_ki_construct      108     3       49M1S   =       20      -138    TCAGGGTCTGTCCCATCGGG
HWI-ST1131:130420:D24VUACXX:4:1107:5240:69691   403     Hdh_Q_ki_construct      113     3       50M     =       3       -160    GTCTGTCCCATCGGGCAGGA
HWI-ST1131:130420:D24VUACXX:4:2214:14055:92976  339     Hdh_Q_ki_construct      118     3       50M     =       1       -167    TCCCATCGGGCAGGAAGCCG
HWI-ST1131:130420:D24VUACXX:4:2212:13297:60796  339     Hdh_Q_ki_construct      119     3       50M     =       53      -116    CCCATCGGGCAGGAAGCCGT

for bamfile in *.bam
do
  samtools view $bamfile  Hdh_Q_ki_construct:1-799 | awk '$7 != "=" {print $0}' >> construct-chimeric-reads.txt
done

bsub.py medium 05:00 "samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam -f 0x100 -c > 28288.0x100.readcount.txt"
bsub.py medium 05:00 "samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam -f 0x04 -c > 28288.0x04.readcount.txt"
bsub.py medium 05:00 "samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam -c > 28288.total.readcount.txt"


samtools view 28288_1_647L_striatum_Q20_WT_M_L4.bam.srtd.bam


cd /data/HD/public/humanbodymap2.0/bams
samtools view GRCh37.HumanBodyMap.adipose.1.bam 4:3076408-3245687 | awk '$10 ~ /(^AAAAA|AAAAA$)/ && $6 ~ /S$/ {print $0}' | less

# see which libs the 7 polyA reads in intron 1 came from
cd /data/HD/analysis/evm/038ea/analysis/1/gsnap/1/srtd
for bamfile in *.bam
do
  samtools view $bamfile 5:34762746-34763234 | awk -v bamfile=$bamfile '$10 ~ /AAAAAAA$/ && $6 ~ /S$/ {print bamfile" "$0}'  >> temp7.txt
done
less temp7.txt

