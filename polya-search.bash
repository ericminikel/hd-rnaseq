for bamfile in *.bam
do
  samtools view $bamfile 5:34761744-34912534 | awk '$10 ~ /(AAAAAAA$)/ && $6 ~ /S$/ {print $4}' >> polya.try4.txt
done
