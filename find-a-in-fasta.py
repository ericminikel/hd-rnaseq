# find all the genomic sites in the Htt transcript that contain at least 7
# consecutive A.  these are potential false positive sites for my other script,
# polya-search.bash, which is trying to find reads which appear to be 
# polyadenylated

import re

with open('c:/sci/034ea/analysis/1/mmu5.htt.fasta.txt',mode='rb') as f:
    fastatxt = f.read()

fastasplit = fastatxt.split('\n')
header = fastasplit[0]
fastatxt = ''.join(fastasplit[1:])

htt_start = 34761739
# print fastatxt[:100]

for match in re.finditer('AAAAAAA',fastatxt):
    print htt_start + match.start()
