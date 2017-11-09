set -ev

git clone https://github.com/wdecoster/nanotest.git

NanoQC -h
NanoQC nanotest/reads.fastq.gz
