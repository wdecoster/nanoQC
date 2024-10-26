set -ev

if [ -d "nanotest" ]; then
    echo "nanotest already cloned"
else
    git clone https://github.com/wdecoster/nanotest.git
fi

nanoQC -h
nanoQC nanotest/reads.fastq.gz
nanoQC nanotest/reads.fastq.gz --minlen 500
