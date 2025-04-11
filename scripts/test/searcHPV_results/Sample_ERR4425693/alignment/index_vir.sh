bwa index MCV.fasta
samtools faidx MCV.fasta
picard CreateSequenceDictionary R=MCV.fasta O=MCV.dictsta
