bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
picard CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.primary_assembly.fa O=Homo_sapiens.GRCh38.dna.primary_assembly.dict
