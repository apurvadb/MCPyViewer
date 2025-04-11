#!/bin/bash
bwa mem -t 1 -R '@RG\tID:hpv\tSM:hpv\tLB:hpv\tPL:ILLUMINA' -M -t 8 MCV.fasta searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.all.fa.cap.contigs > searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.contigToVirus.sam;
samtools view -@ 1 -bhS searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.contigToVirus.sam >  searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.contigToVirus.bam;
samtools sort -@ 1 searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.contigToVirus.bam -o searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToVirus.sort.bam;
samtools index -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToVirus.sort.bam;
samtools faidx -@ 1 searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.all.fa.cap.contigs;
rm searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.contigToVirus.bam;
bwa mem -t 1 -R '@RG\tID:hpv\tSM:hpv\tLB:hpv\tPL:ILLUMINA' -M -t 8 MCV.fasta searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.all.fa.cap.contigs > searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.contigToVirus.sam;
samtools view -@ 1 -bhS searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.contigToVirus.sam >  searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.contigToVirus.bam;
samtools sort -@ 1 searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.contigToVirus.bam -o searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToVirus.sort.bam;
samtools index -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToVirus.sort.bam;
samtools faidx -@ 1 searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.all.fa.cap.contigs;
rm searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.contigToVirus.bam;
