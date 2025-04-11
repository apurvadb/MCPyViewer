#!/bin/bash
bwa mem -t 1 -R '@RG\tID:hpv\tSM:hpv\tLB:hpv\tPL:ILLUMINA' -M -t 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.all.fa.cap.contigs > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.sam;
samtools view -@ 1 -bhS searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.sam > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.bam;
samtools sort -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.bam -o searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.sort.bam;
samtools index -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.sort.bam;
rm searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contigToGenome.bam;
bwa mem -t 1 -R '@RG\tID:hpv\tSM:hpv\tLB:hpv\tPL:ILLUMINA' -M -t 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.all.fa.cap.contigs > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.sam;
samtools view -@ 1 -bhS searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.sam > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.bam;
samtools sort -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.bam -o searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.sort.bam;
samtools index -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.sort.bam;
rm searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contigToGenome.bam;
