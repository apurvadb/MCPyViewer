#!/bin/bash
bwa mem -t 1 -R '@RG\tID:hpv\tSM:hpv\tLB:hpv\tPL:ILLUMINA' -M -t 8 searcHPV_results/Sample_ERR4425693/hg_hpv.fa searcHPV_results/Sample_ERR4425693/assemble//6.20569311/pearOutput//6.20569311.all.fa.cap.contigs > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.sam;
samtools view -@ 1 -bhS searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.sam > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.bam;
samtools sort -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.bam -o searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.sort.bam;
samtools index -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.sort.bam;
rm searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20569311//6.20569311.contig.bam;
bwa mem -t 1 -R '@RG\tID:hpv\tSM:hpv\tLB:hpv\tPL:ILLUMINA' -M -t 8 searcHPV_results/Sample_ERR4425693/hg_hpv.fa searcHPV_results/Sample_ERR4425693/assemble//6.20635172/pearOutput//6.20635172.all.fa.cap.contigs > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.sam;
samtools view -@ 1 -bhS searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.sam > searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.bam;
samtools sort -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.bam -o searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.sort.bam;
samtools index -@ 1 searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.sort.bam;
rm searcHPV_results/Sample_ERR4425693/call_fusion_virus/6.20635172//6.20635172.contig.bam;
