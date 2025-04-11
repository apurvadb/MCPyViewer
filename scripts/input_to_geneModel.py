#!/usr/bin/env python
# coding: utf-8

import argparse
import pysam
import os
import pickle
from collections import Counter
import pandas as pd
import re
import matplotlib.pyplot as plt
import csv


# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="MCV gene analysis and sequence extraction tool")
    parser.add_argument('-w', '--workdir', required=True, help="Working directory path")
    parser.add_argument('-g', '--gtf', required=True, help="Path to the GTF genome file")
    parser.add_argument('-e', '--exon_gtf', required=True, help="Path to the exons only GTF genome file")
    parser.add_argument('-r', '--ref', required=True, help="Path to the reference (Human+MCPV) genome file")
    parser.add_argument('-o', '--output', required=True, help="Output directory for results")
    parser.add_argument('-s', '--sample_list', nargs='+', required=True, help="List of sample names")
    return parser.parse_args()

def run_sh_script(script_path):
    """Function to run the shell script."""
    os.system(f"bash {script_path}")

def main():
    args = parse_args()
    # Ensure output directory exists
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # Define paths and input files
    workdir = args.workdir
    gtf_file = args.gtf
    sample_list = args.sample_list
    ref_hum_mcpv = args.ref
    exon_gtf = args.exon_gtf
    
    MCV_dic = {'VP2': [[465, 1190]],
               'VP1': [[1156, 2427]],
               'large_T': [[2503, 4722], [5154, 5387]],
               'small_T': [[4827, 5387]]}

    # Prepare the MCV DataFrame
    mcvList = []
    for gene in MCV_dic:
        for each in MCV_dic[gene]:
            mcvList.append([gene, each[0], each[1]])

    mcvDf = pd.DataFrame(mcvList, columns=['gene', 'start', 'end'])
    direction = ['+', '+', '-', '-', '-']
    mcvDf['direction'] = direction

    # Adjust start and end coordinates based on direction
    newStart, newEnd = [], []
    for i, each in enumerate(direction):
        if each == '-':
            newStart.append(mcvList[i][2])
            newEnd.append(mcvList[i][1])
        else:
            newStart.append(mcvList[i][1])
            newEnd.append(mcvList[i][2])

    mcvDf['newStart'] = newStart
    mcvDf['newEnd'] = newEnd
    mcvDf = mcvDf.drop(['start', 'end'], axis=1)
    mcvDf.to_csv(os.path.join(output_dir, 'mcvDf.csv'), index = False)

    # Load the site dictionary from pickle files
    siteDic = {}
    for sample in sample_list:
        path = os.path.join(workdir, f'searcHPV_results/Sample_{sample}/call_fusion_virus')
        res = os.path.join(path, 'filteredSelectedContig.pickle')
        with open(res, 'rb') as res_file:
            resDic = pickle.load(res_file)
            if resDic != {}:
                siteDic[sample] = resDic

    # Process gene data from GTF file
    geneDic = {}
    with open(gtf_file) as genome:
        rows = genome.read().rstrip().split('\n')
        for eachRow in rows:
            hasGeneName = False
            for each in eachRow.split('\t')[8].split(';'):
                if each.split('"')[0] == ' gene_name ':
                    geneName = each.split('"')[1]
                    hasGeneName = True
                if each.split('"')[0] == ' transcript_id ':
                    trID = each.split('"')[1]
            if not hasGeneName:
                geneName = trID
            chro = eachRow.split('\t')[0]
            start = eachRow.split('\t')[3]
            end = eachRow.split('\t')[4]
            strand = eachRow.split('\t')[6]
            if trID not in geneDic:
                geneDic[trID] = [geneName, chro, start, end, strand]

    # Process gene sites
    geneSiteList = []
    for sample in siteDic:
        for ins in siteDic[sample]:
            chro = ins.split('.')[0]
            site = ins.split('.')[1]
            geneNameDic = {}
            flag = False
            for gene in geneDic:
                if geneDic[gene][1] == chro and int(geneDic[gene][2]) <= int(site) and int(geneDic[gene][3]) > int(site):
                    start = geneDic[gene][2]
                    end = geneDic[gene][3]
                    dirct = geneDic[gene][4]
                    geneNameDic[geneDic[gene][0]] = [start, end, dirct]
                    flag = True
            for gene in geneNameDic:
                geneSiteList.append([sample, gene, chro, site] + geneNameDic[gene])

    # Count occurrences of genes
    countDic = dict(pd.Series([gene for sample in geneSiteList for gene in sample[1]]).value_counts())
    
    genePlotDfList = []
    for each in geneSiteList:
        genePlotDfList.append(each)
    
    # Create DataFrame for plotting and output
    genePlotDf = pd.DataFrame(genePlotDfList, columns=['sample', 'gene', 'chr', 'ins', 'start', 'end', 'direction'])
    genePlotDf.to_csv(os.path.join(output_dir, 'geneTrack_MCV.csv'), index = False)

    # Write shell scripts
    with open(os.path.join(output_dir, 'extractContig.sh'), 'w') as scriptfile:
        for each in genePlotDf.itertuples():
            chro = each.chr
            ins = each.ins
            sample = each.sample
            gene = each.gene
            contigfile = os.path.join(workdir, f'searcHPV_results/Sample_{sample}/call_fusion_virus/ContigsSequence.fa')
            resfile = os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.fa')
            scriptfile.write(f'grep -A2 {chro}.{ins} {contigfile} >> {resfile};\n')
            
    with open(os.path.join(output_dir, 'alignContig.sh'), 'w') as scriptfile:
        #scriptfile.write('''module load Bioinformatics;
#module load bwa;
#module load samtools;
#''')
        for each in genePlotDf.itertuples():            
            gene = each.gene
            sample = each.sample
            resfile = os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.fa')
            ref = os.path.join(workdir, ref_hum_mcpv)
            scriptfile.write(f'''bwa mem -R '@RG\\tID:mcv\\tSM:mcv\\tLB:mcv\\tPL:ILLUMINA' -M -t 8 {ref} {resfile} > {os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.sam')};
    samtools view -bhS {os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.sam')} > {os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.bam')};
    samtools sort -@ 8 -o {os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.sort.bam')} {os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.bam')};
    samtools index {os.path.join(workdir, f'searcHPV_results/Sample_{sample}/{gene}.sort.bam')};''')

    
    # Run shell scripts
    run_sh_script(os.path.join(output_dir, 'extractContig.sh'))
    run_sh_script(os.path.join(output_dir, 'alignContig.sh'))
 
    exonDic = {}

    with open(exon_gtf) as inputFile:
        inputFile.readline()  # Skip header
        rows = inputFile.read().strip().split('\n')  # Avoid empty lines at the end
        for row in rows:
            if row.strip():  # Skip empty lines
                fields = row.split('\t')
                if len(fields) < 9:  # Skip incomplete rows
                    continue

                attributes = fields[8].split(';')
                hasExonNum = False
                hasGeneName = False
                exonNum, geneName, trID, exonID = None, None, None, None

                for each in attributes:
                    key_value = each.strip().split('"')
                    if len(key_value) < 2:
                        continue  # Skip malformed entries

                    if key_value[0].strip() == 'exon_number':
                        exonNum = key_value[1]
                        hasExonNum = True
                    elif key_value[0].strip() == 'gene_name':
                        geneName = key_value[1]
                        hasGeneName = True
                    elif key_value[0].strip() == 'transcript_id':
                        trID = key_value[1]
                    elif key_value[0].strip() == 'exon_id':
                        exonID = key_value[1]

                if not hasGeneName:
                    geneName = trID

                if hasExonNum:
                    exonStart = fields[3]
                    exonEnd = fields[4]
                    for each in genePlotDfList:
                        if each[1] == geneName:
                            exonDic[exonID] = [geneName,exonStart,exonEnd,exonNum,trID]

    annotation = []
    inserts_annotation = []
    for eachline in genePlotDfList:
        sample = eachline[0]
        gene = eachline[1]
        # if sample == 'Sample_2655-CB-89' and gene == 'KLF13':
        dirct = eachline[6]
        ins = eachline[3]
        path = workdir
        bam = f'{path}/searcHPV_results/Sample_{sample}/{gene}.sort.bam'
        chro = eachline[2]
        genestart = eachline[4]
        geneend = eachline[5]
        genedir = eachline[6]
        samfile = pysam.AlignmentFile(bam,'rb')
        inserts = []

        for read in samfile.fetch(chro,int(genestart),int(geneend)):
            start  = read.pos
            length = read.query_alignment_length
            fullLength = read.infer_read_length()
            name = read.qname
            contig = name.split('.')[2]
            contigins = name.split('.')[1]
            #print(sample,chro,ins)
            for eachContig in siteDic[sample][chro + '.' + contigins]:
                if eachContig[0] == contig:
                    if eachContig[-1] != 'lowConfidence':
                        mcvSite = eachContig[-1]
                    else:
                        mcvSite = eachContig[-2]

            mcvGeneBreak = ''
            for eachMCVgene in MCV_dic:
                for eachSlot in MCV_dic[eachMCVgene]:
                    if int(mcvSite) >= int(eachSlot[0]) and int(mcvSite) <= int(eachSlot[1]):
                        mcvGeneBreak = mcvGeneBreak + eachMCVgene + ';'
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            #add flag
            flag = 'intron'
            exonNumList = ''
            for exon in exonDic:
                exonStart = int(exonDic[exon][1])
                exonEnd = int(exonDic[exon][2])
                exonNum = int(exonDic[exon][3])

                if start <=  int(exonEnd) and start + length >= exonStart:
                    flag = 'exon'
                    exonNumList = exonNumList + exon + ';'




            #resolve contig hum components
            cigar = read.cigar
            infoList = [sample,gene,name,contig,'hum',flag,genedir,fullLength,mcvSite,'',mcvGeneBreak,contigins,genestart,geneend,length,exonNumList,strand]


            if cigar[0][0] == 0:
                #matched part, if strand of contig is different from gene direction, reverse compliment the cigar, switch the start and end
                if infoList[6] == '+':
                    pointStart = 0
                    pointEnd = cigar[0][1]
                else:
                    pointStart = fullLength
                    pointEnd = fullLength - cigar[0][1]

            else:
                pointStart = cigar[1][1]
                if cigar[-1][0] == 4 or cigar[-1][0] == 5:
                    pointEnd = fullLength - cigar[-1][1]
                else:
                    pointEnd = fullLength
                if infoList[6] == '-':
                    pointStart = fullLength - cigar[0][1]
                    pointEnd = fullLength - pointEnd
            inserts.append([pointStart,pointEnd] + infoList)


            if read.has_tag('SA'):
                saTag = read.get_tag('SA')
            else:
                print(read)
            #print(saTag)
            aligns = saTag.split(';')
            for align in aligns:
                if align != '':
                    supchro = align.split(',')[0]
                    if supchro == 'NC_010277.2':
                        suppos = align.split(',')[1]
                        supstrand = align.split(',')[2]
                        supcigar = align.split(',')[3]
                        slots = re.split(r'[A-Z]',supcigar)[:-1]
                        slotCigars = re.split(r'\d+',supcigar)[1:]
                        #print(align,slots,slotCigars)
                        i = 0
                        #add flag
                        if slotCigars[0] == 'M':
                            match = int(slots[0])
                            supEnd = mcvSite - match
                        else:
                            match = int(slots[1])
                            supEnd = mcvSite + match
                        newFlag = []
                        each1 = None
                        hasEach1 = False
                        for each in MCV_dic:
                            for eachone in MCV_dic[each]:
                                if int(suppos) >= int(eachone[0]) and int(suppos) <= int(eachone[1]):
                                    each1 = each
                                    newFlag.append(each1)
                                    hasEach1 = True
                                if int(suppos) + match >= int(eachone[0]) and int(suppos) + match <= int(eachone[1]):
                                    each2 = each
                                    if each2 != each1 and hasEach1:
                                        newFlag.append(each2)
                                    elif not hasEach1:
                                        newFlag.append(each2)
                        if not hasEach1 and each1 is None:
                            each1 = ""
                        if newFlag == []:
                            newFlag = ''
                        else:
                            newFlag = ';'.join(newFlag)
                        #print(newFlag)
                        infoList = [sample,gene,name,contig,'mcpv',newFlag,genedir,fullLength,mcvSite,supEnd,mcvGeneBreak,contigins,genestart,geneend,length,'',supstrand]
                        if slotCigars[0] == 'M':
                            if supstrand == strand: #read converted to consistent cigar to original human alignment
                                if infoList[6] == '+': #read on hum refence was not converted
                                    pointStart = 0
                                    pointEnd = int(slots[0])
                                else: #read on hum refence was not converted
                                    pointStart = fullLength
                                    pointEnd = fullLength - int(slots[0])
                            else: #read have inverted cigar to original human alignment
                                if infoList[6] == '+': #read on hum refence was not converted, converted to human alignments
                                    pointStart = fullLength
                                    pointEnd = fullLength - int(slots[0])
                                else:#read on hum refence was not converted, don't converted
                                    pointStart = 0
                                    pointEnd = int(slots[0])

                            inserts.append([pointStart,pointEnd] + infoList)
                            #print(supcigar,supstrand,strand)
                            #print([pointStart,pointEnd] + infoList)
                        else:

                            if supstrand == strand: #read converted to consistent cigar to original human alignment
                                pointStart = int(slots[0])
                                if slotCigars[-1] == 'S' or slotCigars[-1] == 'H':
                                    pointEnd = fullLength - int(slots[-1])
                                else:
                                    pointEnd = fullLength
                                if infoList[6] == '-': #read on hum refence was converted
                                    pointStart = fullLength - pointStart
                                    pointEnd = fullLength - pointEnd

                            else:

                                if slotCigars[-1] == 'S' or slotCigars[-1] == 'H':
                                    pointStart = int(slots[-1])
                                else:
                                    pointStart = 0
                                pointEnd = fullLength - int(slots[0])
                                #base on direction change order of start and end
                                #print(infoList[6])
                                #print(pointStart)
                                if infoList[6] == '-': #read on hum refence was not converted, converted to human alignments
                                    pointStart = fullLength - pointStart
                                    pointEnd = fullLength - pointEnd
                                inserts.append([pointEnd,pointStart] + infoList)



            annotation.append([sample,gene,sample,contig,chro,start,length,strand,contigins])
            #print(inserts)
        #print(suppos,int(suppos) + match,inserts)
        inserts_annotation += inserts
        #pd.DataFrame(inserts,columns = ['name','contig','start','end','type','feature','geneDirection','fullLength']).to_csv(outputFile2)


    # Save final output files
    annotation_df = pd.DataFrame(annotation, columns=['sampleID', 'gene', 'sampleName', 'contig', 'chr', 'start', 'length', 'strand', 'ins'])
    annotation_df.drop_duplicates(inplace=True)
    annotation_df.to_csv(os.path.join(output_dir, 'annotation.csv'), index=False)

    inserts_df = pd.DataFrame(inserts_annotation, columns=['start', 'end', 'sample', 'gene', 'name', 'contig', 'type', 'feature', 'geneDirection', 'fullLength', 'mcpvSite', 'mcpvEnd', 'mcpvGene', 'ins', 'geneStart', 'geneEnd', 'length', 'exonNum', 'strand'])
    inserts_df.drop_duplicates(inplace=True)
    inserts_df.to_csv(os.path.join(output_dir, 'inserts.csv'), index=False)

if __name__ == '__main__':
    main()

