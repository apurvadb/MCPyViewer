#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import os
import pickle
import pandas as pd
from collections import Counter
import pysam
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Process and analyze MCPyV integration and generate plots.")
    parser.add_argument('-w', '--workdir', required=True, help='Path to the working directory')
    parser.add_argument('-i', '--input_txt', required=True, help='Path to the sample txt file')
    parser.add_argument('-o', '--output_dir', required=True, help="Directory to store the output files.")
    return parser.parse_args()

def read_data(path):
    selected_new_contig = {}
    for sample in os.listdir(path):
        if sample.startswith('Sample_'):
            res_dic_path = os.path.join(path, sample, 'call_fusion_virus/filteredSelectedContig.pickle')
            if os.path.exists(res_dic_path):
                with open(res_dic_path, 'rb') as res_dic_file:
                    selected_new_contig[sample] = pickle.load(res_dic_file)
    return selected_new_contig

def process_samples(mcv_info, all_dic):
    MCV_gene_list = []
    MCV_dic = {'VP2': [[465, 1190]], 'VP1': [[1156, 2427]], 
               'large_T': [[2503, 4722], [5154, 5387]], 'small_T': [[4827, 5387]]}
    
    for sample in mcv_info['Sample']:
        if str(sample) != 'nan':
            sample_key = 'Sample_' + str(sample)
            if sample_key in all_dic:
                for ins in all_dic[sample_key]:
                    for contig in all_dic[sample_key][ins]:
                        MCV_site = contig[-1] if contig[-1] != 'lowConfidence' else contig[-2]
                        for MCV_gene in MCV_dic:
                            for slot in MCV_dic[MCV_gene]:
                                if slot[0] <= MCV_site <= slot[1]:
                                    MCV_gene_list.append(MCV_gene)
    
    count_dic = dict(Counter(MCV_gene_list))
    count_df_mcv = pd.DataFrame.from_dict(count_dic, orient='index')
    return count_df_mcv, count_dic, MCV_dic

def calculate_expected(count_dic, MCV_dic):
    tot_site = sum(count_dic.values())
    mcv_len_dic = {}
    tot_len = 0
    for site in MCV_dic:
        mcv_len = sum([int(slot[1]) - int(slot[0]) for slot in MCV_dic[site]])
        mcv_len_dic[site] = mcv_len
        tot_len += mcv_len
    expect_mcv = {each: round(tot_site * (mcv_len_dic[each] / tot_len)) for each in mcv_len_dic}
    return expect_mcv

def write_integration_files(mcv_info, all_dic, MCV_dic, output_dir):
    mcpv_integration_bed = os.path.join(output_dir, "MCVIntegration.bed")
    genome_integration_bed = os.path.join(output_dir, "genomeIntegration_MCV.bed")
    with open(genome_integration_bed,'w') as outputGenome:
        with open(mcpv_integration_bed,'w') as outputmcpv:
            for sample in mcv_info['Sample']:
                MCV_sample_list = []
                if str(sample) != 'nan':
                    sample_key = 'Sample_' + str(sample)
                    if sample_key in all_dic:
                        for ins in all_dic[sample_key]:
                            chro = ins.split('.')[0]
                            site = ins.split('.')[1]
                            for contig in all_dic[sample_key][ins]:
                                if contig[-1] != 'lowConfidence':
                                    MCV_site = contig[-1]
                                else:
                                    MCV_site = contig[-2]
                            has_gene = False
                            for MCV_gene in MCV_dic:
                                k = 1
                                for slot in MCV_dic[MCV_gene]:
                                    if MCV_site >= slot[0] and MCV_site <= slot[1]:
                                        has_gene = True
                                        if len(MCV_dic[MCV_gene]) > 1:
                                            outputmcpv.write(MCV_gene + f'_{str(k)}' + '\t' + str(MCV_site*40000) + '\t' + str(MCV_site* 40000) + '\n')
                                        else:
                                            outputmcpv.write(MCV_gene + '\t' + str(MCV_site*40000) + '\t' + str(MCV_site* 40000) + '\n')
                                        outputGenome.write('chr' + chro + '\t' + site + '\t' + site + '\n')
                                    k += 1
                            if not has_gene:
                                if MCV_site >= 0 and MCV_site <= 465:
                                    outputmcpv.write('intergenic_1' + '\t' + str(MCV_site*40000) + '\t' + str(MCV_site* 40000) + '\n') 
                                elif MCV_site >= 4722 and MCV_site <= 4827:
                                    outputmcpv.write('intergenic_3' + '\t' + str(MCV_site*40000) + '\t' + str(MCV_site* 40000) + '\n') 
                                elif MCV_site >= 2427 and MCV_site <= 2503: 
                                    outputmcpv.write('intergenic_2' + '\t' + str(MCV_site*40000) + '\t' + str(MCV_site* 40000) + '\n') 
                                else:
                                    print(MCV_site)
                                outputGenome.write('chr' + chro + '\t' + site + '\t' + site + '\n')
    

def get_block2(bam):
    # Function to get alignment blocks from BAM files
    sam = pysam.AlignmentFile(bam, "rb")
    blockList = {}  # blockList[contig] = [start,end,cigarString]
    high_quality_contig = []
    
    for read in sam.fetch():
        start = read.pos
        contig = read.qname
        blocks = read.get_blocks()
        cigar = read.cigartuples
        length = read.infer_read_length()
        
        if read.mapping_quality >= 0:
            # check if cigar is not None
            if cigar is not None:
                # if start with clip
                high_quality_contig.append(contig)
                start = 0
                end = 0
                for eachSlot in cigar:
                    end += eachSlot[1]
                    if not read.is_reverse:
                        if contig not in blockList:
                            blockList[contig] = [[start, end, eachSlot[0]]]
                        else:
                            blockList[contig].append([start, end, eachSlot[0]])
                    else:
                        if contig not in blockList:
                            blockList[contig] = [[length - end + 1, length - start + 1, eachSlot[0]]]
                        else:
                            blockList[contig].append([length - end + 1, length - start + 1, eachSlot[0]])
                    start += eachSlot[1]

    newBlockList = {}
    for contig in blockList:
        for each in blockList[contig]:
            if each[2] == 0:
                if contig not in newBlockList:
                    newBlockList[contig] = [[each[0], each[1]]]
                else:
                    newBlockList[contig].append([each[0], each[1]])
                
    return newBlockList, high_quality_contig

def count_overlap(selected_contig, genome_block, MCV_block):
    # Similar logic as in the previous part
    score_dic = {}
    for each in selected_contig:
        score = 0
        overlap = False
        clean_end = False
        for slot1 in genome_block.get(each, []):
            start1 = slot1[0]
            end1 = slot1[1]
            for slot2 in MCV_block.get(each, []):
                start2 = slot2[0]
                end2 = slot2[1]
                #overlap
                oldScore = score
                if start2 <= end1 and start1 <= end2:
                    overlap = True
                    if end1-start2 > end2-start1:
                        newScore = end2-start1 + 1
                    else:
                        newScore = end1-start2 + 1
                    if newScore > oldScore:
                        score = newScore
        #clean end
        if overlap == False:
            score = 0
            for slot1 in genome_block[each]:
                start1 = slot1[0]
                end1 = slot1[1]
                for slot2 in MCV_block[each]:
                    start2 = slot2[0]
                    end2 = slot2[1]
                    #clean end
                    if start2 == end1 + 1 or start1 == end2 + 1:
                        clean_end = True
                        score += 0
            #gap
            if overlap == False and clean_end == False:
                score = 0
                for slot1 in genome_block[each]:
                    start1 = slot1[0]
                    end1 = slot1[1]
                    for slot2 in MCV_block[each]:
                        start2 = slot2[0]
                        end2 = slot2[1]
                        #gap
                        oldScore = score
                        if end1 < start2 or end2 < start1:
                            if end1 < start2:
                                newScore2 = -(start2 - end1 - 1)
                            elif end2 < start1:
                                newScore2 = -(start1 - end2 - 1)
                            if newScore2 > oldScore and oldScore != 0:
                                score = newScore2
                                #print(newScore2)
                            elif oldScore == 0:
                                score = newScore2
        #filter score > 100 or <-100
#         if score < 100 and score > -100:
        score_dic[each] = score
    return(score_dic)

def generate_link_data(mcv_info, workdir):
    score_list = []
    path = os.path.join(workdir, 'searcHPV_results')

    for sample in mcv_info['Sample']:
        if str(sample) != 'nan':
            sample_name = 'Sample_' + str(sample)
            res_path = os.path.join(path, sample_name, 'call_fusion_virus')

            if os.path.isdir(res_path):
                final_res = os.path.join(res_path, 'filteredSelectedContig.pickle')
                with open(final_res, 'rb') as input_file:
                    res_dic = pickle.load(input_file)

                for each in res_dic:
                    if os.path.isdir(os.path.join(res_path, each)):
                        align_genome = os.path.join(res_path, each, f"{each}.contigToGenome.sort.bam")
                        align_mcv = os.path.join(res_path, each, f"{each}.contigToVirus.sort.bam")

                        genome_block, high_quality_contig = get_block2(align_genome)
                        MCV_block, high_quality_contig2 = get_block2(align_mcv)
                        
                        # Selected Contig
                        selected_contig = [eachSlot[0] for eachSlot in res_dic[each]]
        
                        # Filter selected contigs
                        selected_contig = [each for each in selected_contig if each in high_quality_contig]
                        selected_contig = [each for each in selected_contig if each in high_quality_contig2]
                
                        score_dic = count_overlap(selected_contig, genome_block, MCV_block)
                        for contig in score_dic:
                            score_list.append(score_dic[contig])
                    
                        #print(f"Sample {sample}, Contig {each}, Scores: {score_dic}")

    return score_list

def generate_bed_files(mcv_info, workdir, output_dir):
    found_bam = None
    
    # Iterate over the samples to find an available BAM file.
    for sample in mcv_info['Sample']:
        if str(sample) != 'nan':
            sample_path = os.path.join(workdir, f'searcHPV_results/Sample_{sample}/alignment/alignment.RG.indelre.mkdup.sort.bam')
            if os.path.isfile(sample_path):
                found_bam = sample_path
                break

    if not found_bam:
        raise FileNotFoundError("No BAM file found for any samples in the list.")

    test = pysam.AlignmentFile(found_bam, "rb")
    
    #sample = os.path.join(workdir, 'searcHPV_results/Sample_152-CB-58/alignment/alignment.RG.indelre.mkdup.sort.bam')
    #test = pysam.AlignmentFile(sample, "rb")
    chr_len = {}
    for i in range(0,22):
        chr_len[int(test.references[i])] = test.lengths[i]
        
    ordered_chr_len = {}
    for key in sorted(chr_len):
        ordered_chr_len[key] = chr_len[key]
    ordered_chr_len[test.references[23]] = test.lengths[23]
    ordered_chr_len[test.references[24]] = test.lengths[24]
    
    
    MCV_dic = {
        'VP2': [[465, 1190]], 'VP1': [[1156, 2427]],
        'large_T': [[2503, 4722], [5154, 5387]], 'small_T': [[4827, 5387]]
    }

    bed_file = os.path.join(output_dir, 'MCV_track.bed')
    with open(bed_file,'w') as output:
        for gene in MCV_dic:
            k = 1
            for slot in MCV_dic[gene]:
                start = str((slot[0]-1)*40000)
                end = str((slot[1]-1)*40000)
            
                if len(MCV_dic[gene]) > 1:
                    output.write(f'{gene}_{str(k)}\t{start}\t{end}\n')
                else:
                    output.write(f'{gene}\t{start}\t{end}\n')
                k += 1
        inter1_end = str(464*40000)
        inter2_start = str(2426*40000)
        inter2_end = str(2502*40000)
        inter3_start = str(4721*40000)
        inter3_end = str(4826*40000)
    
        output.write(f'intergenic_1\t0\t{inter1_end}\n')
        output.write(f'intergenic_2\t{inter2_start}\t{inter2_end}\n')
        output.write(f'intergenic_3\t{inter3_start}\t{inter3_end}\n')

        for chro in ordered_chr_len:
            end = str(ordered_chr_len[chro])
            output.write(f'chr{chro}\t0\t{end}\n')
    


    label_file = os.path.join(output_dir, 'label.bed')
    with open(label_file, 'w') as output:
        for gene in MCV_dic:
            k = 1
            for slot in MCV_dic[gene]:
                start = str((slot[0]-1)*40000/2)
                end = str((slot[1]-1)*40000/2)
            
                if len(MCV_dic[gene]) > 1:
                    output.write(f'{gene}_{str(k)}\t{start}\t{end}\t{gene}\n')
                else:
                    output.write(f'{gene}\t{start}\t{end}\t{gene}\n')
                k += 1
        inter1_end = str(464*40000/2)
        inter2_start = str(2426*40000/2)
        inter2_end = str(2502*40000/2)
        inter3_start = str(4721*40000/2)
        inter3_end = str(4826*40000/2)
        output.write(f'intergenic_1\t0\t{inter1_end}\tintergenic\n')
        output.write(f'intergenic_2\t{inter2_start}\t{inter2_end}\tintergenic\n')
        output.write(f'intergenic_3\t{inter3_start}\t{inter3_end}\tintergenic\n')

        for chro in ordered_chr_len:
            end = str(ordered_chr_len[chro]/2)
            output.write(f'chr{chro}\t0\t{end}\tchr{chro}\n')
    
def main():
    args = parse_args()

    # Load data
    mcv_info = pd.read_csv(args.input_txt, sep = "\t")
    all_dic = read_data(f'{args.workdir}/searcHPV_results/')

    # Process samples
    count_df_mcv, count_dic, MCV_dic = process_samples(mcv_info, all_dic)
    expect_mcv = calculate_expected(count_dic, MCV_dic)

    # Create DataFrames
    expect_mcv_df = pd.DataFrame.from_dict(expect_mcv, orient='index')

    # Debug print statements
    #print("Count DF MCV Shape:", count_df_mcv.shape)
    #print("Count DF MCV Head:\n", count_df_mcv.head())
    #print("Expect MCV DF Shape:", expect_mcv_df.shape)
    #print("Expect MCV DF Head:\n", expect_mcv_df.head())

    # Merge DataFrames
    combined_df = pd.merge(count_df_mcv, expect_mcv_df, left_index=True, right_index=True, how='outer')

    # Debug print to understand merged DataFrame structure
    #print("Combined DF Shape:", combined_df.shape)
    #print("Combined DF Head:\n", combined_df.head())

    # Check the expected number of columns before setting new column names
    #if combined_df.shape[1] != 2:
    #    raise ValueError(f"Unexpected number of columns in combined_df: {combined_df.shape[1]}")

    # Set column names
    combined_df.columns = ['observed', 'expected']

    # Save to CSV
    combined_df.to_csv(os.path.join(args.output_dir, 'MCV_gene.csv'))

    # Generate other output files
    write_integration_files(mcv_info, all_dic, MCV_dic, args.output_dir)
    score_list = generate_link_data(mcv_info, args.workdir)
    #print("Scores:", score_list)
    #print("First 10 items of score_list:")
    #for item in score_list[:10]:
    #    print(item)
    
    score_df = pd.DataFrame(score_list, columns=['Scores'])
    score_df.to_csv(os.path.join(args.output_dir, 'score.csv'), index=False)

    # Write score_list to CSV
    #score_csv_path = os.path.join(args.output_dir, 'score_list.csv')
    #with open(score_csv_path, 'w', newline='') as csvfile:
    #    writer = csv.writer(csvfile)
    #    writer.writerow(['Sample', 'Contig', 'Selected Contig', 'Score'])  # CSV header
    #    for score_row in score_list:
    #        writer.writerow(score_row)

    generate_bed_files(mcv_info, args.workdir, args.output_dir)

if __name__ == "__main__":
    main()
