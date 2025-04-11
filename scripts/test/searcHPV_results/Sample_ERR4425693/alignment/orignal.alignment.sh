
    bwa mem -t 1 /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/hg_hpv.fa '<cat ERR4425693_MCC_reads_1.fq' '<cat ERR4425693_MCC_reads_2.fq' > /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/alignment/alignment.sam
    samtools view -@ 1 -bhS /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/alignment/alignment.sam > /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/alignment/alignment.bam
    samtools sort -@ 1 /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/alignment/alignment.bam -o /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/alignment/alignment.sort.bam
    rm /gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/searcHPV_results/Sample_ERR4425693/alignment/alignment.sam
    echo 'alignment done'
    