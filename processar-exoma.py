#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 18:38:29 2026

@author: vitor
"""

import subprocess
import os
from concurrent.futures import ThreadPoolExecutor

caminho = "/home/venus/mar/vtarghetta/sra-toolkit/phs001493/"
lista = sorted(os.listdir(caminho)) #[:2]

def processar_amostra(amostra):
    caminho_amostra = os.path.join(caminho, amostra)
    
    # Processa apenas se for diretório
    if not os.path.isdir(caminho_amostra):
        return

    print(f"-> Iniciando processamento para: {amostra}")

    # FIX: Garantir que os diretórios necessários existam
    subpastas = ["fastqc", "results/trimmed", "results/aligned", "results/qc", "results/variants", "tmp"]
    for pasta in subpastas:
        os.makedirs(os.path.join(caminho_amostra, pasta), exist_ok=True)

    # Definição dos Comandos
    fastqc = f"fastqc {caminho_amostra}/{amostra}_1.fastq.gz {caminho_amostra}/{amostra}_2.fastq.gz -o {caminho_amostra}/fastqc/"
    
    trim = f"trim_galore --paired --quality 20 --length 50 --fastqc --output_dir {caminho_amostra}/results/trimmed/ {caminho_amostra}/{amostra}_1.fastq.gz {caminho_amostra}/{amostra}_2.fastq.gz"
    
    bwa_mem = f'bwa mem -t 12 -M -R "@RG\\tID:{amostra}\\tSM:{amostra}\\tPL:ILLUMINA\\tLB:{amostra}_lib" /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta {caminho_amostra}/results/trimmed/{amostra}_1_val_1.fq.gz {caminho_amostra}/results/trimmed/{amostra}_2_val_2.fq.gz | samtools sort -@ 12 -o {caminho_amostra}/results/aligned/{amostra}.bam'
    
    mark_duplicates = f"gatk MarkDuplicates -I {caminho_amostra}/results/aligned/{amostra}.bam -O {caminho_amostra}/results/aligned/{amostra}_marked_duplicates.bam -M {caminho_amostra}/results/aligned/{amostra}_duplicate_metrics.txt --CREATE_INDEX true --TMP_DIR {caminho_amostra}/tmp/"

    # FIX: Corrigido --known-sities para --known-sites
    base_recalibrator = f"gatk BaseRecalibrator -I {caminho_amostra}/results/aligned/{amostra}_marked_duplicates.bam -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta --known-sites /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /home/venus/mar/vtarghetta/wxs_analysis/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O {caminho_amostra}/results/aligned/{amostra}_recal_data.table"
    
    applyBQSR = f"gatk ApplyBQSR -I {caminho_amostra}/results/aligned/{amostra}_marked_duplicates.bam -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta --bqsr-recal-file {caminho_amostra}/results/aligned/{amostra}_recal_data.table -O {caminho_amostra}/results/aligned/{amostra}_recalibrated.bam"
    
    collectHSmetrics = f"gatk CollectHsMetrics -I {caminho_amostra}/results/aligned/{amostra}_recalibrated.bam -O {caminho_amostra}/results/qc/{amostra}_hs_metrics.txt -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta -BAIT_INTERVALS /home/venus/mar/vtarghetta/bed/20260525-exome_targets.interval_list -TARGET_INTERVALS /home/venus/mar/vtarghetta/bed/20260525-exome_targets.interval_list"
    
    depthOfCoverage = f"gatk DepthOfCoverage -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta -I {caminho_amostra}/results/aligned/{amostra}_recalibrated.bam -O {caminho_amostra}/results/qc/{amostra}_coverage -L /home/venus/mar/vtarghetta/bed/20260525-comprehensiveBED.bed --summary-coverage-threshold 10 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 50 --summary-coverage-threshold 100"
    
    haplotype_caller = f"gatk HaplotypeCaller -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta -I {caminho_amostra}/results/aligned/{amostra}_recalibrated.bam -O {caminho_amostra}/results/variants/{amostra}.g.vcf.gz -ERC GVCF -L /home/venus/mar/vtarghetta/bed/20260525-comprehensiveBED.bed --dbsnp /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    
    genotypeGVCFs = f"gatk GenotypeGVCFs -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta -V {caminho_amostra}/results/variants/{amostra}.g.vcf.gz -O {caminho_amostra}/results/variants/{amostra}_raw_variants.vcf.gz -L /home/venus/mar/vtarghetta/bed/20260525-comprehensiveBED.bed"
    
    variantFiltration = f'gatk VariantFiltration -R /home/venus/mar/vtarghetta/wxs_analysis/reference/Homo_sapiens_assembly38.fasta -V {caminho_amostra}/results/variants/{amostra}_raw_variants.vcf.gz \
        --filter-expression "QD < 2.0" --filter-name "LowQD" \
        --filter-expression "FS > 60.0" --filter-name "HighFS" \
        --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
        --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
        --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum" \
        --filter-expression "DP < 20" --filter-name "LowDepth" \
        --filter-expression "DP > 500" --filter-name "HighDepth" \
        -O {caminho_amostra}/results/variants/{amostra}_filtered.vcf.gz'
        
    bcf_extract = f"bcftools view -f PASS -Oz -o {caminho_amostra}/results/variants/{amostra}_pass.vcf.gz {caminho_amostra}/results/variants/{amostra}_filtered.vcf.gz"
    
    index = f"bcftools index {caminho_amostra}/results/variants/{amostra}_pass.vcf.gz"

    # Pasta onde os logs da amostra serão salvos
    pasta_logs = os.path.join(caminho_amostra, "logs")
    os.makedirs(pasta_logs, exist_ok=True)

    # Adicionamos um terceiro elemento na tupla: o nome base do arquivo de log
    etapas = [
        ("FastQC", fastqc, "01_fastqc"),
        ("Trim Galore", trim, "02_trim_galore"),
        ("BWA MEM + Samtools", bwa_mem, "03_bwa_mem"),
        ("MarkDuplicates", mark_duplicates, "04_mark_duplicates"),
        ("BaseRecalibrator", base_recalibrator, "05_base_recalibrator"),
        ("ApplyBQSR", applyBQSR, "06_apply_bqsr"),
        ("CollectHsMetrics", collectHSmetrics, "07_collect_hs_metrics"),
        ("DepthOfCoverage", depthOfCoverage, "08_depth_of_coverage"),
        ("HaplotypeCaller", haplotype_caller, "09_haplotype_caller"),
        ("GenotypeGVCFs", genotypeGVCFs, "10_genotype_gvcfs"),
        ("VariantFiltration", variantFiltration, "11_variant_filtration"),
        ("BCFTools Extract", bcf_extract, "12_bcf_extract"),
        ("BCFTools Index", index, "13_bcf_index")
    ]

    # Execução das etapas com salvamento de log
    for nome, cmd, id_log in etapas:
        log_file = os.path.join(pasta_logs, f"{amostra}_{id_log}.log")
        print(f"<- Iniciando {nome} para: {amostra} (Log: {log_file})")
        
        # Abre o arquivo de log e redireciona stdout e stderr para ele
        with open(log_file, "w") as f_log:
            subprocess.run(
                cmd, 
                shell=True, 
                check=True, 
                stdout=f_log, 
                stderr=subprocess.STDOUT
            )
            
        print(f"<- Finalizado {nome} para: {amostra}")



# Executa em paralelo (Ajuste max_workers de acordo com a memória/CPUs disponíveis)
if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=6) as executor:
        executor.map(processar_amostra, lista)
