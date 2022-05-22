#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess

def vir_map(ref_genome, fastq_r1, fastq_r2, threads, outdir, log):
    """
    Essa função roda as análises de mapeamento contra os genomas virais

    Argumentos:
        ref_genome = Genoma de referência
        fastq_r1 = Arquivo fastq R1 filtrado e sem os reads do hospedeiro
        fastq_r2 = Arquivo fastq R2 filtrado e sem os reads do hospedeiro
        threads = Número de threads
        outdir = Diretório de output
        log = Arquivo que log
    """
    print(f"bwa index {outdir}/ref_genomes/{ref_genome} &&\nbwa mem -f {threads} {outdir}/ref_genomes/{ref_genome} {outdir}/{fastq_r1} {outdir}/{fastq_r2} -o {outdir}/{ref_genome}.bam &&\nsamtools sort -o {outdir}/{ref_genome}.sorted.bam {outdir}/{ref_genome}.bam  &&\nsamtools index {outdir}/{ref_genome}.sorted.bam\n")

    subprocess.call([f"bwa index {outdir}/ref_genomes/{ref_genome} && \
                       bwa mem -f {threads} {outdir}/ref_genomes/{ref_genome} {outdir}/{fastq_r1} {outdir}/{fastq_r2} -o {outdir}/{ref_genome}.bam && \
                       samtools sort -o {outdir}/{ref_genome}.sorted.bam {outdir}/{ref_genome}.bam  && \
                       samtools index {outdir}/{ref_genome}.sorted.bam"], shell=True, stdin = log, stdout = log, stderr= log)

def get_consensus(ref_genome, count, sorted_bam, prefix, outdir, log):
    """
    Essa função gera os genomas consensos usando a ferramenta iVar

    Argumentos:
        ref_genome = Genoma de referência
        count = valor contador
        sorted_bam = Arquivo sorted.bam
        prefix = Prefixo
        outdir = Diretório de output
        log = Arquivo que log
    """
    print(f"samtools mpileup -aa -d 50000 -aa --reference {outdir}/ref_genomes/{ref_genome} -B {sorted_bam} | ivar consensus -p {outdir}/{prefix}-{ref_genome} -q 20 -t 0 -m 1 -n N\n")

    subprocess.call([f"samtools mpileup -aa -d 50000 -aa --reference {outdir}/ref_genomes/{ref_genome} -B {sorted_bam} | ivar consensus -p {outdir}/{prefix}-{ref_genome} -q 20 -t 0 -m 1 -n N"],
    shell=True, stdin = log, stdout = log, stderr= log)
    return(print("|"+"-"*8+f"Etapa 6.{count}: Montagem por Referência com {ref_genome}- Finalizada"+"-"*8+"|\n"))

class RunBWA():
    def __init__(self,ref_genome, fastq_r1, fastq_r2, threads, outdir, log):
        self.ref_genome = ref_genome
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.threads = threads
        self.outdir = outdir
        self.log = log

    def run_vir_map(self):
        vir_map(self.ref_genome, self.fastq_r1, self.fastq_r2, self.threads, self.outdir, self.log)

class GetConsensus():
    def __init__(self, ref_genome, count, sorted_bam, prefix, outdir, log):
        self.ref_genome = ref_genome
        self.count = count
        self.sorted_bam = sorted_bam
        self.prefix = prefix
        self.outdir = outdir
        self.log = log
    
    def run_get_consensus(self):
        get_consensus(self.ref_genome, self.count, self.sorted_bam, self.prefix, self.outdir, self.log)