#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess

def human_map(ref_genome,fastq_r1,fastq_r2,threads,outdir,prefix,log):
    """
    Essa função roda a ferramenta bowtie2 e gera arquivos fastq com reads nao mapeados

    Argumentos:
        ref_genome = Genoma de referência
        fastq_r1 = Arquivo fastq R1 filtrado
        fastq_r2 = Arquivo fastq R2 filtrado
        threads = Número de threads
        outdir = Diretório de output
        prefix = Prefixo
        log = Arquivo que log
    """
    print(f"bowtie2 --un-conc-gz {outdir} -p {threads} -x {ref_genome} -1 {outdir}/{fastq_r1} -2 {outdir}/{fastq_r2} > {outdir}/{prefix}.human.bam &&\nmv {outdir}/un-conc-mate.1 {outdir}/{prefix}.un_R1.fq.gz &&\nmv {outdir}/un-conc-mate.2 {outdir}/{prefix}.un_R2.fq.gz\n")

    subprocess.call([f"bowtie2 --un-conc-gz {outdir} -p {threads} -x {ref_genome} -1 {outdir}/{fastq_r1} -2 {outdir}/{fastq_r2} > {outdir}/{prefix}.human.bam && \
                       mv {outdir}/un-conc-mate.1 {outdir}/{prefix}.un_R1.fq.gz && \
                       mv {outdir}/un-conc-mate.2 {outdir}/{prefix}.un_R2.fq.gz"], shell=True, stdin = log, stdout = log, stderr= log)
    return(print("|"+"-"*14+"Etapa 2: Remoção de Reads do Hospedeiro - Finalizada"+"-"*15+"|\n"))


class RunBOWTIE2():
    def __init__(self,ref_genome, fastq_r1,fastq_r2,threads, outdir,prefix,log):
        self.ref_genome = ref_genome
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.threads = threads
        self.outdir = outdir
        self.prefix = prefix
        self.log = log
    
    def run_bowtie2(self):
        human_map(self.ref_genome, self.fastq_r1, self.fastq_r2, self.threads, self.outdir, self.prefix,self.log)