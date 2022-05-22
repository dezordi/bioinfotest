#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess

def denovo_assembly(fastq_r1,fastq_r2,threads,mem_gb,outdir,prefix,log):
    """
    Essa função roda a ferramenta spades

    Argumentos:
        fastq_r1 = Arquivo fastq R1 filtrado e sem os reads do hospedeiro
        fastq_r2 = Arquivo fastq R2 filtrado e sem os reads do hospedeiro
        threads = Número de threads
        mem_gb = Valor em Gb de RAM
        outdir = Diretório de output
        prefix = Prefixo
        log = Arquivo que log
    """
    print(f"spades.py --meta -1 {fastq_r1} -2 {fastq_r2} -t {threads} -m {mem_gb} -o {outdir} &&\nmv {outdir}/scaffolds.fasta {outdir}/{prefix}.scaffolds.fa\n")

    subprocess.call([f"spades.py --meta -1 {fastq_r1} -2 {fastq_r2} -t {threads} -m {mem_gb} -o {outdir}/ && \
                       mv {outdir}/scaffolds.fasta {outdir}/{prefix}.scaffolds.fa"], shell=True, stdin = log, stdout = log, stderr= log)
    return(print("|"+"-"*22+"Etapa 3: Montagem denovo - Finalizada"+"-"*22+"|\n"))


class RunSPADES():
    def __init__(self,fastq_r1,fastq_r2,threads, mem_gb, outdir,prefix,log):
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.threads = threads
        self.outdir = outdir
        self.prefix = prefix
        self.mem_gb = mem_gb
        self.log = log
    
    def run_denovo_assembly(self):
        denovo_assembly(self.fastq_r1, self.fastq_r2, self.threads, self.mem_gb,self.outdir, self.prefix,self.log)