#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess

def fastp(fastq_r1,fastq_r2,threads, outdir,prefix,log):
    """
    Essa função roda a ferramenta fastp

    Argumentos:
        fastq_r1 = Arquivo fastq R1 filtrado e sem os reads do hospedeiro
        fastq_r2 = Arquivo fastq R2 filtrado e sem os reads do hospedeiro
        threads = Número de threads
        outdir = Diretório de output
        prefix = Prefixo
        log = Arquivo que log
    """
    print(f"fastp -i {fastq_r1} -I {fastq_r2} -w {threads} -o {outdir}/{prefix}.R1.fq.gz -O {outdir}/{prefix}.R2.fq.gz -h {outdir}/{prefix}.fastp.html -j {outdir}/{prefix}.fastp.json --cut_front --qualified_quality_phred 30\n")
    
    subprocess.call([f"fastp -i {fastq_r1} -I {fastq_r2} -w {threads} -o {outdir}/{prefix}.R1.fq.gz -O {outdir}/{prefix}.R2.fq.gz -h {outdir}/{prefix}.fastp.html -j {outdir}/{prefix}.fastp.json --cut_front --qualified_quality_phred 30"], shell=True, stdin = log, stdout = log, stderr= log)
    return(print("|"+"-"*19+"Etapa 1: Controle de Qualidade - Finalizada"+"-"*19+"|\n"))


class RunFASTP():
    def __init__(self, fastq_r1,fastq_r2,threads, outdir,prefix,log):
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.threads = threads
        self.outdir = outdir
        self.prefix = prefix
        self.log = log
    
    def run_fastp(self):
        fastp(self.fastq_r1, self.fastq_r2, self.threads, self.outdir, self.prefix,self.log)