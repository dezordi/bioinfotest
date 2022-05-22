#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess
import pandas as pd

def blastn(database_file, scaffolds, threads, outdir, prefix,log):
    """
    Essa função realiza a análise de BLAST

    Argumentos:
        database_file = arquivo com genomas de referencia usados na analise de blast
        scaffolds = arquivo de scaffols
        threads = número de threads
        outdir = Diretório de output
        prefix = Prefixo
        log = Arquivo que log
    """
    print(f"makeblastdb -in {database_file} -dbtype nucl &&\nblastn -task blastn -query {scaffolds} -db {database_file} -out {outdir}/{prefix}.blast.tsv -outfmt '6 qseqid sseqid evalue bitscore length pident qcovhsp' -evalue 1e-6 -max_target_seqs 1 -num_threads {threads}\n")

    subprocess.call([f"makeblastdb -in {database_file} -dbtype nucl && \
                       blastn -task blastn -query {scaffolds} -db {database_file} -out {outdir}/{prefix}.blast.tsv -outfmt '6 qseqid sseqid evalue bitscore length pident qcovhsp' -evalue 1e-6 -max_target_seqs 1 -num_threads {threads}"], shell=True, stdin = log, stdout = log, stderr= log)

def fillter_blast_output(blast_output, database_table , len_cutoff, outdir, prefix):
    """
    Essa função realiza o filtro dos resultados de blast

    Argumentos:
        blast_output = arquivo com os resultados do blast
        database_table = arquivo com as informações dos genomas utilizados na analise de blast
        len_cutoff = tamanho minimo para os hits serem considerados
        outdir = Diretório de output
        prefix = Prefixo
        log = Arquivo que log
    """
    headers = ['query','Accession','evalue','bitscore','length','pident','qcovhsp']
    blast_dataframe = pd.read_csv(blast_output, sep='\t', names = headers)
    blast_dataframe = blast_dataframe.loc[(blast_dataframe['pident'] >= 80) & (blast_dataframe['qcovhsp'] >= 80) & (blast_dataframe['length'] >= int(len_cutoff))]
    vir_dataframe = pd.read_csv(database_table)
    blast_annot_dataframe = pd.merge(blast_dataframe, vir_dataframe, on="Accession", how='left') #isso sera retornado
    with open(f"{outdir}/{prefix}.blast.filtred.tsv",'w') as blast_filtred_writer:
        blast_annot_dataframe.to_csv(blast_filtred_writer, index=False, header=True)
    return(print("|"+"-"*21+"Etapa 4: Análise de BLAST - Finalizada"+"-"*22+"|\n"))

class RunBLASTN():
    def __init__(self,database_file, scaffolds, threads,outdir,prefix,log):
        self.database_file = database_file
        self.scaffolds = scaffolds
        self.threads = threads
        self.outdir = outdir
        self.prefix = prefix
        self.log = log
    
    def run_blastn(self):
        blastn(self.database_file, self.scaffolds, self.threads, self.outdir, self.prefix,self.log)

class FilterBLASTN():
    def __init__(self, blast_output, database_table , len_cutoff, outdir, prefix):
        self.blast_output = blast_output
        self.database_table = database_table
        self.len_cutoff = len_cutoff
        self.outdir = outdir
        self.prefix = prefix

    def run_blast_filter(self):
        fillter_blast_output(self.blast_output, self.database_table, self.len_cutoff, self.outdir, self.prefix)