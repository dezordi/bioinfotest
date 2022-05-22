#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from Bio import SeqIO

def get_refs(blast_filtred_output,database_file, outdir):
    """
    Essa função recupera os genomas de referência reportados no blast filtrado

    Argumentos:
        blast_filtred_output = Arquivo tsv de blast filtrado
        database_file = Arquivo com genomas de referencia usados na analise de blast
        outdir = Diretório de output
    """
    print("Rodando funções internas para recuperação dos genomas...\n")

    blast_filtred_dataframe = pd.read_csv(blast_filtred_output)
    os.mkdir(f"{outdir}/ref_genomes/")
    acc_list = blast_filtred_dataframe['Accession'].tolist()
    for accession in acc_list:
        for record in SeqIO.parse(database_file,'fasta'):
            if accession == record.id:
                with open(f"{outdir}/ref_genomes/{accession}.fa",'w') as output_accession:
                    SeqIO.write(record, output_accession,'fasta')
    return(print("|"+"-"*11+"Etapa 5: Recuperação de Genomas de Referência - Finalizada"+"-"*12+"|\n"))

class GetREFS():
    def __init__(self,blast_filtred_output,database_file, outdir):
        self.blast_filtred_output = blast_filtred_output
        self.database_file = database_file
        self.outdir = outdir

    def run_get_refs(self):
        get_refs(self.blast_filtred_output, self.database_file, self.outdir)