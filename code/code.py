#!/usr/bin/python3
# -*- coding: utf-8 -*-

##################################>INFORMACOES GERAIS<##################################

__author__ = {"Nomes": "Filipe Dezordi",
               "email": "zimmer.filipe@gmail.com",
               "git": "https://github.com/dezordi"}
__license__ = "GPL"
__version__ = "0.1dev"
__maintainer__ = "Filipe Dezordi"
__status__ = "Prototype"

##################################>IMPORTANDO MODULOS<##################################

import argparse, re, os, glob, sys
from pathlib import Path
from modulos.script_info import ScriptInfo
from modulos.quality_control import RunFASTP
from modulos.remove_human import RunBOWTIE2
from modulos.denovo_assembly import RunSPADES
from modulos.blast_vir import RunBLASTN, FilterBLASTN
from modulos.get_refs import GetREFS
from modulos.ref_map import RunBWA, GetConsensus

######################################>ARGUMENTOS<######################################

parser = argparse.ArgumentParser(description='Esse script realiza a montagem de genomas virais em análises de metagenômica')
parser.add_argument("-r", "--ref_host",
                    help="Arquivo Fasta do genoma do hospedeiro", required=True)
parser.add_argument("-vd", "--vir_database",
                    help="Aquivo Fasta do banco de dados de genomas virais",required=True)
parser.add_argument("-mt", "--vir_metadata",
                    help="Arquivo CSV do banco de dados de genomas virais.",required=True)
parser.add_argument("-f1", "--fastq_r1",
                    help="Arquivo FASTQ R1",required=True)
parser.add_argument("-f2", "--fastq_r2",
                    help="Arquivo FASTQ R2",required=True)
parser.add_argument("-p", "--threads",
                    help="Número de Threads, default = 1.", default=1)
parser.add_argument("-m", "--mem_gb",
                    help="Quantos GB de memória RAM para rodar a análise, default = 2.", default=2)
parser.add_argument("-c", "--cut_off",
                    help="Tamanho minimo do alinhamento para ser considerado nas analises de BLAST, default = 1000.", default=1000)
parser.add_argument("-od", "--outdir",
                    help="Caminho e nome do diretorio de output")
parser.add_argument("-pr", "--prefix",
                    help="Prefixo a ser usado no nome dos arquivos de output")

# Armazenando argumentos em variaveis
args = parser.parse_args()
ref_host = args.ref_host
database_file = args.vir_database
database_table = args.vir_metadata
fastq_r1 = args.fastq_r1
fastq_r2 = args.fastq_r2
threads = args.threads
mem_gb = args.mem_gb
len_cutoff = args.cut_off
out_dir = args.outdir
prefix = args.prefix

##################################>CHECANDO ARGUMENTOS<###############################

list_required_args = [ref_host,database_file,database_table,fastq_r1,fastq_r2]
def check_files(file):
    if Path(file).is_file() == False:
        sys.exit(f"\nO arquivo {file} não é valido. Interrompendo execução do código.\n")
for required_arg in list_required_args:
    check_files(required_arg)

# Se não foi setado um prefixo, usar o nome da amostra
if prefix == None:
    prefix = fastq_r1
    prefix = re.sub("\.\..*\/", "", prefix).rstrip('\n')
    prefix = re.sub('.*/', '', prefix)
    prefix = re.sub('.*/', '', prefix)
    prefix = re.sub('_.*', '', prefix)

# Checar se o diretorio de output foi passado, se não, gerar o nome de acordo com o prefixo, e criar o diretório caso ele não exista
if out_dir != None:
    if "/" in out_dir:
        out_dir = re.sub("/$", "", out_dir)
    else:
        out_dir = './'+out_dir
else:
    out_dir = re.sub("/.*", "/", prefix).rstrip('\n')
    out_dir = './'+out_dir

if os.path.isdir(out_dir) == False:
    os.mkdir(out_dir)

# Criando arquivo de log
log_file = open(f"{out_dir}/run.log.txt","w+")

# Iniciando a análise
if __name__ == '__main__':
    print_info = ScriptInfo(ref_host, database_file, database_table, fastq_r1, fastq_r2, threads, mem_gb, len_cutoff, out_dir, prefix)
    print_info.print_message()
    print("|"+"-"*20+"Etapa 1: Controle de Qualidade - Iniciada"+"-"*20+"|\n")
    run_fastp = RunFASTP(fastq_r1,fastq_r2, threads, out_dir, prefix, log_file)
    run_fastp.run_fastp()
    print("|"+"-"*15+"Etapa 2: Remoção de Reads do Hospedeiro - Iniciada"+"-"*16+"|\n")
    run_bowtie2 = RunBOWTIE2(ref_host,
                             f"{prefix}.R1.fq.gz",f"{prefix}.R2.fq.gz",
                             threads, out_dir, prefix, log_file)
    run_bowtie2.run_bowtie2()
    print("|"+"-"*23+"Etapa 3: Montagem denovo - Iniciada"+"-"*23+"|\n")
    run_spades = RunSPADES(f"{out_dir}/{prefix}.un_R1.fq.gz",
                           f"{out_dir}/{prefix}.un_R2.fq.gz",
                           threads, mem_gb, out_dir, prefix, log_file)
    run_spades.run_denovo_assembly()
    print("|"+"-"*22+"Etapa 4: Análise de BLAST - Iniciada"+"-"*23+"|\n")
    run_blastn = RunBLASTN(database_file,f"{out_dir}/{prefix}.scaffolds.fa",
                           threads, out_dir, prefix, log_file)
    run_blastn.run_blastn()
    filter_blastn = FilterBLASTN(f"{out_dir}/{prefix}.blast.tsv",
                                    database_table,
                                    len_cutoff,
                                    out_dir,
                                    prefix)
    filter_blastn.run_blast_filter()
    print("|"+"-"*12+"Etapa 5: Recuperação de Genomas de Referência - Iniciada"+"-"*13+"|\n")
    get_refs = GetREFS(f"{out_dir}/{prefix}.blast.filtred.tsv", database_file, out_dir)
    get_refs.run_get_refs()
    print("|"+"-"*19+"Etapa 6: Montagem por Referência - Iniciada"+"-"*19+"|\n")
    count = 1
    for ref_genome in glob.glob(f"{out_dir}/ref_genomes/*.fa"):
        ref_genome = re.sub(r".*/","",ref_genome).rstrip('\n')
        print("|"+"-"*9+f"Etapa 6.{count}: Montagem por Referência com {ref_genome}- Iniciada"+"-"*9+"|\n")
        un_map_vir = RunBWA(ref_genome,
                            f"{prefix}.un_R1.fq.gz",
                            f"{prefix}.un_R2.fq.gz",
                            threads, out_dir, log_file)
        un_map_vir.run_vir_map()
        get_consensus = GetConsensus(ref_genome, count,
                                     f"{out_dir}/{ref_genome}.sorted.bam",
                                     prefix, out_dir, log_file)
        get_consensus.run_get_consensus()
        count += 1
    print("|"+"-"*18+"Etapa 6: Montagem por Referência - Finalizada"+"-"*18+"|\n")