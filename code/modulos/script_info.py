#!/usr/bin/python3
# -*- coding: utf-8 -*-

def print_script_info(ref_human, database_file, database_table, fastq_r1, fastq_r2, threads, mem_gb, len_cutoff, out_dir, prefix):
    print("\n")
    print("|"+"-"*32+"INICIANDO ANÁLISE"+"-"*32+"|\n|")
    print(f"| Genoma de Referência para filtro:  {ref_human}")
    print(f"| Banco de Genomas Virais:           {database_file}")
    print(f"| Metadados Virais:                  {database_table}")
    print(f"| Arquivo FASTQ R1:                  {fastq_r1}")
    print(f"| Arquivo FASTQ R2:                  {fastq_r2}")
    print(f"| Número de Threads:                 {threads}")
    print(f"| Memória RAM:                       {mem_gb}Gb")
    print(f"| Tamanho minimo para scaffolds:     {len_cutoff}")
    print(f"| Diretório de Output:               {out_dir}")
    print(f"| Prefixo dos arquivos output:       {prefix}\n|")
    print("|"+"-"*81+"|")
    print("\n")

class ScriptInfo():
    def __init__(self, ref_human, database_file, database_table, fastq_r1, fastq_r2, threads, mem_gb, len_cutoff, out_dir, prefix):
        self.ref_human = ref_human
        self.database_file = database_file
        self.database_table = database_table
        self.fastq_r1 = fastq_r1
        self.fastq_r2 = fastq_r2
        self.threads = threads
        self.mem_gb = mem_gb
        self.len_cutoff = len_cutoff
        self.out_dir = out_dir
        self.prefix = prefix
    def print_message(self):
        print_script_info(self.ref_human, self.database_file, self.database_table, self.fastq_r1, self.fastq_r2, self.threads, self.mem_gb, self.len_cutoff,self.out_dir,self.prefix)