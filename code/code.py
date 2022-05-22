#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv

from numpy import append
from modulos.get_kmer import GetKmers
from modulos.get_nodes import GetNodes

input_file = argv[1]
k_length = argv[2]

#Primeira etapa, obter os kmers
kmer_set = []
with open(input_file, 'r') as input_file_reader:
    for line in input_file_reader:
        line = line.rstrip('\n')
        line = line.rstrip(' ')
        get_kmers = GetKmers(line, int(k_length),1)
        kmer_set.append(get_kmers.run_get_kmers())
kmer_set = [kmers for kmer_list in kmer_set for kmers in kmer_list]

#Segunda etapa, definir os nós e edges
get_nodes = GetNodes(kmer_set, int(k_length))
nodes_dict = get_nodes.run_get_nodes()
print(nodes_dict)
#Terceira etapa, definir o caminho de conexoes dos nós