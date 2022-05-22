#!/usr/bin/python3
# -*- coding: utf-8 -*-

def get_nodes(kmer_set, k):
    node_dict = {}
    for kmer in kmer_set:
        prefix = kmer[:k-1]
        sufix = kmer[1:]
        node_dict.setdefault(prefix, []).append(sufix)
    return node_dict

class GetNodes():
    def __init__(self, kmer_set, k):
        self.kmer_set = kmer_set
        self.k = k
    
    def run_get_nodes(self):
        return get_nodes(self.kmer_set, self.k)