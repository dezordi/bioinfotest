#!/usr/bin/python3
# -*- coding: utf-8 -*-

def get_kmers(sequence, k, step):
    return [sequence[x:x+k].upper() for x in range(0, len(sequence) - k +1, step)]

class GetKmers():
    def __init__(self, sequence, k, step):
        self.sequence = sequence
        self.k = k
        self.step = step
    
    def run_get_kmers(self):
        return get_kmers(self.sequence, self.k, self.step)