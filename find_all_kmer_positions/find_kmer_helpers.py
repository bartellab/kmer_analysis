import itertools as it
import numpy as np

def make_trie(kmer_length, seq_length):
    """Recursive make a full trie based on the given kmer length.
    Each node has a child called "stop" that stores position counts.
    """
    if kmer_length == 0:
        return {'stop': np.zeros(seq_length)}
    else:
        new_dict = {'stop': np.zeros(seq_length)}
        for nt in ['A','C','G','T']:
            new_dict[nt] = make_trie(kmer_length - 1, seq_length) 
        return  new_dict

def add_to_trie(seq, trie, pos):
    """Loops through characters of a sequence and recursively update trie counts."""
    trie['stop'][pos] += 1
    if len(seq) > 0:
        char = seq.pop(0)
        add_to_trie(seq, trie[char], pos)

def count_trie(trie, current_kmer):
    """Traverses trie and returns a dictionary with all the substring and position counts"""
    kmer_dict = {current_kmer: trie['stop']}
    if len(trie) > 1:
        for nt in ['A','C','G','T']:
            kmer_dict.update(count_trie(trie[nt], current_kmer + nt))

    return kmer_dict

def count_kmers(seqs, kmer_length, seq_length):
    """Given a list of sequences:
    1. Make a full trie
    2. Add each sequence to the trie
    3. Aggregate counts into a dictionary and return dict
    """
    mytrie = make_trie(kmer_length, seq_length)

    for seq in seqs:
        for pos in range(len(seq)):
            subseq = list(seq[pos:min(len(seq),pos+kmer_length)])
            add_to_trie(subseq, mytrie, pos)

    return count_trie(mytrie, '')


def sequential(kmers, seqs, kmer_length, seq_length):
    """Naive method that just loops through all the sequences and updates counts
    """
    master_kmer_dict = {kmer: np.zeros(seq_length) for kmer in kmers}
    for seq in seqs:
        for i in range(seq_length):
            for j in range(i+1, min(seq_length, i+kmer_length) + 1):
                subseq = seq[i:j]
                if subseq in kmers:
                    master_kmer_dict[subseq][i] += 1

    return master_kmer_dict

### some utility functions ###

def dicts_equal(dict1, dict2):
    """Checks if two dictionaries are identical"""
    keys1 = sorted(dict1.keys())
    keys2 = sorted(dict2.keys())

    if keys1 != keys2:
        print "not same keys"
        return False

    for key in keys1:
        if len(dict1[key] > 0):
            if min(dict1[key] == dict2[key]) == 0:
                print dict1[key], dict2[key], key
                return False

    return True


def generate_random_seq(length):
    """Generates a random nucleotide sequence"""
    nts = ["A","C","G","T"]
    seq = np.random.choice(nts, size=length, replace=True)
    return ''.join(seq)


def generate_kmers(kmer_length):
    """Generates all kmers of length 'kmer_length' or shorter"""
    kmers = []
    for length in range(1, kmer_length+1):
        kmers += ["".join(kmer) for kmer in list(it.product(["A","C","G","T"],repeat=length))]

    return kmers
