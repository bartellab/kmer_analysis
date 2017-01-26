import sys
import time

import find_kmer_helpers as helpers

if __name__ == '__main__':
    """Generates random sequences and tests speed against naive method.
    Checks that both methods produce the same results.
    """
    kmer_length = 8
    seq_length = 25
    reps = 10

    kmers = helpers.generate_kmers(kmer_length)
    seqs = [helpers.generate_random_seq(seq_length) for i in range(seq_length)]

    time naive method
    t0 = time.time()
    for i in range(reps):
        kmer_dict_sequential = helpers.sequential(kmers, seqs, kmer_length, seq_length)
    print 'average time for sequential: {}'.format((time.time() - t0)/reps)

    # time trie method
    t0 = time.time()
    for i in range(reps):
        kmer_dict_trie = helpers.count_kmers(seqs, kmer_length, seq_length)
    print 'average time for trie: {}'.format((time.time() - t0)/reps)
    kmer_dict_trie.pop('')

    # checks if they return the same dictionary
    print helpers.dicts_equal(kmer_dict_sequential, kmer_dict_trie)

