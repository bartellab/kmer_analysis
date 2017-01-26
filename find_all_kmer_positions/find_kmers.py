from optparse import OptionParser
import os
import sys
import time

import concurrent.futures
import numpy as np

import find_kmer_helpers as helpers


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="infile",
                      help="FILE with sequence data")
    parser.add_option("-c", "--concurrent",
                      action="store_true", dest="concurrent", default=False,
                      help="run concurrently")

    (options, args) = parser.parse_args()

    print "Analyzing sequences from {}".format(options.infile)

    len_data_chunk = 500
    max_lines = 2
    seq_length = 25
    kmer_length = 8

    # make a list of kmers up to 8 long
    kmers = helpers.generate_kmers(kmer_length)

    # make a dictionary to store counts
    kmer_dict = {kmer: np.zeros(seq_length) for kmer in kmers}

    # run in parallel if specified
    print "Running in parallel"
    if options.concurrent:
        executor = concurrent.futures.ProcessPoolExecutor()
        futures = []

        # iterate through lines in files and split off jobs
        seqs = []
        i = 0
        with open(options.infile, 'r') as f:
            for line_num, line in enumerate(f):
                seqs.append(line[:seq_length])
                i += 1
                if i == len_data_chunk:
                    # add job with collected sequences
                    futures.append(executor.submit(helpers.count_kmers,
                                                   seqs,
                                                   kmer_length,
                                                   seq_length))

                    # add a sleep to prevent the executor from getting clogged
                    time.sleep(0.0001)

                    # reset sequences
                    seqs = []
                    i = 0

                # if i == max_lines:
                #     break

            # analyze last chunk
            futures.append(executor.submit(helpers.count_kmers,
                                           seqs,
                                           kmer_length,
                                           seq_length))

        print "all jobs submitted"
        # as jobs are completed, add the information
        for future in concurrent.futures.as_completed(futures):
            new_dict = future.result()
            for kmer in kmers:
                kmer_dict[kmer] += new_dict[kmer]

        executor.shutdown()

    # otherwise, run normally
    else:
        seqs = []
        i = 0
        with open(options.infile, 'r') as f:
            for line_num, line in enumerate(f):
                seqs.append(line[:seq_length])
                i += 1
                if i == len_data_chunk:
                    new_dict = helpers.count_kmers(seqs, kmer_length, seq_length)
                    for kmer in kmers:
                        kmer_dict[kmer] += new_dict[kmer]

                    seqs = []
                    i = 0

                # if i == max_lines:
                #     break

            # analyze last chunk
            new_dict = helpers.count_kmers(seqs, kmer_length, seq_length)
            for kmer in kmers:
                kmer_dict[kmer] += new_dict[kmer]

    print kmer_dict
