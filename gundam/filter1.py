from gundam import Gundam
import argparse
import csv

P_CUTOFF = 1e-6
INFO_CUTOFF = 10.0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simple information-content and p-val based filter')
    parser.add_argument('kmer_file')
    parser.add_argument('pos_seqs', help='fasta format')
    parser.add_argument('neg_seqs', help='fasta format')

    args = parser.parse_args()

    with open(args.kmer_file) as kmer_f:
        kmers = [[int(i), int(j), int(k), float(p)]
                 for i, j, k, p in csv.reader(kmer_f)
                 if float(p) <= P_CUTOFF]
        G = Gundam( kmers, args.pos_seqs, args.neg_seqs )
        means_idx = []
        orig_len = len(G)
        for i in range(orig_len):
            new_idx = G.add_mean(i)
            if G.info_content(new_idx) >= INFO_CUTOFF:
                means_idx.append( (i, new_idx) )
        print '# original length:', orig_len
        for _, new_idx in means_idx:
            motif = G[new_idx]
            print 'motif: {}, seq ct: {}'.format(motif['repr'], len(motif['pos_seqs']))

