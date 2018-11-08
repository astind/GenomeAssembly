import format_data
import gen_contigs
import sys
import math


def do_sample(filename, limit):
    reads = format_data.extact_from_fasta(filename)
    assembly_not_found = True
    read_length = len(reads[0])
    i = read_length
    best_contigs_len = math.inf
    best_contigs = None
    while assembly_not_found:
        if i == 1:
            break
        kmers = format_data.split_into_kmers(reads, i)
        final_kmers = check_errors(kmers, limit)
        graph = gen_contigs.db_graph(final_kmers)
        degrees = gen_contigs.calc_degrees(graph)
        contigs = gen_contigs.find_non_branching(graph, degrees)
        if len(contigs) < best_contigs_len:
            best_contigs_len = len(contigs)
            best_contigs = contigs
        i -= 1
    return best_contigs, best_contigs_len


def check_errors(kmers, limit):
    counts = {}
    for kmer in kmers:
        if kmer not in counts.keys():
            counts[kmer] = 1
        else:
            counts[kmer] = counts[kmer] + 1
    new_kmers = []
    for key, count in counts.items():
        if count > limit:
            new_kmers.append(key)
    return new_kmers


if __name__ == '__main__':
    contigs, contigs_len = do_sample(sys.argv[1], sys.argv[2])
    print(contigs)
    print(contigs_len)

