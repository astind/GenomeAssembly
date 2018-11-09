import format_data
import gen_contigs
import sys
import math


def do_sample(filename, limit):
    all_contigs = {}
    reads = format_data.extact_from_fasta(filename)
    assembly_not_found = True
    errors = True
    if limit == 'na':
        errors = False
    read_length = len(reads[0])
    i = read_length
    #best_contigs_len = math.inf
    #best_contigs = None
    while assembly_not_found:
        if i == 1:
            break
        kmers = format_data.split_into_kmers(reads, i)
        if errors:
            final_kmers = check_errors(kmers, int(limit))
            graph = gen_contigs.db_graph(final_kmers)
        else:
            graph = gen_contigs.db_graph(set(kmers))
        degrees = gen_contigs.calc_degrees(graph)
        contigs = gen_contigs.find_non_branching(graph, degrees)
        all_contigs[i] = contigs
        #if len(contigs) < best_contigs_len and len(contigs) != 0:
            #best_contigs_len = len(contigs)
            #best_contigs = contigs
        i -= 1
    return all_contigs


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


def run_tests(filename, error_limit):
    contigs = do_sample(filename, error_limit)
    longest_contig = 0
    best_set = None
    kmer_size = None
    for kmer, contig in contigs.items():
        for seq in contig:
            length = len(seq)
            if length > longest_contig:
                longest_contig = length
                best_set = contig
                kmer_size = kmer
    return longest_contig, best_set, kmer_size


def format_output(file, longest_contig, best_set, kmer_size, error_limit):
    file.write('Kmer size: ' + str(kmer_size) + '\n')
    if error_limit != 'na':
        file.write('Error limit: ' + str(error_limit) + '\n')
    file.write('Longest contig size: ' + str(longest_contig) + '\n')
    file.write('Number of contigs: ' + str(len(best_set)) + '\n')
    average = 0
    file.write('Best Set:' + '\n')
    for x in best_set:
        file.write(x + '\n')
        average += len(x)
    average = average / len(best_set)
    file.write('Average Contig size: ' + str(average) + '\n')


def print_results():
    with open('student_assembler_results.txt', "w") as file:
        file.write('Gennome Assembly Project Results: \n \n')
        file_list = ['synthetic.example.noerror.small.fasta', 'synthetic.noerror.small.fasta', 'synthetic.noerror.large.fasta',
                     'real.error.small.fasta', 'real.error.large.fasta']
        error_list = ['na', 'na', 'na', 2, 3]
        i = 0
        for f in file_list:
            file.write(f + '\n')
            l_c, b_s, k_s = run_tests(f, error_list[i])
            format_output(file, l_c, b_s, k_s, error_list[i])
            i += 1
            file.write('\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run_all':
        print_results()
    else:
        contigs = do_sample(sys.argv[1], sys.argv[2])
        longest_contig = 0
        best_set = None
        kmer_size = None
        for kmer, contig in contigs.items():
            for seq in contig:
                length = len(seq)
                if length > longest_contig:
                    longest_contig = length
                    best_set = contig
                    kmer_size = kmer

        print('Longest contig size: ' + str(longest_contig))
        print('Number of contigs: ' + str(len(best_set)))
        average = 0
        for x in best_set:
            average += len(x)
        average = average / len(best_set)
        print('Average Contig size: ' + str(average))
        print('Kmer size: ' + str(kmer_size))
        print('Best Set:')
        print(best_set)

