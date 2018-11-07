import format_data
import gen_contigs


def do_sample():
    reads = format_data.extact_from_fasta('example.data.fasta')
    assmebly_not_found = True
    read_length = len(reads[0])
    i = read_length
    while assmebly_not_found:
        if i == 1:
            return None
        kmers = format_data.split_into_kmers(reads, i)
        graph = gen_contigs.db_graph(kmers)
        degrees = gen_contigs.calc_degrees(graph)
        contigs = gen_contigs.find_non_branching(graph, degrees)
        if len(contigs) > 1:
            i -= 1
        else:
            return contigs


if __name__ == '__main__':
    genome = do_sample()
    print(genome)
