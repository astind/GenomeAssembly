

def extact_from_fasta(filename):
    with open(filename, "r") as file:
        reads = []
        for line in file:
            if line[0] != '>':
                reads.append(line.strip())
    reads.sort()
    return reads


def split_into_kmers(reads, k_size):
    read_length = len(reads[0])
    if k_size > read_length:
        return None
    else:
        kmers = []
        r = (read_length - k_size) + 1
        for read in reads:
            for i in range(r):
                kmers.append(read[i:i+k_size])
        kmers.sort()
        return set(kmers)


# i don't think this works
def check_for_errors(kmers):
    counts = {}
    for kmer in kmers:
        if kmer not in counts.keys():
            counts[kmer] = 0
        else:
            counts[kmer] = counts[kmer] + 1
    for _, count in counts.items():
        if count == 0:
            return False, counts
    return True, counts


if __name__ == '__main__':
    # test
    reads = extact_from_fasta('synthetic.noerror.large.fasta')
    kmers = split_into_kmers(reads, 10)
    print(kmers)
