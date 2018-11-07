

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
    if k_size > read_length or read_length % k_size != 0:
        return None
    else:
        kmers = []
        for read in reads:
            split_set = [read[i:i+k_size] for i in range(0, len(read), k_size)]
            for k in split_set:
                kmers.append(k)
        kmers.sort()
        return kmers

# i don't think this works
def check_for_errors(kmers):
    counts = {}
    for kmer in kmers:
        if kmer not in counts.keys():
            counts[kmer] = 0
        else:
            counts[kmer] = counts[kmer] + 1
    base = None
    for _, count in counts.items():
        if base == None:
            base = count
        else:
            if base == count:
                continue
            else:
                return False, counts
    return True, counts


if __name__ == '__main__':
    # test
    reads = extact_from_fasta('synthetic.noerror.large.fasta')
    no_split, _ = check_for_errors(reads)
    if no_split:
        print('no_split')
    else:
        errors = True
        i = 2
        while errors:
            if i == len(reads[0]):
                break
            kmers = split_into_kmers(reads, i)
            if kmers != None:
                res, counts = check_for_errors(kmers)
                if res:
                    errors = False
            i += 1
    print(kmers)
