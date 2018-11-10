import sys

def main():
    with open("AbyssResults/synthex-3.fa", "r") as file:
        reads = []
        cnt = 0
        first = True
        for line in file:
            if line[0] == '>' and first is False:
                reads.append(cnt)
                cnt = 0
            elif line[0] != '>':
                first = False
                cnt = cnt + len(line.strip())
        reads.append(cnt)
    average = sum(reads) / len(reads)
    print("length list: " + str(len(reads)))
    print("average: " + str(average))
    print("max: " + str(max(reads)))
    return reads


if __name__ == '__main__':
    main()
