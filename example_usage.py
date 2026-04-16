import hamming
import submatrix
import global_aln as gb
import local_aln as lc


def main():
    print("HAMMING DISTANCE EXAMPLE")
    hamming_pair = hamming.AlnSeq("ACT", "AGT")
    print("Sequence 1:", hamming_pair.seq1)
    print("Sequence 2:", hamming_pair.seq2)
    print("Hamming distance:", hamming_pair.compute_hamming_distance())
    print("Hamming similarity:", hamming_pair.compute_hamming_similarity())
    print()

    sub_matrix = submatrix.SubstitutionMatrix("TTM.txt")

    print("GLOBAL ALIGNMENT EXAMPLE")
    global_pair = gb.SeqPair("TCA", "TA")
    aln1, aln2, score = global_pair.nw(sub_matrix, -2)
    print("Aligned sequence 1:", aln1)
    print("Aligned sequence 2:", aln2)
    print("Global alignment score:", score)
    print()

    print("LOCAL ALIGNMENT EXAMPLE")
    local_pair = lc.SeqPair("TGA", "GA")
    aln1, aln2, score = local_pair.sw(sub_matrix, -2)
    print("Aligned sequence 1:", aln1)
    print("Aligned sequence 2:", aln2)
    print("Local alignment score:", score)


if __name__ == "__main__":
    main()