from itertools import product
from matplotlib import pyplot as plt

protein_dic = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}



def determine_likelihood_of_mutated_protein(p):
    total_likelihood = 0.0
    l = list(product('ACGT', repeat=3))
    for orig in l:
        for mutated in l:
            orig_str = "".join(orig)
            mut_str = "".join(mutated)
            if protein_dic[orig_str] != protein_dic[mut_str]:
                num_same = 0
                if orig_str[0] == mut_str[0]: num_same += 1
                if orig_str[1] == mut_str[1]: num_same += 1
                if orig_str[2] == mut_str[2]: num_same += 1
                num_diff = 3 - num_same
                total_likelihood += (1 - p) ** num_same * (p/3) ** num_diff * (1.0 / 64)
    return total_likelihood


if __name__ == "__main__":
    ps = [1.0 * i / 100 for i in range(101)]
    ls = [determine_likelihood_of_mutated_protein(p) for p in ps]
    plt.plot(ps, ls)
    plt.grid(linestyle='--')
    plt.xlabel('Nucleotide mutation rate (given)')
    plt.ylabel('AA mutation rate in a codon (computed)')
    plt.savefig('p_aa_vs_p_nt.pdf')
