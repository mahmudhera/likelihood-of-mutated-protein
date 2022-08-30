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

def determine_likelihood_of_mutated_protein_given_codon(codon, p):
    # given: codon: 3-mer string
    # p: mutation rate
    l = list(product('ACGT', repeat=3))
    likelihood_nonsyn = 0.0
    likelihood_syn = 0.0
    for mutated in l:
        mut_codon = "".join(mutated)

        num_same = 0
        if codon[0] == mut_codon[0]: num_same += 1
        if codon[1] == mut_codon[1]: num_same += 1
        if codon[2] == mut_codon[2]: num_same += 1
        num_diff = 3 - num_same

        if protein_dic[codon] != protein_dic[mut_codon]:
            likelihood_nonsyn += (1 - p) ** num_same * (p/3) ** num_diff
        elif codon != mut_codon:
            likelihood_syn += (1 - p) ** num_same * (p/3) ** num_diff

    return likelihood_syn, likelihood_nonsyn

if __name__ == "__main__":
    """
    ps = [1.0 * i / 100 for i in range(101)]
    ls = [determine_likelihood_of_mutated_protein(p) for p in ps]
    plt.plot(ps, ls)
    plt.grid(linestyle='--')
    plt.xlabel('Nucleotide mutation rate (given)')
    plt.ylabel('AA mutation rate in a codon (computed)')
    plt.savefig('p_aa_vs_p_nt.pdf')
    """

    p = 0.03
    all_codons = [ "".join(x) for x in list(product('ACGT', repeat=3)) ]
    for codon in all_codons:
        syn_mut_prob, nonsyn_mut_prob = determine_likelihood_of_mutated_protein_given_codon(codon, p)
        print(codon, syn_mut_prob, nonsyn_mut_prob, syn_mut_prob+nonsyn_mut_prob, 1.0-(1.0-p)**3)
