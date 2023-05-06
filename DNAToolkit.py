#DNA toolkit file
import collections
from typing import Counter
from structures import *
# Check the sequence to make sure it is a DNA String
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

#count the nucleotide frequency
def countNucFrequency(seq):
    tmpFreqDict = {'A':0,'C':0,'G':0,'T':0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

#########################################################
#TRANSCRIPTION
#DNA -> RNA Transcription
def transcription(seq):
    return seq.replace('T', 'U')

#DNA Reverse compliment
def reverse_compliment(seq):
    #return ''.join([DNA_ReverseCompliment[nuc] for [nuc] in seq])

    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)

#GC content calculation
def gc_content(seq):
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

#GC content in subset
def gc_content_subset(seq, k=20):
    #####GC content in a DNA/RNA sub-sequence length k, k=20 by default###
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res

#TRANSLATION
#translate DNA sequence into AA sequence, init_pos indicates the reading frame
def translate_seq(seq, init_pos=0): #init_pos=0 means the default value is set to the start of the sequence
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]

#Codon usage statistics
def codon_usage(seq, AA): #You have to specify which AA you are looking for
    ###provides the freq of each codon encoding a given AA in a DNA sequence
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == AA:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict

# Generate the 6 reading frames of a DNA sequence, including reverse compliment
def gen_reading_frames(seq):
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_compliment(seq), 0))
    frames.append(translate_seq(reverse_compliment(seq), 1))
    frames.append(translate_seq(reverse_compliment(seq), 2))
    return frames

#Protein search within a reading frame
def proteins_from_rf(aa_seq):
    ### Compute all possible proteins in an amino acid seq and return a list of possible proteins
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "STOP":
            # STOP accumulating aa if stop codon was found
            if current_prot:
                 for p in current_prot:
                       proteins.append(p)
                 current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                 current_prot.append("")
            for i in range(len(current_prot)):
                 current_prot[i] += aa
    return proteins

def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)

    if ordered:
        return sorted(res, key=len, reverse=True) # sort by length, from the longest string to the shortest string
    
    res = [p.replace("[", "").replace("]", "").replace("'", "") + "\n" for p in res]
    return res

def longest_protein_from_orfs(seq, startReadPos=0, endReadPos=0):
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)

    longest = ""
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            if len(p) > len(longest):
                longest = p

    longest = longest.replace("[", "").replace("]", "").replace("'", "")
    return longest

# Generate all RFs
# Extract proteins
# Return a list sorted/unsorted







