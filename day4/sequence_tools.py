# Reverse complementation of sequences.
# Dictionary mapping bases to their Watson-Crick pairing partners.
COMPLEMENTS = {'N': 'N',
               'A': 'T',
               'C': 'G',
               'G': 'C',
               'T': 'A'}
# Get the set of base characters allowed from the keys of the dict.
ALLOWED_BASES = set(COMPLEMENTS.keys())


def reverse_comp(seq):
    """Reverse complement a sequence."""
    upper_seq = seq.upper()  # upper case
    upper_bases = set(upper_seq)
    if not ALLOWED_BASES.issuperset(upper_bases):
        print 'Unknown bases in', seq
        return
    
    complemented = [COMPLEMENTS[c] for c in upper_seq]
    complemented.reverse()
    return ''.join(complemented)
    
    
def parse_fasta(filename):
    """Parse a FASTA file given the filename."""
    current_key = None
    current_seq = []
    fasta_dict = {}
    with open(filename) as fastaf:    
        for line in fastaf:
            # If we encounter '>' then we've hit a new sequence.
            if line.startswith('>'):
                # If we've been accumulating data, save it.
                if current_key is not None:
                    fasta_dict[current_key] = ''.join(current_seq)
                
                # The key is the current line stripped of
                # newlines and the '>' character.
                current_key = line.strip().strip('>')
                current_seq = []  # Reset the list for accumulating sequence.
            else:
                # Not a new sequence name, so just accumulate the current sequence.
                current_seq.append(line.strip())
        else:
            # Make sure we get the last sequence into the dict.
            if current_key is not None:
                fasta_dict[current_key] = ''.join(current_seq)
    return fasta_dict


def find_longest_orf(seq):
    """Finds the longest ORF in the sequence.
    
    ORF goes from ATG to an in-frame TAA, TAG or TGA.
    
    Returns a 2-tuple of the index of the first base of ATG 
    and the first base of the stop codon. 
    """
    # Could you edit this code to find all of the ORFs in the 
    # sequence instead of just the longest one?
    seq_upper = seq.upper()
    start_idx = seq_upper.find('ATG')
    if start_idx == -1:
        return -1, -1
    start_frame = start_idx % 3
    
    taa_idx = seq_upper.rfind('TAA')
    taa_frame = taa_idx % 3
    
    tag_idx = seq_upper.rfind('TAG')
    tag_frame = tag_idx % 3
    
    tga_idx = seq_upper.rfind('TGA')
    tga_frame = tga_idx % 3
    
    end_idx = -1
    end_idxs = [taa_idx, tag_idx, tga_idx]
    for idx in end_idxs:
        frame = idx % 3
        if frame == start_frame and idx > end_idx:
            end_idx = idx
    
    if start_idx > end_idx:
        return -1, -1
    return start_idx, end_idx


def get_longest_orf(seq):
    """Gets the sequence of the ORF from this sequence.
    
    Calls find_longest_orf to get the index of the start 
    and stop codons. 
    """
    start_idx, end_idx = find_longest_orf(seq)
    if start_idx == -1 or end_idx == -1:
        return ''
    return seq[start_idx: end_idx+3]


def get_orf_dict(fasta_dict):
    """Converts a dictionary of sequences into a dictionary of ORFs."""
    d = {}
    for gene_name, gene_seq in fasta_dict.iteritems():
        d[gene_name] = get_longest_orf(gene_seq)
    return d
