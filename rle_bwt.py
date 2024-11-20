#!/usr/bin/env python3

"""
DESCRIPTION:
    Template code for the BWT assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Name and student ID here!>
"""

import argparse



# Implement the following functions.
# Replace "raise NotImplementedError" with your own code.

def read_fasta(filename):
    sequence = []
    header_found = False
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header_found: 
                    break
                header_found = True 
                continue
            if not header_found: 
                raise ValueError("Invalid FASTA: missing header line")
            if not line.isalpha():
                raise ValueError(f"Invalid FASTA: non-sequence line detected: {line}")
            sequence.append(line)
    
    if not sequence:
        raise ValueError("Invalid FASTA: no sequence data found")
    
    seq = ''.join(sequence)
    
    return seq

##########################################################

def string_rotations(seq):
    rotations = [seq[i:] + seq[:i] for i in range(len(seq))]
    return rotations
    
def bwt(rotations):
    sorted_rotations = sorted(rotations)
    last_column = [rotation[-1] for rotation in sorted_rotations]
    bwt_seq = ''.join(last_column)
    return bwt_seq  

##########################################################

def rle(seq):
    if not seq:
        return ""
    
    encoded = []
    count = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            count += 1
        else:
            encoded.append(seq[i - 1] + str(count))
            count = 1

    encoded.append(seq[-1] + str(count))
    rle_encoded = ''.join(encoded)
    return rle_encoded
    
def rle_invert(rle_seq):
    if not rle_seq:
        return ""
    
    decoded = []
    i = 0

    while i < len(rle_seq):
        char = rle_seq[i]
        i += 1
        count = 0

        while i < len(rle_seq) and rle_seq[i].isdigit():
            count = count * 10 + int(rle_seq[i])
            i += 1

        decoded.append(char * count)
    rle_decoded = ''.join(decoded)
    return rle_decoded
        
###########################################################

def compute_rank_vector(bwt_seq):
    rank = {}
    rank_vector = []

    for char in bwt_seq:
        if char not in rank:
            rank[char] = 0
        rank_vector.append(rank[char])
        rank[char] += 1

    return rank_vector

def compute_f_map(bwt_seq):
    f_column = sorted(bwt_seq)
    f_map = {}

    for i, char in enumerate(f_column):
        if char not in f_map:
            f_map[char] = i

    return f_map


def bwt_invert(bwt_seq, rank, f_map):
    n = len(bwt_seq)
    original = [''] * n 

    idx = bwt_seq.index('$')
    
    for i in range(n):
        original[n - i - 1] = bwt_seq[idx] 
        idx = f_map[bwt_seq[idx]] + rank[idx]

    return ''.join(original)
##########################################################



# Code for testing:

def main():
    parser = argparse.ArgumentParser()
    file_or_stdin = parser.add_mutually_exclusive_group(required=True)
    file_or_stdin.add_argument('file', nargs='?', type=str, help='FASTA file to load')
    file_or_stdin.add_argument('--stdin', action='store_true', help='Read string as a single line from STDIN')
    file_or_stdin.add_argument(
            '--test',
            choices=['read_fasta', 'string_rotations', 'bwt',
                     'rle', 'rle_invert',
                     'compute_rank_vector', 'compute_f_map', 'bwt_invert',
                     'full'],
            help='Single function to test (codegrade)')
    args = parser.parse_args()

    if not args.test:
        rle_inv_seq = inv_seq = None
        try:
            # Calling the functions
            if args.stdin:
                print('Reading string from STDIN')
                seq = input() + '$'
            else:
                print('Reading ', args.file)
                seq = read_fasta(args.file) + '$'
            print('Sequence:')
            print(seq)

            rotations = string_rotations(seq)
            print('Computed rotations')

            bwt_seq = bwt(rotations)
            print('BWT of sequence:')
            print(bwt_seq)

            rle_seq = rle(bwt_seq)
            print('RLE of BWT:')
            print(rle_seq)

            rle_inv_seq = rle_invert(rle_seq)
            print('RLE inverted:')
            print(rle_inv_seq)

            rank = compute_rank_vector(rle_inv_seq)
            print('Computed rank vector')

            f_map = compute_f_map(rle_inv_seq)
            print('Computed F column map')

            inv_seq = bwt_invert(rle_inv_seq, rank, f_map)
            print('BWT inverted:')
            print(inv_seq)
        except NotImplementedError:
            print('(Exercise is unfinished.)')
        finally:
            # Check the resulting strings, even if not all tasks were finished
            if rle_inv_seq is not None:
                if rle_inv_seq == bwt_seq:
                    print('RLE inversion matches original BWT.')
                else:
                    print('RLE inversion does NOT match original BWT.')

            if inv_seq is not None:
                if inv_seq == seq:
                    print('BWT inversion matches original sequence.')
                else:
                    print('BWT inversion does NOT match original sequence.')
    else:
        # DO NOT CHANGE CODE BELOW -- NECESSARY FOR CODEGRADE

        inp = input()
        if args.test == 'read_fasta':
            print(read_fasta(inp))
        elif args.test == 'string_rotations':
            for rot in string_rotations(inp):
                print(rot)
        elif args.test == 'bwt':
            print(bwt(string_rotations(inp)))
        elif args.test == 'rle':
            print(rle(inp))
        elif args.test == 'rle_invert':
            print(rle_invert(inp))
        elif args.test == 'compute_rank_vector':
            for r in compute_rank_vector(inp):
                print(r)
        elif args.test == 'compute_f_map':
            for c, v in compute_f_map(inp).items():
                print(c, v)
        elif args.test == 'bwt_invert':
            ranks = compute_rank_vector(inp)
            f_map = compute_f_map(inp)
            print(bwt_invert(inp, ranks, f_map))
        elif args.test == 'full':
            inp = rle_invert(rle(bwt(string_rotations(inp))))
            print(bwt_invert(inp, compute_rank_vector(inp), compute_f_map(inp)))

if __name__ == '__main__':
    main()
