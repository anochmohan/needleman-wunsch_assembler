#!/usr/bin/env python3

import argparse
import magnumopus

#Read args
parser = argparse.ArgumentParser(description="Perform in-silico PCR on two assemblies and align the amplicons")
parser.add_argument('-1',
                    required=True,
                    dest="ASSEMBLY1",
                    help=("Path to the first assembly file"))
parser.add_argument('-2',
                    required=True,
                    dest="ASSEMBLY2",
                    help="Path to the second assembly file")
parser.add_argument('-p',
                    required=True,
                    dest="PRIMERS",
                    help="Path to the primer file")
parser.add_argument('-m',
                    required=True, 
                    type=int,
                    dest="MAX_AMPLICON_SIZE",
                    help="maximum amplicon size for isPCR")
parser.add_argument('--match',
                    required=True,
                    type=int,
                    metavar="MATCH",
                    help="match score to use in alignment")
parser.add_argument('--mismatch',
                    required=True,
                    type=int,
                    metavar="MISMATCH",
                    help="mismatch penalty to use in alignment")
parser.add_argument('--gap',
                    required=True,
                    type=int,
                    metavar="GAP",
                    help="gap penalty to use in alignment")
#parser.parse_args(['-h', '--help'])

args = parser.parse_args()

# Variable Setup
assembly_1 = args.ASSEMBLY1
assembly_2 = args.ASSEMBLY2
primer_file = args.PRIMERS
maximum_amplicon_size = args.MAX_AMPLICON_SIZE
match = args.match
mismatch = args.mismatch
gap = args.gap

# Call ispcr
assembly1_run = magnumopus.ispcr(primer_file, assembly_1, maximum_amplicon_size)
#print(assembly1_run)

assembly2_run = magnumopus.ispcr(primer_file, assembly_2, maximum_amplicon_size)
#print(assembly2_run)

# Remove header (1 line only) from assembly1_run and assembly2_run as we only need the sequence for Needleman-Wunch
sequence_1 = "\n".join(assembly1_run.split("\n")[1:])
#print(sequence_1)

sequence_2 = "\n".join(assembly2_run.split("\n")[1:])
#print(sequence_2)

# Call Neewdle-Wunch
aligned_1, score_1 = magnumopus.needleman_wunsch(sequence_1, sequence_2, match, mismatch, gap)

# Reverse compliment check
reverse_com_sequence_2 = []
for i in range(len(sequence_2)):
    if sequence_2[i] == 'A':
        reverse_com_sequence_2.append('T')
    elif sequence_2[i] == 'T':
        reverse_com_sequence_2.append('A')
    elif sequence_2[i] == 'C':
        reverse_com_sequence_2.append('G')
    else:
        reverse_com_sequence_2.append('C')

reverse_com_sequence_2 = "".join(reverse_com_sequence_2)[::-1]
#print(reverse_com_sequence_2)

aligned_2, score_2 = magnumopus.needleman_wunsch(sequence_1, reverse_com_sequence_2, match, mismatch, gap)

if score_1 > score_2:
    print("\n".join(aligned_1) + "\n")
    print(score_1)
else:
    print("\n".join(aligned_2) + "\n")
    print(score_2+2)