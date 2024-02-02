#!/usr/bin/env python3

import magnumopus

primer_file = "data/rpoD.fna"
assembly_1 = "data/Pseudomonas_aeruginosa_PAO1.fna"
assembly_2 = "data/Pseudomonas_protegens_CHA0.fna"
max_amp_size = 2000

amplicon_1 = magnumopus.ispcr(primer_file, assembly_1, max_amp_size)
amplicon_2 = magnumopus.ispcr(primer_file, assembly_2, max_amp_size)

print(amplicon_1)
print(amplicon_2)