# Run this script to count the codons in the sequences extracted from "12. Identifying stop codons and reading frames in intron retention and A3S events.Rmd"
```
$ module load samtools/ bedtools/
```
# insert your .tsv file (containing coordinates of interest) below:
```
$ tail -n +2 A3S_coordinates.tsv > A3S_coordinates_noheader.tsv
```
```
$ awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $1, ".", $3}' A3S_coordinates_noheader.tsv > A3S_coordinates.bed
```

# Use bedtools slop to extend the region by 1 base on each side and then extract the sequence
```
$ bedtools slop -i A3S_coordinates.bed -g genome/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai -l 1 -r 1 -s | \
$ bedtools getfasta -fi genome/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa -bed - -fo A3S_sequences.fa

# eliminate the final character of each sequence since we "slopped" by one on each side:
$ sed '/^>/!s/.$//' A3S_sequences.fa > A3S_sequences.tsv
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
