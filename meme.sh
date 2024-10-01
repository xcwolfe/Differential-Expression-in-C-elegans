# Running MEME on HPCC:

## Download MEME and associated scripts using installation guidelines found on https://meme-suite.org/meme/doc/download.html and then unpack MEME using guidelines on https://meme-suite.org/meme/doc/install.html?man_type=web

## MEME command lines (example using AVM cassette deltaPSIs > 5) in discriminative mode:

```
$ export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.4:$PATH
$ module load bedtools/
```
```
$ bed2fasta -s -both -o AVM_cassette_MEME.tsv.fa AVM_cassette_MEME.tsv /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
```
```
$ bed2fasta -s -both -o AVM_cassette_MEME.tsv.fa AVM_cassette_MEME.tsv /data/apache-fasta-shuffle-letters AVM_cassette_MEME.tsv.fa shuffled_AVM_cassette_MEME.tsv
```
```
$ meme shuffled_AVM_cassette_MEME.tsv.fa -dna -oc . -nostatus -time 14400 -mod anr -nmotifs 50 -minw 6 -maxw 8 -objfun de -revcomp -markov_order 0
```
```
$ mast meme.xml shuffled_AVM_cassette_MEME.tsv.fa -oc . -nostatus
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
