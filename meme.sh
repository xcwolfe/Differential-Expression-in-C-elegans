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
