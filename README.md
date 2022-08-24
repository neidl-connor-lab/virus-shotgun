# characterize-virus
Characterize virus stocks

## Setup

There are several tools you will need to run `pipeline.sh`. By first running `setup.sh`, we'll set up the project directory!

User MUST provide:
- virus FASTA (all one line?)
- virus GTF/CSV

Setup steps in no particular order:
- Set up Kraken2 -- install if necessary; always install standard db under pipeline/kraken2
- make Bowtie2 index -- created under pipeline/bowtie2
- install vodka -- created under pipeline/vodka

## Pipeline steps

1. check for contaminants (Kraken2)
2. Align to virus (Bowtie2)
3. Assess coverage (samtools)
4. Call SNVs (LoFreq)
5. Annotate SNVs (R; from GTF?)
6. assess DVGs (VODKA)