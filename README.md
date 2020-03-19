## Filter a .ab1 or .abi file on phred quality scores

Very small script which takes an input .ab1 *directory* and filters then concatenates all sequences into a single fasta file.

## Usage

`filter_ab1.py [-h] [--phred PHRED] [--replacement REPLACEMENT] path`

Mandatory argument:
path                  Input file(s) path.
Optional arguments:
-h, --help            show this help message and exit
--phred PHRED, -p PHRED
		      Phred quality score cutoff, default 5
--replacement REPLACEMENT, -r REPLACEMENT
                        The replacement letter for nucleotides below quality
                        threshold. Default=N
