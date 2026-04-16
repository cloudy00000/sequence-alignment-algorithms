# Sequence Alignment Algorithms

This repository contains Python implementations of basic sequence alignment algorithms developed for learning and portfolio purposes in bioinformatics.

The project includes:

* Hamming distance and similarity
* Substitution matrix parsing
* Needleman-Wunsch algorithm for global alignment
* Smith-Waterman algorithm for local alignment

## Algorithms included

### Hamming distance

Implemented in `hamming.py`.

This script defines a class for comparing two sequences of equal length and computing:

* Hamming distance
* a simple similarity score based on matches and mismatches

### Substitution matrix

Implemented in `submatrix.py`.

This script defines a class that reads a substitution matrix from a text file and stores it in a dictionary for score lookup.

### Global alignment

Implemented in `global_aln.py`.

This script defines a class for pairwise global sequence alignment using the Needleman-Wunsch algorithm.

### Local alignment

Implemented in `local_aln.py`.

This script defines a class for pairwise local sequence alignment using the Smith-Waterman algorithm.

## Repository structure

```text
.
├── hamming.py
├── submatrix.py
├── global_aln.py
├── local_aln.py
└── README.md
```

## Requirements

This project uses standard Python and does not require external libraries.

## Basic usage

The global and local alignment scripts require:

* two input sequences
* a substitution matrix
* a gap penalty

## Example workflow

* Load a substitution matrix from a text file
* Create a sequence pair
* Run the desired alignment method
* Retrieve the aligned sequences and alignment score

## Notes

* Hamming distance can only be computed for sequences of equal length.
* The alignment methods assume valid sequence strings and a correctly formatted substitution matrix.
* The project is intended for educational use and algorithm practice.

## Context

This repository was developed as part of my learning in bioinformatics and includes cleaned and reorganized implementations of classical sequence analysis methods.

## Suggested reading order

1. `README.md`
2. `hamming.py`
3. `submatrix.py`
4. `global_aln.py`
5. `local_aln.py`
6. `example_usage.py`