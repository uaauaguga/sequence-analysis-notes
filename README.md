# Toy implementations for some sequence analysis task
This repo contains toy implemantations for some bioinformatics algorithms
### Sequence comparison
- Exact matching
  - KMP algorithm
  - Z algorithm
    - See animation [here](https://personal.utdallas.edu/~besp/demo/John2010/z-algorithm.htm)
  - Finite state automaton and regular expression
- Gapped alignment
  - Longest common substring
  - Local and global pairwise sequence alignment
  - BWT and FM index
    - <https://people.unipmn.it/manzini/papers/mfcs99x.pdf>
  - Suffix tree / Suffix trie / Suffix array
  - Profile HMM
### Sequence assembly
- DeBruijn
### RNA analysis
- Sequence folding
- SCFG 
### Motif finding
- Random projection
```bash
  g++ src/em-motif.cpp src/utils.cpp src/motifs.cpp  -o bin/em-motif
```
- Expectation Maximization
- Gibbs sampling



- Reference
- <https://web.stanford.edu/class/cs97si/>
