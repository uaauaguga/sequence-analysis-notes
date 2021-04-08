# Toy implementations for some sequence analysis task
This repo contains toy implemantations for some bioinformatics algorithms
- Exact matching of strings
  - KMP algorithm
  - Z algorithm
    - See animation [here](https://personal.utdallas.edu/~besp/demo/John2010/z-algorithm.htm)
- Finite state automaton and regular expression
- Some dynamic programming examples
  - Longest common substring
  - Local and global pairwise sequence alignment
  - Profile HMM
  - RNA folding 
- Motif finding
  - Random projection
  ```bash
  g++ src/em-motif.cpp src/utils.cpp src/motifs.cpp  -o bin/em-motif
  ```
  - Expectation Maximization
  - Gibbs sampling



- Reference
- <https://web.stanford.edu/class/cs97si/>
