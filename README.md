# Some sequence analysis task
- This repo contains **toy implemantations** for some bioinformatics algorithms, not for practical usage

### Sequence comparison
- Exact matching
  - KMP algorithm (learn from mismatch)
  - Z algorithm
    - See animation [here](https://personal.utdallas.edu/~besp/demo/John2010/z-algorithm.htm)
  - [AC-automaton](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm)
    - <https://github.com/cjgdev/aho_corasick>
    - <http://web.stanford.edu/class/archive/cs/cs166/cs166.1166/lectures/02/Small02.pdf>
  - Finite state automaton and regular expression
  - trie (re**trie**val)
    - also called digital tree or prefix tree
    - Check whether query word appears in a set of words
- [radix tree](https://en.wikipedia.org/wiki/Radix_tree)
    - Also called radix trie or compact prefix tree
  - Suffix tree / Suffix trie / Suffix array
  - [BWT](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform) and FM index
    - <https://people.unipmn.it/manzini/papers/mfcs99x.pdf>
- Gapped alignment
  - Longest common substring
  - Local and global pairwise sequence alignment

  - Profile HMM

### Sequence assembly
- DeBruijn Graph

### RNA analysis
- Folding
- SCFG 
- Structural alignment

### Motif finding
- Random projection
```bash
  g++ src/em-motif.cpp src/utils.cpp src/motifs.cpp  -o bin/em-motif
```
- Expectation Maximization
- Gibbs sampling



- Reference
- <https://web.stanford.edu/class/cs97si/>
