# Some sequence analysis task
- This repo contains **toy** implemantations for some bioinformatics algorithms, **not** tested and optimized for practical usage

## Sequence comparison
### Exact matching
- KMP algorithm (learn from mismatch)
- Z algorithm
  - See animation [here](https://personal.utdallas.edu/~besp/demo/John2010/z-algorithm.htm)

- Finite state automaton and regular expression

- trie (re**trie**val)
  - also called digital tree or prefix tree
  - Check whether query word appears in a set of words
  
```bash
  g++ src/utils.cpp src/exact-match.cpp src/trie.cpp -o bin/trie
```
- [AC-automaton](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm)
  - <https://github.com/cjgdev/aho_corasick>
  - <http://web.stanford.edu/class/archive/cs/cs166/cs166.1166/lectures/02/Small02.pdf>  

- [radix tree](https://en.wikipedia.org/wiki/Radix_tree)
  - A compressed version of trie
  - Also called radix trie or compact prefix tree

- Suffix tree 
  - Very verstile data structure in sequence analysis
  - <https://web.stanford.edu/~mjkay/suffix_tree.pdf>

- Suffix array
  - Simialr functionality as siffix, with less memory overhead
  - <https://www.cs.unc.edu/~prins/Classes/555/Media/Lec17.pdf>
  - Mutiple method for suffix array construction <http://www.cas.mcmaster.ca/~bill/best/algorithms/07Taxonomy.pdf>


-  FM index
  - Suffix array is still memory intensive
  - [BWT](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform)
    - Burrows-Wheeler block-sorting transform
    - Orinially designed for data compression 
      - <https://www.cs.unc.edu/~prins/Classes/555/Media/Lec17.pdf>
  - <https://people.unipmn.it/manzini/papers/mfcs99x.pdf>
  - <https://www.labri.fr/perso/ruricaru/bioinfo_master2/cours3.pdf>
  
- Compressed suffix array
  - Similar to FM index
  - <https://www.cs.cmu.edu/~dga/csa.pdf>

### Gapped alignment
- Longest common substring
  - common substring: allow gaps
  - common subsequence: not allow gaps
- Local and global pairwise sequence alignment

- Profile HMM

## Sequence assembly

- DeBruijn Graph

## RNA analysis
- See <http://rna.informatik.uni-freiburg.de/Teaching/index.jsp>
- MFE Folding
  - Without pseudoknot
    - Nussinov Folding: maximize number of base pairing

   ```bash
   g++ src/nussinov-folding.cpp -o bin/nussinov-folding
   ```

  - With pseudoknot
- Partition fucntion
  - McCaskill's algorithm
    - <https://onlinelibrary.wiley.com/doi/abs/10.1002/bip.360290621>
  - Pairing probability
  - MEA structure
- SCFG 
- Structural alignment

## Motif finding
- Random projection
```bash
  g++ src/em-motif.cpp src/utils.cpp src/motifs.cpp  -o bin/em-motif
```
- Expectation Maximization
- Gibbs sampling



- Reference
- <https://web.stanford.edu/class/cs97si/>
- <https://www.comp.nus.edu.sg/~ksung/algo_in_bioinfo/>
