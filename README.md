# Some sequence analysis task
- This repo contains **toy** implemantations for some bioinformatics algorithms in C++, only for practice, **not** well tested and optimized for practical usage
- Here list some related papers, notes, and other meterials

## Sequence comparison
### Exact matching
- [KMP algorithm](https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm) (learn from mismatch)
- Z algorithm
  - See animation [here](https://personal.utdallas.edu/~besp/demo/John2010/z-algorithm.htm)
  - <https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/zalg.pdf>

- Finite state automaton and regular expression
  - Pattern matching FSA

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
  - 2014, *BIB*, [A bioinformatician’s guide to the forefront of suffix array construction algorithms](https://academic.oup.com/bib/article/15/2/138/212729)
  - Mutiple method for suffix array construction <http://www.cas.mcmaster.ca/~bill/best/algorithms/07Taxonomy.pdf>


-  FM index
  - Suffix array is still memory intensive
  - [BWT](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform)
    - Burrows-Wheeler block-sorting transform
    - Orinially designed for data compression 
    - Some reference
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
  - Scoring of sequence alignment
  - [Needleman-Wunsch global alignment](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) 
  ```bash
  # with affine gap penalty, implemented as a finite state automaton
  g++  src/needleman-wunsch.cpp -o bin/needleman-wunsch
  ```
  - [Smith-Waterman local alignment](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)

- [Dynamic time warping (DTW)](https://en.wikipedia.org/wiki/Dynamic_time_warping)
  - Highly similar to dynamic programming for sequence alignment

### HMM in sequence analysis
- Inference
  - Given the model, calculate probability of an observed sequence
  - Forward algorithm
- Decoding
  - Given model and observed sequence, calculate most probable hidden states
  - viterbi/global decoding and posterior/local decoding
- Training
  - Given obervations, and HMM model structure, calculate model parameters
  - BM training and viterbi training
- Profile HMM

## Sequence assembly

- [DeBruijn Graph](https://en.wikipedia.org/wiki/De_Bruijn_graph)
  - 2008, Genome Research [Velvet: Algorithms for de novo short read assembly using de Bruijn graphs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2336801/pdf/821.pdf)

## RNA analysis
- See <http://rna.informatik.uni-freiburg.de/Teaching/index.jsp>


### Nussinov Folding
- maximize number of base pairing
- 1980,PNAS,[Fast algorithm for predicting the secondary structure of single-stranded RNA](https://www.pnas.org/content/77/11/6309)
- <https://en.wikipedia.org/wiki/Nussinov_algorithm>
```bash
  g++ src/nussinov-folding.cpp -o bin/nussinov-folding
```



### Partition fucntion of RNA structure
- McCaskill's algorithm
  - <https://onlinelibrary.wiley.com/doi/abs/10.1002/bip.360290621>
- Pairing probability

### MEA (**m**aximum **e**xpected **a**ccuracy) folding
- 2006, Chuong B. Do et al, *Bioinformatics*, [CONTRAfold: RNA secondary structure prediction without physics-based models](https://academic.oup.com/bioinformatics/article/22/14/e90/228433)
- 2007, Hisanori Kiryu et al, *Bioinformatics*, [Robust prediction of consensus secondary structures using averaged base pairing probability matrices](https://academic.oup.com/bioinformatics/article/23/4/434/182043)
- 2009, Zhi John Lu et al, *RNA*, [Improved RNA secondary structure prediction by maximizing expected pair accuracy](https://rnajournal.cshlp.org/content/15/10/1805.long)

### [Context free grammer](https://en.wikipedia.org/wiki/Context-free_grammar) (CFG) and sttochastic context free grammer (SCFG)
- CFG and [CYK parsing](https://en.wikipedia.org/wiki/CYK_algorithm)


### Structural alignment
- The RIBOSUM scoring schema
  - params/RIBOSUM90.mat is copied from rsearch-1.1/matrices/RIBOSUM90.mat in <http://eddylab.org/software/rsearch/rsearch.tar.gz>
  - Also see  2003, *BMC Bioinformatics*, [RSEARCH: Finding homologs of single structured RNA sequences](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-4-44)
- Sankoff's simultaneous alignment and folding
  - Dynamic programming for structural alignment
  - 1985, David Sankoff, Society for Industrial and Applied Mathematics, [Simultaneous solution of the RNA folding, alignment and protosequence problems](https://epubs.siam.org/doi/10.1137/0145048)
  - Dynalign
  - Foldalign
- [PMcomp](https://www.tbi.univie.ac.at/RNA/PMcomp/) like algorithm
  - Align structure based on pre-computed base pairing probability matrix
  - 2004, Ivo L. Hofacker et al, *Bioinformatics*, [Alignment of RNA base pairing probability matrices](https://academic.oup.com/bioinformatics/article/20/14/2222/214007)
  ```bash
  g++ src/align-dotplot.cpp src/utils.cpp -o bin/align-dotplot
  bin/align-dotplot --dotplot-A test/test1_dp.ps  --dotplot-B test/test_dp.ps
  ```
  - without pseudoknot
    - [LocaRNA](http://www.bioinf.uni-freiburg.de/Software/LocARNA/)
      - 2007, *Plos Computational Biology*, [Inferring Noncoding RNA Families and Classes by Means of Genome-Scale Structure-Based Clustering](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030065)
    - Murlet
      - 2007, *Bioinformatics*, [Murlet: a practical multiple alignment tool for structural RNA sequences](https://academic.oup.com/bioinformatics/article/23/13/1588/221758)
    - FoldalignM
      - 2007, *Bioinformatics*, [Multiple structural alignment and clustering of RNA sequences](https://academic.oup.com/bioinformatics/article/23/8/926/198434)
  - with pseudoknot
    - 2012, NAR, [CARNA—alignment of RNA structure ensembles](https://academic.oup.com/nar/article/40/W1/W49/1078095), also see their publication in 2010 [A Propagator for Maximum Weight String Alignment with Arbitrary Pairwise Dependencies](https://link.springer.com/content/pdf/10.1007/978-3-642-15396-9_16.pdf)
    - 

- SCFG based
  - Consan
  - Stemloc


## Motif finding
- Random projection
```bash
  g++ src/em-motif.cpp src/utils.cpp src/motifs.cpp  -o bin/em-motif
```
- Expectation Maximization
- Gibbs sampling



## Reading
### Algorithms for General Programming
- <https://web.stanford.edu/class/cs97si/>
### Algorithms in Bioinformatics
- <https://www.comp.nus.edu.sg/~ksung/algo_in_bioinfo/>
- [CS262, computational genomics](https://web.stanford.edu/class/cs262/cgi-bin/index.php)
