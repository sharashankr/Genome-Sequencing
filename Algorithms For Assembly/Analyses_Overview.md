# Analysis of Questions on Shortest Common Superstrings and Genome Assembly

## Problem 1: Shortest Common Superstring Length  
The first problem investigates the concept of the shortest common superstring (SCS) for a set of strings. A common superstring is a string that contains each of the input strings as a substring. Multiple shortest common superstrings can exist for the same set of input strings.

### Input Strings:
- CCT
- CTT
- TGC
- TGG
- GAT
- ATT

### Result:  
The length of the shortest common superstring for the above strings was calculated, confirming the total length.

---

## Problem 2: Number of Different Shortest Common Superstrings  
The second question examines how many different shortest common superstrings can be generated from the same input strings used in the previous question.

### Result:  
The number of distinct shortest common superstrings was calculated, providing insight into the combinatorial nature of superstring assembly.

---

## Problem 3: Genome Assembly  
The third problem focuses on assembling a viral genome from a given set of reads, which are exact copies of substrings from the virus genome. Each read is the same length (100 bases), and it is specified that no sequencing errors or variations from the reverse strand need to be considered.

### Assembly Method:  
A greedy approach for finding the shortest common superstring was applied to assemble the reads. Given that the virus genome is exactly 15,894 bases long, the assembly process was optimized for efficiency, leveraging previously discussed methods.

### Result:  
The assembled genome was saved to a file, and the number of adenine (A) nucleotides in the final assembled genome was counted and reported.

---

## Problem 4: Counting Thymine Nucleotides  
The fourth question continues the analysis of the assembled viral genome by counting the number of thymine (T) nucleotides.

### Result:  
The total number of T nucleotides in the assembled genome was confirmed.

---

## Summary  
These problems explore the foundational concepts of sequence assembly and superstring theory within bioinformatics. The first two questions assess the properties of common superstrings, while the last two questions focus on practical applications in genomic assembly, including the efficient handling of sequencing reads and nucleotide counting. Each step contributes to understanding the complexity and significance of genomic data analysis in biological research.
