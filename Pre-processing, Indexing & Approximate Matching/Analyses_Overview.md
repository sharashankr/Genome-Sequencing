# Analysis of String Matching Algorithms in Genomic Sequences

This document provides a brief note on each of the questions analyzed in the code, summarizing the methods used and the significance of the findings.

## Question 1
**How many alignments does the naive exact matching algorithm try when matching the string `GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG` to the excerpt of human chromosome 1? (Don't consider reverse complements.)**

The naive exact matching algorithm checks for all possible starting positions in the genome where the pattern could match. It counts every attempt to align the pattern with the genome sequence, regardless of whether the alignment is successful or not. The result provides insight into the algorithm's efficiency, highlighting how many times the algorithm had to check for a match, which could be computationally expensive, especially for longer patterns and larger sequences.

## Question 2
**How many character comparisons does the naive exact matching algorithm try when matching the string `GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG` to the excerpt of human chromosome 1? (Don't consider reverse complements.)**

This question focuses on the number of character comparisons made by the naive algorithm. Each character comparison is a check between the corresponding characters in the pattern and the genome at a specific alignment position. This number indicates the algorithm's performance, as a high count of character comparisons suggests inefficiencies in the algorithm, especially when there are few matches, leading to many unnecessary checks.

## Question 3
**How many alignments does the Boyer-Moore algorithm try when matching the string `GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG` to the excerpt of human chromosome 1? (Don't consider reverse complements.)**

The Boyer-Moore algorithm is generally more efficient than the naive algorithm due to its use of heuristics that skip sections of the text based on mismatches. This question examines how many alignments the Boyer-Moore algorithm attempted to find matches for the same pattern in the genome. A lower number of alignments compared to the naive algorithm indicates the effectiveness of Boyer-Moore in reducing the search space, which is particularly beneficial for large datasets.

## Question 4
**How many times does the string `GGCGCGGTGGCTCACGCCTGTAAT`, which is derived from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human chromosome 1? (Don't consider reverse complements here.)**

This question investigates the occurrence of the specified pattern in the genome, allowing for a maximum of 2 substitutions (mismatches). This is crucial in genomics, where variations can occur in sequences. The count of occurrences provides information on the pattern's frequency, helping understand its biological significance and potential functional roles in the genome, particularly for sequences related to Alu elements, which are known for their variability.

## Question 5
**Using the instructions given in Question 4, how many total index hits are there when searching for occurrences of `GGCGCGGTGGCTCACGCCTGTAAT` with up to 2 substitutions in the excerpt of human chromosome 1?**

This question is similar to Question 4 but focuses specifically on the total index hits, which may include duplicate occurrences across the search. It emphasizes the efficiency of using an index for searching genomic sequences, demonstrating the benefit of data structures that improve search performance by limiting redundant checks. This metric is significant for assessing the efficiency of genomic searches in broader applications, such as variant detection.

## Question 6
**When using this function, how many total index hits are there when searching for `GGCGCGGTGGCTCACGCCTGTAAT` with up to 2 substitutions in the excerpt of human chromosome 1? (Again, don't consider reverse complements.)**

This question applies a subsequence indexing approach to search for the same pattern with a focus on how many total hits are identified. The use of a subsequence index allows for faster searching through the genome, demonstrating the utility of indexing techniques in bioinformatics. The results from this search can inform researchers about the potential variants or mutations present in the genome and their biological implications, particularly in understanding the genomic architecture related to Alu sequences.

## Conclusion
These questions and their corresponding analyses highlight different aspects of string matching algorithms in the context of genomic sequence analysis. They reveal insights into algorithm efficiency, the importance of substitutions in biological sequences, and the utility of indexing in improving search performance, all critical elements in computational biology and genomics.
