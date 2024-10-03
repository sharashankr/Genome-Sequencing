# Computational Methods for Sequencing Alignment and Genome Assembly

This project focuses on analyzing DNA sequencing data using various algorithms and data structures. Below is a breakdown of the key questions answered by the project code, demonstrating its functionalities.

## Test Dataset

**Test with "phix.fa" and pattern "ATTA":**
- The program identifies the occurrences and leftmost offset of the pattern "ATTA" in the `phix.fa` genomic sequence.
- This serves as a basic functionality test of the `findPattern` class.

## Lambda Virus Genome Analysis (Questions 1-6)

### Q1: AGGT or ACCT Occurrences
- The program calculates how many times the pattern "AGGT" or its reverse complement "ACCT" appears in the lambda virus genome.
- Demonstrates the ability to search for both a pattern and its reverse complement.

### Q2: TTAA Occurrences
- Finds the number of occurrences of the pattern "TTAA" or its reverse complement in the lambda virus genome.
- As the reverse complement is the same as the original pattern, care is taken not to double count.

### Q3: Leftmost Occurrence of ACTAAGT or Reverse Complement
- Identifies the leftmost (smallest) offset where the pattern "ACTAAGT" or its reverse complement appears in the genome.
- Highlights the ability to find the earliest occurrence.

### Q4: Leftmost Occurrence of AGTCGA or Reverse Complement
- Finds the leftmost occurrence of the pattern "AGTCGA" or its reverse complement in the genome.
- Further tests the functionality for detecting the offset of a pattern.

### Q5: TTCAAGCC with up to 2 Mismatches
- Searches for the pattern "TTCAAGCC" in the lambda virus genome, allowing up to 2 mismatches.
- Demonstrates the program's capability of handling inexact matches and tolerating a defined number of mismatches.

### Q6: Leftmost Occurrence of AGGAGGTT with 2 Mismatches
- Identifies the leftmost occurrence of "AGGAGGTT" in the lambda virus genome, allowing up to 2 mismatches.
- Combines both the mismatch feature and offset detection.

### Q7: Sequencing Cycle with the Poorest Quality
- The `checkQuality` class examines a FASTQ sequencing file to determine which sequencing cycle has the most frequent poor-quality scores.
- Identifies the cycle with the lowest-quality scores, helping to pinpoint sequencing issues.
