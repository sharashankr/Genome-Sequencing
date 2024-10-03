# -*- coding: utf-8 -*-
"""
This project includes two classes: findPattern and checkQuality. 
The findPattern class identifies the occurrence and position of a given pattern 
in a genomic sequence. The checkQuality class examines sequencing quality for each cycle.

@author: Sharath
"""

class findPattern:
    """
    This class finds the occurrence and position of a given pattern in a genomic sequence.
    """
    def __init__(self, pattern, filename):
        # Initialize pattern and filename
        self.pattern = pattern
        self.filename = filename

    def readGenome(self):
        """
        Reads the genomic DNA sequence from a file into a string.
        """
        genome = ""
        with open(self.filename, "r") as f:
            for line in f:
                # Skip the header
                if not line.startswith(">"):
                    genome += line.rstrip()
        return genome

    def reverseComplement(self):
        """
        Generates the reverse complement of the given DNA sequence.
        """
        complement = {"A": "T", "C": "G", "T": "A", "G": "C"}
        revComPattern = ""
        for nt in self.pattern:
            revComPattern = complement[nt] + revComPattern  # Build reverse complement
        return revComPattern

    def match(self, string1, string2, numOfMismatch):
        """
        Returns True if string1 and string2 match under the allowed number of mismatches.
        """
        if len(string1) != len(string2):
            return False
        mismatches = sum(1 for i in range(len(string1)) if string1[i] != string2[i])
        return mismatches <= numOfMismatch

    def patternIdentifier(self, numOfMismatch):
        """
        Finds positions of the given pattern and its reverse complement in the genome.
        """
        genome = self.readGenome()
        revComPattern = self.reverseComplement()
        patternLength = len(self.pattern)
        occurrences = []

        for i in range(patternLength):
            for j in range(i, len(genome), patternLength):
                genomeMotif = genome[j: j + patternLength]
                # Compare genomic motif with the pattern and reverse complement
                if (self.match(genomeMotif, self.pattern, numOfMismatch) or
                    self.match(genomeMotif, revComPattern, numOfMismatch)) and j not in occurrences:
                    occurrences.append(j)
        return occurrences


class checkQuality:
    """
    This class examines the quality of sequencing for each cycle.
    """
    def __init__(self, filename):
        self.filename = filename

    def readFastq(self):
        """
        Reads DNA sequences and quality scores from a FASTQ file into lists.
        """
        sequences, qualities = [], []
        with open(self.filename, "r") as f:
            while True:
                f.readline()  # Skip name line
                seq = f.readline().rstrip()  # Read sequence
                f.readline()  # Skip strand line
                qual = f.readline().rstrip()  # Read quality scores
                if len(seq) == 0:
                    break
                sequences.append(seq)
                qualities.append(qual)
        return sequences, qualities

    def phred33ToQ(self, qualString):
        """
        Converts a quality string to a list of Phred+33 quality scores.
        """
        return [ord(qual) - 33 for qual in qualString]

    def findPoorQuality(self):
        """
        Finds the index of the poorest quality score in each sequence.
        """
        _, qualities = self.readFastq()
        lowestQScoreIndex = []
        for qualString in qualities:
            qScore = self.phred33ToQ(qualString)
            lowestQScoreIndex.append(qScore.index(min(qScore)))
        return lowestQScoreIndex

    def countPoorQuality(self):
        """
        Counts the number of poorest quality scores for each cycle.
        """
        import collections
        return collections.Counter(self.findPoorQuality())

    def plotHist(self):
        """
        Plots a histogram showing the distribution of the poorest quality scores.
        """
        import matplotlib.pyplot as plt
        data = self.countPoorQuality()
        plt.bar(data.keys(), data.values())
        plt.show()


if __name__ == "__main__":
    # Test Dataset
    filename = "../data/phix.fa"
    pattern = "ATTA"
    patterns = findPattern(pattern, filename)
    print("Test dataset results - Occurrences and leftmost offset:")
    print(len(patterns.patternIdentifier(0)), min(patterns.patternIdentifier(0)), "\n")

    # Question 1-6 (Lambda Virus Genome Analysis)
    filename1 = "../data/lambda_virus.fa"

    # Q1: How many times does AGGT or its reverse complement (ACCT) occur in the lambda virus genome?
    pattern = "AGGT"
    patterns = findPattern(pattern, filename1)
    print(f"Q1: The 'AGGT' or 'ACCT' occurs {len(patterns.patternIdentifier(0))} times \n")

    # Q2: How many times does TTAA or its reverse complement occur in the lambda virus genome?
    pattern = "TTAA"
    patterns = findPattern(pattern, filename1)
    print(f"Q2: The 'TTAA' occurs {len(patterns.patternIdentifier(0))} times \n")

    # Q3: What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement?
    pattern = "ACTAAGT"
    patterns = findPattern(pattern, filename1)
    print(f"Q3: The offset of the leftmost occurrence of ACTAAGT is {min(patterns.patternIdentifier(0))} \n")

    # Q4: What is the offset of the leftmost occurrence of AGTCGA or its reverse complement?
    pattern = "AGTCGA"
    patterns = findPattern(pattern, filename1)
    print(f"Q4: The offset of the leftmost occurrence of AGTCGA is {min(patterns.patternIdentifier(0))} \n")

    # Q5: How many times does TTCAAGCC occur in the Lambda virus genome allowing up to 2 mismatches?
    pattern = "TTCAAGCC"
    patterns = findPattern(pattern, filename1)
    print(f"Q5: The 'TTCAAGCC' occurs {len(patterns.patternIdentifier(2))} times with up to 2 mismatches \n")

    # Q6: What is the offset of the leftmost occurrence of AGGAGGTT allowing up to 2 mismatches?
    pattern = "AGGAGGTT"
    patterns = findPattern(pattern, filename1)
    print(f"Q6: The offset of the leftmost occurrence of AGGAGGTT with up to 2 mismatches is {min(patterns.patternIdentifier(2))} \n")

    # Q7: Identify which sequencing cycle has the most poor-quality scores.
    filename2 = "../data/ERR037900_1.first1000.fastq"
    qualities = checkQuality(filename2)
    counters = qualities.countPoorQuality()
    print(f"Q7: The cycle with the most frequent poor quality is {max(counters, key=lambda x: counters[x])}")
