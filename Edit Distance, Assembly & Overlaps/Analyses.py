# -*- coding: utf-8 -*-
"""
This project includes the findPatternV3 class. The findPatternV3
class finds the edit distance of a given pattern in a given genomic 
sequence and constructs overlap graphs.

Author: Sharath
"""

from itertools import permutations

class findPatternV3:
    """
    This class finds the edit distance of a given pattern in a given genomic 
    sequence and constructs overlap graphs.
    """

    def __init__(self, filename, pattern=None):
        # Initialize parameters
        self.pattern = pattern
        self.filename = filename
        
    def readGenome(self):
        """
        Read genomic DNA sequence from a file and return it as a string.
        """
        genome = ""
        with open(self.filename, "r") as f:
            for line in f:
                # Skip header
                if not line.startswith(">"):
                    genome += line.strip()
        return genome
        
    def readFastq(self):
        """
        Read DNA sequences and quality scores from a FASTQ sequencing file 
        and return them as lists.
        """
        with open(self.filename, "r") as f:
            sequences = []
            qualities = []
            while True:
                f.readline()  # Skip name line
                seq = f.readline().strip()  # Read sequence line
                f.readline()  # Skip strand line
                qual = f.readline().strip()  # Read quality line
                if not seq:  # Finish read
                    break
                # Add sequence and quality information to list
                sequences.append(seq)
                qualities.append(qual)
        return sequences, qualities
        
    def editDistance(self):
        """
        Implement dynamic programming algorithm to calculate the edit distance 
        between a given pattern and a given genome.
        """
        pattern = self.pattern
        genome = self.readGenome()
        pattern_length = len(pattern) + 1
        genome_length = len(genome) + 1
        
        # Generate matrix
        matrix = [[0] * genome_length for _ in range(pattern_length)]
        
        # Initialize the first column
        for i in range(pattern_length):
            matrix[i][0] = i
            
        for i in range(1, pattern_length):
            for j in range(1, genome_length):
                dist_hor = matrix[i][j - 1] + 1
                dist_vel = matrix[i - 1][j] + 1
                dist_diag = matrix[i - 1][j - 1] + 1 if pattern[i - 1] != genome[j - 1] else matrix[i - 1][j - 1]
                matrix[i][j] = min(dist_hor, dist_vel, dist_diag)
                
        return min(matrix[-1])                
            
    def phraseReads(self, k_mer):
        """
        Construct a dictionary of prefixes and suffixes of reads with reads as
        keys and prefixes/suffixes as values.
        """
        reads, _ = self.readFastq()
        reads_dict = {}
        for read in reads:
            for i in range(len(read) - k_mer + 1):
                substring = read[i:i + k_mer]
                reads_dict.setdefault(substring, set()).add(read)
        return reads_dict

    def overlap(self, read1, read2, k_mer):
        """
        Find the leftmost overlap offset between two reads.
        """
        start = 0
        while True:
            start = read1.find(read2[:k_mer], start)
            if start == -1:
                return 0  # No overlap
            if read2.startswith(read1[start:]):
                return len(read1) - start
            start += 1
            
    def overlapGraph(self, k_mer):
        """
        Construct a graph with reads as nodes and values as other reads that 
        overlap with the given read.
        """
        reads_dict = self.phraseReads(k_mer)
        reads, _ = self.readFastq()
        graph = {}
        
        for read1 in reads:
            k_mer_string = read1[len(read1) - k_mer:]
            if k_mer_string in reads_dict:
                edges = set()
                reads_set = reads_dict[k_mer_string]
                
                for read2 in reads_set:
                    if read1 != read2:  # Skip self-comparison
                        offset = self.overlap(read1, read2, k_mer)           
                        if offset > 0:  # Skip non-overlapped pairs
                            edges.add(read2)  # Add overlapped reads to be values
                            graph[read1] = edges           
        return graph
                        
    def naive_overlap_map(self, k_mer):
        """
        Construct a graph with keys as pairs of reads with overlap and 
        values as the leftmost offset of the overlap.
        """
        graph = {}
        reads, _ = self.readFastq()
        
        for read1, read2 in permutations(reads, 2):
            # Skip non-overlapped reads
            if read1[len(read1) - k_mer:] in read2:
                offset = self.overlap(read1, read2, k_mer)
                # Check if reads[i] overlapped with reads[j]
                if offset != 0:
                    graph[(read1, read2)] = offset  
        return graph

if __name__ == "__main__":
    
    # Q1: What is the edit distance of the best match between pattern 
    # GCTGATCGATCGTACG and the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    pattern = "GCTGATCGATCGTACG"
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    patterns = findPatternV3(filename, pattern)
    edit_dist = patterns.editDistance()
    print(f"Q1: The edit distance of the best match between the pattern and the genome is {edit_dist}\n")
           
    # Q2: What is the edit distance of the best match between pattern 
    # GATTTACCAGATTGAG and the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    pattern = "GATTTACCAGATTGAG"
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    patterns = findPatternV3(filename, pattern)
    edit_dist = patterns.editDistance()
    print(f"Q2: The edit distance of the best match between the pattern and the genome is {edit_dist}\n")
    
    # Q3: Picture the overlap graph corresponding to the overlaps just calculated. 
    # How many edges are in the graph? In other words, how many distinct pairs 
    # of reads overlap?
    # Q4: Picture the overlap graph corresponding to the overlaps computed for 
    # the previous question. How many nodes in this graph have at least one 
    # outgoing edge? (In other words, how many reads have a suffix involved in 
    # an overlap?)
    
    import time
    
    t1 = time.time()
    filename = "../data/ERR266411_1.for_asm.fastq"
    patterns = findPatternV3(filename)
    k_mer = 30
    graph = patterns.naive_overlap_map(k_mer)
    t2 = time.time()
    
    print(f"Running time for naive overlap mapping: {t2 - t1:.2f} sec\n")
    
    reads = patterns.phraseReads(k_mer)
    t3 = time.time()
    
    print(f"Running time for phrase reads: {t3 - t2:.2f} sec\n")
    
    graph = patterns.overlapGraph(k_mer)
    t4 = time.time()
    
    print(f"Running time for optimized algorithm: {t4 - t3:.2f} sec\n")
    
    numberOfNodes = len(graph)
    numberOfEdges = sum(len(edges) for edges in graph.values())
    
    print(f"Q3: The total number of edges is {numberOfEdges}\n")
    print(f"Q4: The total number of nodes is {numberOfNodes}")

