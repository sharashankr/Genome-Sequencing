# -*- coding: utf-8 -*-
"""
This project includes findPatternV2 class. The findPatternV2
class finds the occurrence and position of a given pattern in a given genomic 
sequence.
@author: Sharath
"""
import bisect

class findPatternV2:
    """
    This class finds the occurrence and position of a given pattern in a given 
    genomic sequence in a file.
    """
    def __init__(self, pattern, filename=False):
        # initiate parameters
        self.pattern = pattern
        self.filename = filename
        
    def readGenome(self):
        """
        Read genomic DNA sequence to a string
        """
        genome = ""
        with open(self.filename, "r") as f:
            for line in f:
                # skip header
                if line[0] != ">":
                    genome += line.rstrip()
        return genome
        
    def naiveMatch(self, numberOfMismatch, text=None):
        """
        This is naive match to find the index of matched patterns in a genome
        and calculate number of total character comparisons and alignments
        """
        genome = text if text else self.readGenome()
        pattern = self.pattern
        occurrences = []
        alignments = 0
        comparisons = 0

        for i in range(len(genome) - len(pattern) + 1):
            match = True
            counter = 0
            for j in range(len(pattern)):
                comparisons += 1
                if pattern[j] != genome[i+j]:
                    counter += 1
                if counter > numberOfMismatch:
                    match = False
                    break
            if match:
                occurrences.append(i)
            alignments += 1
            
        return occurrences, alignments, comparisons
        
    def boyerMoore(self, numberOfMismatch, bm, text=None):
        """
        This is the Boyer-Moore match to find the index of matched patterns in a genome
        and calculate number of total character comparisons and alignments
        """
        i = 0
        genome = text if text else self.readGenome()
        pattern = self.pattern
        occurrences = []
        alignments = 0
        comparisons = 0
        
        while i < len(genome) - len(pattern) + 1:
            shift = 1
            match = True
            
            for j in range(len(pattern) - 1, -1, -1):
                comparisons += 1
                if pattern[j] != genome[i+j]:
                    badCharacterSkip = bm.bad_character_rule(j, genome[i+j])
                    goodSuffixSkip = bm.good_suffix_rule(j)
                    shift = max(shift, badCharacterSkip, goodSuffixSkip)
                    match = False
                    break
                    
            if match:
                occurrences.append(i)
                goodSuffixSkip = bm.match_skip()
                shift = max(shift, goodSuffixSkip)
            i += shift
            alignments += 1
            
        return occurrences, alignments, comparisons
        
    def matchedIndex(self, index, k_mer, pattern, isSubseqIndex):
        """
        Find number of hits, occurrences and time of occurrence for a given pattern
        using string index in genome
        """
        genome = self.readGenome()
        occurrences_match = []
        hit_index = []
        occurrence_genome = []
        counter = 0
        
        length = len(pattern) - k_mer + 1 if not isSubseqIndex else isSubseqIndex
        
        for i in range(length):  # loop over to generate kmers
            pattern_q = pattern[i:i+k_mer] if not isSubseqIndex else pattern[i:]
            hits = index.query(pattern_q)  # query each kmer
            
            for hit in hits:
                counter += 1  # count total number of hits
                text = genome[hit - i: hit + len(pattern) - i]
                
                if hit - i not in hit_index:  # avoid duplicated counts
                    hit_index.append(hit - i)
                    occurrence, _, _ = self.naiveMatch(2, text)
                    occurrences_match.extend(occurrence)
                    
                if len(occurrence) != 0 and hit - i not in occurrence_genome:
                    occurrence_genome.append(hit - i)
                    
        return occurrence_genome, len(occurrences_match), counter
        
class Index:
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
            
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
            
        return hits
    
    def genome_index(self):
        return self.index       

class SubseqIndex:
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
            
        return hits

if __name__ == "__main__":
    from bm_preproc import BoyerMoore
    # Questions 1-3
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    
    # Q1: How many alignments does the naive exact matching algorithm try when 
    # matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG 
    # (derived from human Alu sequences) to the excerpt of human chromosome 1? 
    # (Don't consider reverse complements.)
    pattern = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
    patterns = findPatternV2(pattern, filename)
    print("Q1: The alignments for naive match algorithm is {}\n".format(patterns.naiveMatch(0)[1]))
    patterns = findPatternV2("GGCGCGGTGGCTCACGCCTGTAAT", filename)    
    
    # Q2: How many character comparisons does the naive exact matching algorithm 
    # try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG 
    # (derived from human Alu sequences) to the excerpt of human chromosome 1? 
    # (Don't consider reverse complements.)
    print("Q2: The characters comparisons for naive match algorithm is {}\n".format(patterns.naiveMatch(0)[2]))
    
    # Q3: How many alignments does Boyer-Moore try when matching the string 
    # GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG 
    # (derived from human Alu sequences) to the excerpt of human chromosome 1? 
    # (Don't consider reverse complements.)
    print("Q3: The alignments for Boyer-Moore algorithm is {}\n".format(patterns.boyerMoore(0, 
          BoyerMoore(pattern, "ACGT"))[1]))
    
    # Q4: How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, 
    # which is derived from a human Alu sequence, occur with up to 2 
    # substitutions in the excerpt of human chromosome 1? 
    # (Don't consider reverse complements here.)
    k_mer = 8      
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    genome = patterns.readGenome()
    index = Index(genome, k_mer)
    occurrences, numberOfOccurs, numberOfhits = patterns.matchedIndex(index, 
                                               k_mer, pattern, isSubseqIndex=False)
    print("Q4: Within 2 mismatches, the string occurs {} times\n".format(numberOfOccurs))
    
    # Q5: Using the instructions given in Question 4, how many total index hits 
    # are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with 
    # up to 2 substitutions in the excerpt of human chromosome 1?
    print("Q5: Within 2 mismatches, the total index hits are {}\n".format(numberOfhits))
    
    # Q6: When using this function, how many total index hits are there when 
    # searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the 
    # excerpt of human chromosome 1? (Again, don't consider reverse complements.)
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    k_mer = 8
    vial = 3
    index = SubseqIndex(genome, k_mer, vial)
    occurrences, numberOfOccurs, numberOfhits = patterns.matchedIndex(index, 
                                               k_mer, pattern, isSubseqIndex=vial)

    print("Q6: Within 2 mismatches, the hits are {}\n".format(numberOfhits))
