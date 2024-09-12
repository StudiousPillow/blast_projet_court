import numpy as np
import pandas as pd 
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio import SeqIO
import re

class Sequence(Seq):
    bases = ['A','C','T','G']
    dna_score_matrix = pd.DataFrame(np.array([[5,-4,-4,-4],[-4,5,-4,-4],[-4,-4,5,-4],[-4,-4,-4,5]]), index=["A","C","T","G"], columns = ["A","C","T","G"])
    def __init__(self, seq):
        Seq.__init__(self, seq)
        self.variants = [] # Sequence object
        self.kmers = {} # position : Sequence object
        self.seeds = [] # (queryposition,dbposition)
        self.hits = [] # (queryposition1,queryposition2,dbposition1,dbposition2)
        self.jhits = [] # (queryposition1,queryposition2,dbposition1,dbposition2)
        self.HSP = [] # (queryidx, seqidx, strip_query, strip_seq, score)

    def kmer(self, idx, k):
        """Creates a kmer from a sequence
        
        Args:
            idx : the index to start the kmer
            k : the size of the kmer
        
        Returns : a Sequence object
        """
        return(Sequence(MutableSeq(self[idx:idx+k])))
    
    def calc_score(self, seq):
        """Calculates the score for a sequence (ungapped)
        
        Args:
            seq : the sequence to compare to self
        
        Returns : a score (int)
        """
        score = 0
        for i in range(len(self)):
            score += self.dna_score_matrix.loc[self[i], seq[i]]
        return(score)
    
    def man_get_muted(self, idx, replacement): ## !! doesn't keep the variants and kmers
        """Creates a mutated version of a sequence
        
        Args:
            idx : the index to mutate
            replacement : the base to put
        
        Returns : a Sequence object
        """
        lseq = [letter for letter in str(self)] 
        lseq[idx] = replacement
        seq = ''.join(lseq)
        return(Sequence(seq))
    
    def get_muted(self, idx, replacement): ## !! doesn't keep the variants and kmers
        """Creates a mutated version of a sequence (using MutableSeq)
        
        Args:
            idx : the index to mutate
            replacement : the base to put
        
        Returns : a Sequence object
        """
        inter_seq = MutableSeq(self)
        inter_seq[idx] = replacement
        return(Sequence(inter_seq))
    
    def generate_variants(self, T, refseq, n=0): ## ne pas utiliser sur des grosses séquences
        """Adds variants to self.variants by mutating the sequence
        
        Args:
            T : the score threshold 
            refseq : the original sequence, for recursivity, equals self at first iteration
            n : the index to start mutating, for recursivity
        
        Returns : nothing, fills self.variants, applied to kmers
        """
        for base in self.bases:
            if (n>=len(self)):
                return(None)
            if (self[n]==base):
                continue
            mut = self.get_muted(n, base)
            score = refseq.calc_score(mut)
            if score>=T:
                refseq.variants.append(mut)
                for i in range(n+1,len(mut)):
                    mut.generate_variants(T, refseq, n=i)
                    
    def generate_kmer(self, it, k=6):
        """Adds kmers to self.kmers
        
        Args:
            it : the space between the kmers
            k : the size of the kmers
        
        Returns : nothing, fills self.kmers, applied to the query sequence
        """
        self.kmers = {}
        for idx in range(0,len(self),it):
            self.kmers[idx] = self.kmer(idx, k)

    def find_close_seeds(self, margin=10):
        """finds the seeds close to each other and puts them in self.hits
        
        Args:
            margin : how close the sequence can be
        
        Returns : nothing, fills self.hits from self.seeds
        """
        for x in range(len(self.seeds)):
            for y in range(x+1,len(self.seeds)):
                lenseed = self.seeds[x][1] - self.seeds[x][0]
                I=self.seeds[x][0] #position in query of first seed
                J=self.seeds[y][0] #position in query of second seed
                Lquery= np.abs(I-J)-lenseed
                i=self.seeds[x][1] #position in dbseq of first seed
                j=self.seeds[y][1] #position in dbseq of second seed
                Lseq = np.abs(j-i)-lenseed
                if np.abs(Lseq-Lquery)<=margin:
                    self.hits.append((self.seeds[x][0], self.seeds[y][0], self.seeds[x][1], self.seeds[y][1]))

    def alignment(self, query):
        """Aligns locally (Waterman) two sequences in one direction
        
        Args:
            query : the sequence to align to self
        
        Returns :  
            aln_seq : the alignement of the self
            aln_query : the alignement of the query
            score : the score calculated
        """
        gap_add_cost = -2
        gap_open_cost = -8
        score_mat = np.zeros((len(query) + 1, len(self) + 1))
        path_mat = np.zeros((len(query) + 1, len(self) + 1), dtype=int)
        right, down, diag = 1, 2, 3
        max_score = 0
        max_pos = (0,0)

        for i in range(1, len(query) + 1):
            for j in range(1, len(self) + 1):
                sub = score_mat[i-1, j-1] + self.dna_score_matrix.loc[str(query[i-1]),str(self[j-1])]
                deletion = score_mat[i,j-1] + (gap_add_cost if path_mat[i,j-1]==3 else gap_open_cost)
                addition = score_mat[i-1,j] + (gap_add_cost if path_mat[i-1,j]==2 else gap_open_cost)

                score_mat[i][j] = max(sub, deletion, addition)

                if score_mat[i][j] == sub:
                    path_mat[i][j] = diag
                elif score_mat[i][j] == deletion:
                    path_mat[i][j] = down
                elif score_mat[i][j] == addition:
                    path_mat[i][j] = right

                if score_mat[i][j] >= max_score:
                    max_score = score_mat[i][j]
                    max_pos = (i, j)

        aln_query = []
        aln_seq = []
        score = score_mat[i][j]
        i,j=max_pos
        while (i != 0) and (j != 0):
            if path_mat[i][j] == diag:
                aln_query.append(query[i - 1])
                aln_seq.append(self[j - 1])
                i -= 1
                j -= 1
            elif path_mat[i][j] == down:
                aln_query.append(query[i - 1])
                aln_seq.append('-')
                i -= 1
            elif path_mat[i][j] == right:
                aln_query.append('-')
                aln_seq.append(self[j - 1])
                j -= 1
        aln_query = aln_query[::-1]
        aln_seq = aln_seq[::-1]
        return(aln_seq, aln_query, score)
    
    def join_hits(self, query, pos_hit, k):
        """joins two hits
        
        Args:
            query : the query sequence
            pos_hit : the positions of the hits
                - query_hit1, query_hit2, target(self)_hit1, target(self)_hit2 
        
        Returns :  
            aln_seq : the alignement on the query
            aln_query : the alignement on the target
            x : the start index of the joined hit on query
            X : the start index of the joined hit on target 
            score : the score calculated
        """
        x, y, X, Y = pos_hit
        length = Y-X+k
        score = sum([self.dna_score_matrix[query[x+i], self[X+i]] for i in range(0, length)])
        self.jhits.append((query[x:x+length], self[X:X+length], x, X, score))
        return(query[x:x+length], self[X:X+length], x, X, score)
    
    def extend_joined_hits(self, query, joined_hit, k=6):
        """extend a joined hit in two directions
        
        Args:
            query : the query sequence
            joined_hit : the output of join_hits, found in self.jhits
        Returns :  
            index of hit in query
            index of hit in target
            aligned query sequence
            aligned target sequence
            score

            adds these as a tuple in self.HSP
        """
        queryhit, seqhit, queryidx, seqidx, joinscore = joined_hit
        # reverse (left)
        rev_query = query[max(queryidx-1, 0)::-1]
        rev_seq = self[max(seqidx-1, 0)::-1]
        raln_query, raln_seq, rscore = rev_query.alignment(rev_seq)
        # normal (right)
        mod_query = query[queryidx+len(queryhit):]
        mod_seq = self[seqidx+len(seqhit):]
        aln_query, aln_seq, score = mod_query.alignment(mod_seq)
        # combine
        f_query = rev_query[::-1] + queryhit + mod_query
        f_seq = rev_seq[::-1] + seqhit + mod_seq
        f_score = score + rscore + joinscore
        strip_query = raln_query.replace('_', '')
        strip_seq = raln_query.replace('_', '')
        queryidx = queryidx - len(strip_query)
        seqidx = seqidx - len(strip_seq)
        self.HSP.append((queryidx, seqidx, strip_query, strip_seq, score))
        return queryidx, seqidx, f_query, f_seq, score

def read_fasta(filename):
    """reads a fasta file and outputs a Sequence object
    
    Args:
        filename : the path to the fasta file
    Returns :  
        Sequence object
    """
    records = SeqIO.parse(filename, "fasta")
    return([Sequence(rec.seq) for rec in records])
    
if __name__ == "__main__":
    K = 6
    database = read_fasta("data/Sars-CoV-2_sequences.fasta")
    query_sequence = database[0]
    database = database[1:2]
    ## cut the query sequence into smaller fragments, kmers
    query_sequence.generate_kmer(4,K)
    # for each kmers, create close sequences with an alignement score less than T, stored in each kmer in .variants attribute
    for kmer in query_sequence.kmers.values():
        kmer.generate_variants(30,kmer)
        kmer.variants.append(kmer)
    ## for each kmers, for each variant, in each database sequence, look if there is a match, store it in each database sequence .seeds
    print('search')
    for kmeridx in range(len(query_sequence.kmers)):
        print(kmeridx, '/', len(query_sequence.kmers), end='\r')
        for variant in list(query_sequence.kmers.values())[kmeridx].variants:
            for db_seq in database:
                matches = re.finditer(str(variant), str(db_seq))
                for match_idx in matches:
                    queryidx = list(query_sequence.kmers.keys())[kmeridx]
                    db_seq.seeds.append((queryidx, match_idx.span()[0]))
    ## filter and keep only the seeds at the same distance on the query and the db sequence, put it in .hits
    print('filter seeds -> hits')
    for db_seq in database:
        db_seq.find_close_seeds()
        print('hit',db_seq.hits)
    ## join the hits, in jhits
    for db_seq in database:
        for hit in db_seq.hits:
            db_seq.join_hits(query_sequence, hit, K)
    ## extend the jhits, create the HSP
    for db_seq in database:
        for jhit in db_seq.jhits:
            db_seq.extend_joined_hits(jhit, k = K)
    ## print the result
    for db_seq in database:
        for HSP in db_seq.HSP:
            print(f'query {HSP[2]}: {HSP[0]}\n target {HSP[3]}: {HSP[1]}\n {HSP[4]}')