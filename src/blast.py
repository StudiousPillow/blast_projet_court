import numpy as np
import pandas as pd ## alternative ?

class Sequence:
    bases = ['A','C','T','G']
    dna_score_matrix = pd.DataFrame(np.array([[5,-4,-4,-4],[-4,5,-4,-4],[-4,-4,5,-4],[-4,-4,-4,5]]), index=["A","C","T","G"], columns = ["A","C","T","G"])
    def __init__(self, seq):
        self.seq = seq
        self.variants = {} # Sequence object : score
        self.kmers = {} # position : Sequence object
    def __str__(self):
        return(f'Seq({self.seq})')
    def __repr__(self):
        return(f'Seq({self.seq})')
    def generate_kmer(self, idx, k):
        return(self.seq[idx:idx+k])
    def calc_score(self, seq):
        score = 0
        for i in range(len(self.seq)):
            score += self.dna_score_matrix.loc[self.seq[i], seq.seq[i]]
        return(score)
    def get_muted(self, idx, replacement):
        lseq = [letter for letter in self.seq] 
        lseq[idx] = replacement
        seq = ''.join(lseq)
        return(Sequence(seq))
    def generate_kuplets(self, T, refseq, n=0):
        for base in self.bases:
            print(base, ":::::::::::::")
            mut = self.get_muted(n, base)
            score = refseq.calc_score(mut)
            if score>=T:
                print(mut)
                refseq.variants[mut]=float(score)
                if n==len(mut.seq):
                    return()
                for i in range(n+1,len(mut.seq)):
                    mut.generate_kuplets(T, refseq, n=i)
                    
            

    
if __name__ == "__main__":
    S1 = Sequence("ACTT")
    S2 = Sequence("ACTG")
    print(S1.calc_score(S2))
    print(S1)
    S1.generate_kuplets(20, S1)
    print(S1.variants)

