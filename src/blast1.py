import numpy as np
import pandas as pd ## alternative ?

class Sequence:
    bases = ['A','C','T','G']
    dna_score_matrix = pd.DataFrame(np.array([[5,-4,-4,-4],[-4,5,-4,-4],[-4,-4,5,-4],[-4,-4,-4,5]]), index=["A","C","T","G"], columns = ["A","C","T","G"])
    def __init__(self, seq):
        self.seq = seq
    def generate_kmer(self, idx, k):
        return(self.seq[idx:idx+k])
    def generate_kuplet(self, score_lim):
        pass
    def calc_score(self, seq):
        score = 0
        for i in range(len(self.seq)):
            score += self.dna_score_matrix.loc[self.seq[i], seq.seq[i]]
        return(score)
    def get_muted(self, idx, replacement):
        lseq = self.seq.split('') 
        lseq[idx] = replacement
        seq = ''.join(lseq)
        return(Sequence(seq))
        
    def generate_kuplet(self, T):
        kuplets = []

            

    

S1 = Sequence("ACTTTGCCGC")
S2 = Sequence("AGGGGGCCGC")
print(S1.calc_score(S2))

