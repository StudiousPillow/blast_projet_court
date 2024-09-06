import numpy as np
import pandas as pd ## alternative ?

class Sequence:
    dna_score_matrix = pd.DataFrame(np.array([[5,-4,-4,-4],[-4,5,-4,-4],[-4,-4,5,-4],[-4,-4,-4,5]]), index=["A","C","T","G"], columns = ["A","C","T","G"])
    def __init__(self, seq):
        self.seq = seq
    def generate_kmer(self, idx, k):
        return(self.seq[idx:idx+k])
    def generate_kuplet(self, score_lim):
        pass
    # def calc_score(self, ):

S = Sequence("ACT")
print(S.dna_score_matrix)