'''This file contains the implementation of the Smith-Waterman algorithm for local sequence alignment. 
The class SeqPair inherits from the class defined in global_aln.py, and contains two methods: _compute_local_matrices, 
which computes the score and traceback matrices for local alignment, and sw, which performs the traceback to obtain the aligned sequences and the final score.'''

import global_aln as nw

# Extending the definition of the class SeqPair using inheritance
class SeqPair(nw.SeqPair):

    def _compute_local_matrices(self, sub_matrix, gap_penalty):
        # Initialize the matrices
        F, T = [], []
        for i in range(len(self.seq2) + 1):
            row_f = [0] * (len(self.seq1) + 1)
            F.append(row_f)
            row_t = ['stop'] * (len(self.seq1) + 1)
            T.append(row_t)
        
        # Filling the rest in the matrices
        for i in range(1, len(F)):
            for j in range(1, len(F[0])):
                
                diagonal = F[i-1][j-1] + sub_matrix[self.seq1[j-1], self.seq2[i-1]], 0, 'diagonal'
                up = F[i-1][j] + gap_penalty, 0, 'up'
                left = F[i][j-1] + gap_penalty, 0, 'left'
                zero = 0, 1, 'stop' 
                F[i][j], _, T[i][j] = max(diagonal, left, up, zero)
        return F, T
    
    def sw(self, sub_matrix, gap_penalty):
        F, T = self._compute_local_matrices(sub_matrix, gap_penalty)
        
        # Trace back from the best-scoring cell
        score = F[1][1]
        start_i, start_j = 1, 1
        for i in range(1, len(self.seq2) + 1):
            for j in range(1, len(self.seq1) + 1):
                if F[i][j] > score:
                    score = F[i][j]
                    start_i, start_j = i, j
        i = start_i
        j = start_j
        
        # This part of the algorithm is identical to the NW
        aln1, aln2 = "", ""
        while not(T[i][j] == 'stop'): # Only correct condition for stopping
            if T[i][j] == 'diagonal':
                aln1 = self.seq1[j-1] + aln1
                aln2 = self.seq2[i-1] + aln2
                i -= 1
                j -= 1
            elif T[i][j] == 'up':
                aln1 = '-' + aln1
                aln2 = self.seq2[i-1] + aln2
                i -= 1
            else:
                aln1 = self.seq1[j-1] + aln1
                aln2 = '-' + aln2
                j -= 1
        
        return aln1, aln2, score
        

