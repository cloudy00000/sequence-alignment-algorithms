'''Implementation of the Needleman-Wunsch algorithm for global sequence alignment. 
The class SeqPair contains two methods: _compute_matrices, which computes the score and traceback matrices, and nw, which performs the traceback to obtain the aligned sequences and the final score. 
The user should call the nw method to perform the alignment, providing a substitution matrix and a gap penalty as arguments.'''

class SeqPair:
    def __init__(self, seq1, seq2):
        if not isinstance(seq1, str):
            raise TypeError("seq1 must be a string")
        if not isinstance(seq2, str):
            raise TypeError("seq2 must be a string")
        self.seq1 = seq1
        self.seq2 = seq2
    
    def _compute_matrices(self, sub_matrix, gap_penalty):
        F = [[0 for j in range(len(self.seq1) + 1)] for i in range(len(self.seq2) + 1)]
        T = [['stop' for j in range(len(self.seq1) + 1)] for i in range(len(self.seq2) + 1)]
    
        # Initializing the first row ...
        for j in range(1, len(F[0])):
            F[0][j] = gap_penalty * j
            T[0][j] = 'left'
        
        # ... and the first column of both matrices
        for i in range(1, len(F)):
            F[i][0] = gap_penalty * i
            T[i][0] = 'up'
        
        # Filling the rest of the cells
        for i in range(1, len(F)):
            for j in range(1, len(F[i])):
                up = F[i-1][j] + gap_penalty, 1, 'up'
                left = F[i][j-1] + gap_penalty, 2, 'left'
                diagonal = F[i-1][j-1] + sub_matrix[self.seq1[j-1], self.seq2[i-1]], 3, 'diagonal'
                
                F[i][j], _, T[i][j] = max(up, left, diagonal)
        return F, T

    def nw(self, sub_matrix, gap_penalty):
        '''Perform Needleman-Wunsch global alignment and return the aligned sequences and score.'''
        # First computing the matrices with the hidden method
        F, T = self._compute_matrices(sub_matrix, gap_penalty)
        # The final score is extracted from the final cell of the matrix F
        score = F[-1][-1]
        # Starting the iteration from the final cell of the matrix T
        i = len(T) - 1
        j = len(T[0]) - 1

        
        aln1, aln2 = "", "" #Initializing the aligned sequences as empty strings
        
        while T[i][j] != 'stop':
            if T[i][j] == 'diagonal':
                aln2 = self.seq2[i - 1] + aln2
                aln1 = self.seq1[j - 1] + aln1
                i -= 1
                j -= 1
            elif T[i][j] == 'up':
                aln1 = '-' + aln1
                aln2 = self.seq2[i - 1] + aln2
                i -= 1
            else: # T[i][j] == 'left'
                aln1 = self.seq1[j - 1] + aln1
                aln2 = '-' + aln2
                j -= 1
                
        return aln1, aln2, score
        

