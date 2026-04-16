class AlnSeq:
    def __init__(self, seq1, seq2):
        if (type(seq1) != str):
            raise TypeError("seq1 must be a string")
        if (type(seq2) != str):
            raise TypeError("seq2 must be a string")
        if len(seq1) != len(seq2):
            raise ValueError("strings of different length")
        #we check all possible errors before doing the initialization
        # We want to initialize all possible attributes in the __init__ for readability
        self.seq1 = seq1
        self.seq2 = seq2
        self.distance = None
        self.similarity = None
    
    def compute_hamming_distance(self):
        ''' Function to compute the Hamming distance between the two sequences. '''
        #The distance is saved as an attribute of the class, to avoid computing it twice if we call both this function and the compute_hamming_similarity function. 
        self.distance = 0
        for i in range(len(self.seq1)):
            if self.seq1[i] != self.seq2[i]:
                self.distance += 1
        return self.distance
    
    def compute_hamming_similarity(self):
        ''' Function to compute the Hamming similarity between the two sequences. '''
        # We compute the similarity by subtracting twice the distance from the total length
        if self.distance is None:
            self.compute_hamming_distance()
        self.similarity = len(self.seq1) - 2 * self.distance
        return self.similarity
            
            
