
class SubstitutionMatrix:
    def __init__(self, filename):
        '''Store substitution scores from a matrix file in a dictionary.'''
        self.matrix = {}
        with open(filename, 'r') as file:
            # The first line of the file is different, we need to read it separately and save the keys it contains
            keys_1 = file.readline().split()
            for line in file:
                # For each other line, we separate the row key in the first position from the values in the other positions
                splitted_line = line.split()
                key_2 = splitted_line[0]
                values = splitted_line[1:]
                
                for key_1, value in zip(keys_1, values):
                    # For each matched pair of column key / corresponding value, build the combined key (col + row) and save the value in the dictionary
                    key = key_1 + key_2
                    self.matrix[key] = int(value)
                    

    def __getitem__(self, key):
        '''This class 'extends' the old one defined hamming.py'''
        a, b = key
        k = a.upper() + b.upper()
        if not k in self.matrix:
            raise KeyError(f"Key {key} not in the Substitution Matrix")
        return self.matrix[k]



