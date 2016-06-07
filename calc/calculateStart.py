class CalculateStart:
    #calculates the start of the codons based on the start of the ribo seq read.

    #the start calculated from the 5 prime of the sequence

    def calc5Prime(self, start):
        return (int(start)+1) % 3

    #the start calculated from the 3 prime of the sequence
    def calc3Prime(self, start, sequence):
        return ((int(start)+1) + len(sequence)-1) % 3

    #the start calculated from the middel of the sequence
    def calcMid(self, start, sequence):
        return ((int(start)+1) + (len(sequence)/2)) % 3