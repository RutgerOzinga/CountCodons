import pysam
from calculateStart import CalculateStart

#This script counts the amount of codon occurances in each codon position. As well
#as looking at the distribution of starting locations of the ribo-seq reads (so if the start is divided by 3 what would be left)
class CountStarts(object):
    dict5Prime = {0: 0, 1: 0, 2: 0}
    codonDict2 = {}
    codonDict3 = {}
    codonDict4 = {}
    codonDict5 = {}
    codonDict6 = {}
    path = ""

    def __init__(self):
        #path to the bam file opened with pysam.
        self.path = pysam.AlignmentFile("/home/rutgero/Desktop/BN_fasta_files/ERR653341.bam", "rb")
        #path the a new file.
        self.newPath = open("/home/rutgero/Desktop/tabsepCodonFile1_version2.txt","w")

    def readFile(self):
        count = 0
        calc = CalculateStart()
        #for each line i the bam file
        for line in self.path:
            count += 1
            #each time the counter can be divided by 100000 without leaving a rest value the count is printed (to see progress).
            if count % 100000 == 0:
                print(count)
            splittedLine = str(line).split("\t")
            #only take the lines of the bam file that do not have the flag 4 which means that the read is not mapped.
            if splittedLine[1] != "4":
                self.count5Starts(calc.calc5Prime(splittedLine[3]))
                #if the start modulo by 3 is 1 and the cigar is 29M (meaning that 29 of the bases of the read matched to 29 of the bases in the reference)
                if calc.calc5Prime(splittedLine[3]) == 1 and splittedLine[5] == "29M":
                    codon2 = splittedLine[9][2:5]
                    codon3 = splittedLine[9][5:8]
                    codon4 = splittedLine[9][8:11]
                    codon5 = splittedLine[9][11:14]
                    codon6 = splittedLine[9][14:17]
                    #adds the second codon to the dict containing the second codons if its not in there and otherwise
                    # adds one to the codon count.
                    if (codon2 in self.codonDict2):
                        self.codonDict2[codon2] += 1
                    else:
                        self.codonDict2.update({codon2:1})

                    # adds the third codon to the dict containing the third codons if its not in there and otherwise
                    # adds one to the codon count.
                    if (codon3 in self.codonDict3):
                        self.codonDict3[codon3] += 1
                    else:
                        self.codonDict3.update({codon3: 1})

                    # adds the fourth codon to the dict containing the fourth codons if its not in there and otherwise
                    # adds one to the codon count.
                    if (codon4 in self.codonDict4):
                        self.codonDict4[codon4] += 1
                    else:
                        self.codonDict4.update({codon4: 1})

                    # adds the fifth codon to the dict containing the fifth codons if its not in there and otherwise
                    # adds one to the codon count.
                    if (codon5 in self.codonDict5):
                        self.codonDict5[codon5] += 1
                    else:
                        self.codonDict5.update({codon5: 1})
                    # adds the sixth codon to the dict containing the sixth codons if its not in there and otherwise
                    # adds one to the codon count.
                    if (codon6 in self.codonDict6):
                        self.codonDict6[codon6] += 1
                    else:
                        self.codonDict6.update({codon6: 1})


            else:
                pass
        #writes a header in the new file.
        self.newPath.write("codon_name\tcodon_2\tcodon_3\tcodon_4\tcodon_5\tcodon_6")
        #for each codon in the codonDict2 keyset.
        for codon in self.codonDict2.keys():
            #if the codon contains a N skip it
            if "N" in str(codon):
                pass
            #else get the values of the codon from each dict and add them to the file along with the codon
            else:
                self.newPath.write(codon + "\t" + str(self.codonDict2[codon])
                               + "\t" + str(self.codonDict3[codon])
                               + "\t" + str(self.codonDict4[codon])
                               + "\t" + str(self.codonDict5[codon])
                               + "\t" + str(self.codonDict6[codon]) + "\n")
        print(self.codonDict2, self.codonDict3, self.codonDict4, self.codonDict5, self.codonDict6 )

    #adds one to the key that corresponds with the start
    def count5Starts(self, start):
        newTotal = self.dict5Prime[start] + 1
        self.dict5Prime[start] = newTotal


if __name__ == "__main__":
    a = CountStarts()
    a.readFile()
