import sys

# I use this to read in all my sequences.
# posPairs and negPairs are the relevant lists of tuples of sequence strings.

def readFasta(filename):
	'''
	reads in a protein sequence from a fasta file.
	returns a string.
	'''
	seq = ""
	for line in open(filename, "r"):
		if line[0] == ">":
			continue
		seq = seq + line.rstrip().upper()
	return seq

allSeqFilenames = []
for line in open("BMI203_HW3_alignment/allSequences.txt", "r"):
	allSeqFilenames.append(line.rstrip())

allSeqs = {}
for filename in allSeqFilenames:
	allSeqs[filename] = readFasta("BMI203_HW3_alignment/"+filename)
posPairFs = [line.rstrip().split() for line in open("BMI203_HW3_alignment/Pospairs.txt")]
posPairs = [(allSeqs[fileA], allSeqs[fileB]) for fileA, fileB in posPairFs]

negPairFs = [line.rstrip().split() for line in open("BMI203_HW3_alignment/Negpairs.txt")]
negPairs = [(allSeqs[fileA], allSeqs[fileB]) for fileA, fileB in negPairFs]