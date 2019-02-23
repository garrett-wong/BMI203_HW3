import sys

def readFasta(filename):
	'''
	reads in a protein sequence from a fasta file
	'''
	seq = ""
	for line in open(filename, "r"):
		if line[0] == ">":
			continue
		seq = seq + line.rstrip().upper()
	return seq

allSeqFilenames = []
for line in open("allSequences.txt", "r"):
	allSeqFilenames.append(line.rstrip())

allSeqs = {}
for filename in allSeqFilenames:
	allSeqs[filename] = readFasta(filename)

posPairFs = [line.rstrip().split() for line in open("Pospairs.txt")]
posPairs = [(allSeqs[fileA], allSeqs[fileB]) for fileA, fileB in posPairFs]

negPairFs = [line.rstrip().split() for line in open("Negpairs.txt")]
negPairs = [(allSeqs[fileA], allSeqs[fileB]) for fileA, fileB in negPairFs]