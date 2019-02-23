import sys

def readSubstitutionMatrix(filename):
	'''
	reads in a substitution matrix
	'''
	matrix = {}
	readHeader = False
	header = []
	rowi = 0
	for line in open("BMI203_HW3_alignment/" + filename, "r"):
		# skip comments
		if line[0] == "#":
			continue
		# read in line and separate on whitespace
		fields = line.rstrip().split()
		# read in header
		if readHeader == False:
			header = fields
			colLabels = fields
			readHeader = True
			continue
		# make sure there's the right number of fields
		assert len(fields) == len(header)
		row = header[rowi]
		# step through the fields and add matrix entries
		for col, score in zip(header, fields):
			key = frozenset((row, col))
			if key not in matrix:
				matrix[key] = int(score)
			else:
				assert matrix[key] == int(score)
		# we're reading the next row now...
		rowi += 1
	return matrix

blosum50 = readSubstitutionMatrix("BLOSUM50")
matio = readSubstitutionMatrix("MATIO")

matrixNames = ["BLOSUM50", "BLOSUM62", "MATIO", "PAM100", "PAM250"]
matrices = map(readSubstitutionMatrix, matrixNames)