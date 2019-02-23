import sys
import numpy as np

def readSubstitutionMatrix(filename):
	'''
	reads in a substitution matrix
	'''
	matrix = {}
	readHeader = False
	header = []
	rowi = 0
	for line in open(filename, "r"):
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

def smithWaterman(a, b, subMat, gapStartCost, gapExtendCost):
	'''
	uses smith-waterman to find a maximal local alignment.
	'''
	# initialize our matrices and trace matrices
	match = np.zeros(shape=(len(a) + 1, len(b) + 1))
	traceMatch = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)]
	Y = np.zeros(shape=(len(a) + 1, len(b) + 1))
	traceY = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)]
	X = np.zeros(shape=(len(a) + 1, len(b) + 1))
	traceX = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)]
	# fill first rows/columns of matrices to avoid aligns that run into weird edges
	for i in range(len(a) + 1):
		Y[i, 0] = float("-inf")
		if i>0: 
			X[i, 0] = gapStartCost + gapExtendCost*(i-1)
	for j in range(len(b) + 1):
		X[0, j] = float("-inf")
		if j>0:
			Y[0, j] = gapStartCost + gapExtendCost*(j-1)
	# fill in matrices
	maxSeen = (float("-inf"), (-1, -1)) # (score, (i,j))
	for i in range(1, len(a) + 1):
		for j in range(1, len(b) + 1): 
			extendX = X[i-1, j] + gapExtendCost  # extend gap in b
			openX = match[i-1, j] + gapStartCost # open gap in b
			if openX > extendX:
				X[i,j] = openX
				traceX[i][j] = "openX"
			else:
				X[i,j] = extendX
				traceX[i][j] = "extendX"
			extendY = Y[i, j-1] + gapExtendCost  # extend gap in a
			openY = match[i, j-1] + gapStartCost # open gap in a
			if openY > extendY: 
				Y[i,j] = openY
				traceY[i][j] = "openY"
			else:
				Y[i,j] = extendY
				traceY[i][j] = "extendY"
			match[i,j] = match[i-1, j-1] + subMat[frozenset((a[i-1], b[j-1]))]
			traceMatch[i][j] = "match"
			if Y[i,j] > match[i,j]:
				match[i,j] = Y[i,j]
				traceMatch[i][j] = "closeY"
			if X[i,j] > match[i,j]:
				match[i,j] = X[i,j]
				traceMatch[i][j] = "closeX"
			if match[i,j] > maxSeen[0]:
				maxSeen = (match[i,j], (i, j))
	# traceback
	i,j = maxSeen[1]
	end = (i,j)
	matrix, traceMatrix = match, traceMatch
	score = matrix[i,j]
	bestScore = score
	matchedA = ""
	matchedB = ""
	while score > 0:
		if traceMatrix[i][j] == "match":
			matchedA = a[i-1] + matchedA
			matchedB = b[j-1] + matchedB
			i,j = i-1, j-1
		elif traceMatrix[i][j] == "closeX":
			matrix, traceMatrix = X, traceX
		elif traceMatrix[i][j] == "closeY":
			matrix, traceMatrix = Y, traceY
		elif traceMatrix[i][j] == "openY":
			matchedA = "-" + matchedA
			matchedB = b[j-1] + matchedB
			i,j = i, j-1
			matrix, traceMatrix = match, traceMatch
		elif traceMatrix[i][j] == "extendY":
			matchedA = "-" + matchedA
			matchedB = b[j-1] + matchedB
			i,j = i, j-1
		elif traceMatrix[i][j] == "openX":
			matchedA = a[i-1] + matchedA
			matchedB = "-" + matchedB
			i,j = i-1, j
			matrix, traceMatrix = match, traceMatch
		elif traceMatrix[i][j] == "extendX":
			matchedA = a[i-1] + matchedA
			matchedB = "-" + matchedB
			i,j = i-1, j
		else: print("we should never get here.")
		score = matrix[i,j]
	start = (i, j)
	#print(start, end, bestScore)
	#print(matchedA)
	#print(matchedB)
	return start, end, bestScore, matchedA, matchedB



