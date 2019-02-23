import sys
from sequences import *
from sw import smithWaterman
import random
import numpy as np

def getAlignments(scoringMatrix, gapStart, gapExtend):
	tpAlignments = []
	for a,b in posPairs:
		score, matchedA, matchedB = smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2:]
		tpAlignments.append((matchedA, matchedB))
	tnAlignments = []
	for a,b in negPairs:
		score, matchedA, matchedB = smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2:]
		tnAlignments.append((matchedA, matchedB))
	return tpAlignments, tnAlignments

def scoreAlignment(a, b, scoringMatrix, gapStart, gapExtend):
	'''
	given an alignment, just score it.
	'''
	openGapA = False
	openGapB = False
	score = 0
	assert len(a) == len(b)
	for A, B in zip(a, b):
		if A == "-" and openGapA == False:
			score += gapStart
			openGapA = True
		elif A == "-" and openGapA == True:
			score += gapExtend
		elif B == "-" and openGapB == False:
			score += gapStart
			openGapB = True
		elif B == "-" and openGapB == True:
			score += gapExtend
		else:
			score += scoringMatrix[frozenset((A, B))]
			openGapA = False
			openGapB = False
	return score

scoreMatrixD = {}
def scoreMatrix(tpAlignments, tnAlignments, scoringMatrix, gapStart, gapExtend):
	'''
	score our alignments using scoringMatrix.
	our objective function is the sum of TP rates at FP rates of 0.0,
	0.1, 0.2, and 0.3.
	'''
	# First, find the cutoff for each FP rate.
	k = frozenset(scoringMatrix.items())
	if k in scoreMatrixD.keys():
		return scoreMatrixD[k]
	negScores = [scoreAlignment(a,b, scoringMatrix, gapStart, gapExtend) for a,b in tnAlignments]
	cutoffs = np.percentile(negScores, [100, 90, 80, 70])
	# then, find the TP rates for each cutoff.
	posScores = [scoreAlignment(a,b, scoringMatrix, gapStart, gapExtend) for a,b in tpAlignments]
	score = sum([sum(map(lambda x: x > cutoff, posScores))/len(posScores) for cutoff in cutoffs])
	scoreMatrixD[k] = score
	return score

def mutateMatrix(scoringMatrix, mutationChance, mutationAmount):
	'''
	mutate some of the values of the matrix.
	'''
	for key, value in scoringMatrix.items():
		if random.random() < mutationChance:
			scoringMatrix[key] = value + random.gauss(0, mutationAmount)
	return scoringMatrix

def selection(pop, weights):
	'''
	creates a new generation of matrices by sampling with replacement with weights
	'''
	return list(np.random.choice(pop, size=len(pop), p=weights))

def scaleScores(scores, selectivePressure):
	'''
	scales scores so the largest/smallest is 10^selectivePressure and the sum is 1.
	'''
	# Is there really not a more elegant way to do this
	oldMin = min(scores)
	oldMax = max(scores)
	newMin = 1
	newMax = 10**selectivePressure
	newScores =  [(newMax - newMin)*(x - oldMin)/(oldMax - oldMin) + newMin for x in scores]
	return [s/sum(newScores) for s in newScores]

def optimizeMatrix_geneticAlg(scoringMatrix, mutationChance, mutationAmount,
	selectivePressure, N, totalStepsToStop, stepsWithNoImprovement,
	librarySize, gapStart, gapExtend, tpAlignments, tnAlignments):
	'''
	optimize an alignment score matrix using a genetic algorithm.

	stops if we hit totalItersToStop iterations, or if we don't see
	a new objective function value in the top 10 in
	stepsWithNoImprovement steps.
	'''
	# initialize population
	pop = [scoringMatrix.copy() for i in range(N)]
	library = {}
	libraryMin = float("-inf")
	# loop, keeping track of change in objective function
	objectiveMeans = []
	noImprovementCounter = 0
	while True:
		noImprovementCounter += 1
		# mutate
		pop = [mutateMatrix(m, mutationChance, mutationAmount) for m in pop]
		# score
		scores = [scoreMatrix(tpAlignments, tnAlignments, m, gapStart, gapExtend) for m in pop]
		# bank good ones
		for score, mat in zip(scores, pop):
			if score > libraryMin:
				if score not in library.keys():
					noImprovementCounter = 0
					library[score] = mat.copy()
					if len(library) > librarySize:
						libraryMin = min(library.keys())
						del library[libraryMin]
						libraryMin = min(library.keys())
		# assess.
		objectiveMean = sum(scores)/N
		# print (len(objectiveMeans), objectiveMean, max(scores), max(library.keys()))
		objectiveMeans.append(objectiveMean)
		# Should we stop?
		if len(objectiveMeans) > totalStepsToStop or noImprovementCounter > stepsWithNoImprovement:
			break
		# we're continuing. Let's make a new generation.
		# Rescale objective to get weights
		weights = scaleScores(scores, selectivePressure)
		# select a new generation
		pop = selection(pop, weights)[:-len(library)] + [x.copy() for x in library.values()]
		# seed with good ones
	return pop, scores, library, objectiveMeans
