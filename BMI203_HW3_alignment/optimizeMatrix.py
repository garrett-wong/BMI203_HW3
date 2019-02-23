import sys
from .sequences import *
from .sw import smithWaterman
import random
import numpy as np

def getAlignments(scoringMatrix, gapStart, gapExtend):
	'''
	fetches true positive and true negative alignments using the scoring matrix, gapStart,
	and gapExtend costs.

	returns a list of true positive alignments (tuples of strings) and the same for true negatives.
	'''
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

	takes as input a and b; aligned strings,
	as well as a scoring matrix and gap/extend parameters.
	returns a score for this alignment.
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

	takes as input true positive and true negative tuples of aligned strings, as
	well as the scoring matrix and gap/extend costs for alignment.

	returns the objective function.
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
	each value has mutationChance probability [0, 1] of mutating,
	and will change by adding gaussian noise with mean 0 and standard deviation mutationAmount.

	returns the SAME matrix, modified.
	'''
	for key, value in scoringMatrix.items():
		if random.random() < mutationChance:
			scoringMatrix[key] = value + random.gauss(0, mutationAmount)
	return scoringMatrix

def selection(pop, weights):
	'''
	creates a new generation of matrices by sampling with replacement with weights.

	takes as input a list of objects and a list of selection probabilities, and
	returns a listof the same size.
	'''
	return list(np.random.choice(pop, size=len(pop), p=weights))

def scaleScores(scores, selectivePressure):
	'''
	scales scores so the largest/smallest is 10^selectivePressure and the sum is 1.

	takes as input a list of scores and returns a new list of scaled scores.
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
	a new objective function value in the top (librarySize) in
	stepsWithNoImprovement steps.

	takes as input:
	scoringMatrix
	mutationChance		probability each entry in a scoring matrix will mutate
	mutationAmount		stdev of mutation gaussian
	selectivePressure	most fit individual is scaled to 10^selectivePressure * least fit
	n                   pop size
	totalStepsToStop	runtime step cutoff
	stepsWithNoImprovement
						step cutoff for seeing no new entries in our library of best matrices
	librarySize			number of best matrices to remember and re-seed population with
	gapStart			alignment params
	gapExtend
	tpAlignments		true positives
	tnAlignments		true negatives

	returns the final population of matrices, the scores for those matrices at the last step,
	the library of best matrices, and a list of the mean objective function value at each step.
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
