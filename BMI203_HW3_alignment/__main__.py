import sys
import optimizeGaps
import rocPlot
from matrices import matrices, matrixNames, blosum50, matio
from optimizeMatrix import getAlignments, optimizeMatrix_geneticAlg

# 1.1: find optimal gap penalties to minimize FP when TP is set to 0.7.
# record, graph.
if sys.argv[1][0:2] == '-A':
	optimizeGaps.optimizeGaps()

# 1.2: comparing scoring matrices with ROC curves
if sys.argv[1][0:2] == '-B':
	gapStart = -7
	gapExtend = -4
	for name, matrix in zip(matrixNames, matrices):
		rocPlot.makeRocPlot(matrix, gapStart, gapExtend, name)

# 1.3: trying length-normalized scores
if sys.argv[1][0:2] == '-C':
	rocPlot.makeRocPlot_normScores(blosum50, gapStart, gapExtend, "blosum50")

# 2.2: optimizing blosum50
# get alignments
if sys.argv[1][0:2] == '-D':
	gapStart = -7
	gapExtend = -4
	scoringMatrix = blosum50
	tpAlignments, tnAlignments = getAlignments(scoringMatrix=blosum50,
		gapStart=-7, gapExtend=-4)
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		scoringMatrix=blosum50, mutationChance=0.5, mutationAmount=0.1,
		selectivePressure=1, N=100, totalStepsToStop=1000,
		stepsWithNoImprovement=100, librarySize=10, gapStart=-7,
		gapExtend=-4, tpAlignments, tnAlignments)
	bestMat_blosum50 = library[max(library.keys())]
	# make ROC curve with given alignments
	negScores = [scoreAlignment(a,b, bestMat_blosum50, gapStart, gapExtend) for a,b in tnAlignments]
	posScores = [scoreAlignment(a,b, bestMat_blosum50, gapStart, gapExtend) for a,b in tpAlignments]
	rocPlot.makeRocPlot_givenScores(posScores, negScores, gapStart, gapExtend, "geneticOptimized_blosum50")
	# make ROC curve by rescoring.
	rocPlot.makeRocPlot(bestMat_blosum50, gapStart, gapExtend, "geneticOptimized_blosum50_rescored")

# 2.3: optimizing matio
	# get alignments
if sys.argv[1][0:2] == '-E':
	gapStart = -7
	gapExtend = -4
	scoringMatrix = blosum50
	tpAlignments, tnAlignments = getAlignments(scoringMatrix=blosum50,
		gapStart=-7, gapExtend=-4)
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		scoringMatrix=matio, mutationChance=0.5, mutationAmount=0.1,
		selectivePressure=1, N=100, totalStepsToStop=1000,
		stepsWithNoImprovement=100, librarySize=10, gapStart=-7,
		gapExtend=-4, tpAlignments, tnAlignments)
	bestMat_matio = library[max(library.keys())]
	# make ROC curve with given alignments
	negScores = [scoreAlignment(a,b, bestMat_matio, gapStart, gapExtend) for a,b in tnAlignments]
	posScores = [scoreAlignment(a,b, bestMat_matio, gapStart, gapExtend) for a,b in tpAlignments]
	rocPlot.makeRocPlot_givenScores(posScores, negScores, gapStart, gapExtend, "geneticOptimized_matio")
	# make ROC curve by rescoring.
	rocPlot.makeRocPlot(bestMat_matio, gapStart, gapExtend, "geneticOptimized_matio_rescored")

