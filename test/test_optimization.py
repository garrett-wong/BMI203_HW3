import sys
from sequences import *
from optimizeMatrix import *
from matrices import blosum50

tpAlignments, tnAlignments = getAlignments(scoringMatrix=blosum50,
	gapStart=-10, gapExtend=-1)

def test_get_alignments():
	# I am gonna track that example from before. I align them all
	# already; I make sure it still matches.
	assert tpAlignments[5] == ("RAECIQR-GVSPSQAQGLGSNLVTE", "RKRKIDRDAVLNMWQQGLGASHISK")

def test_scoring():
	# I redo the scoring from SW
	testA, testB = tpAlignments[5]
	assert scoreAlignment(testA, testB, blosum50, -10, -1) == 35

# I skipped testing my objective function code: it combs through a TON of data
# that I can't make toy examples of easily.

def test_mutate_matrix():
	# I want to make sure that I don't make more copies than necessary:
	M = blosum50.copy()
	M2 = mutateMatrix(M, 1, 1)
	assert M2 == M
	# But that I do change the matrix:
	assert M2 != blosum50

def test_selection():
	# Very random and hard to test, but I want to make sure that my weights
	# are working and that I'm sampling with replacement: 
	L = ["a", "b", "c"]
	w = [0, 0, 1]
	assert selection(L, w) == ["c", "c", "c"]

def test_scale_scores():
	# just make sure I'm doing what I said: largest/smallest is 10^selectivePressure
	# and the sum is 1.
	a,b = scaleScores([1, 2], 1)
	assert a+b == 1
	assert b/a == 10

# For the actual meat of the function, I mostly can only use the running output to
# verify that it's running correctly. Here, I make sure that all our data stays
# the right length:
def test_geneticAlg_skeleton():
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		blosum50, 1, 1, 1, 10, 5, 5, 3, -10, -1, tpAlignments, tnAlignments)
	assert len(pop) == 10
	assert len(scores) == 10
	assert len(library) == 3
	assert len(objectiveMeans) == 5 + 1 # there's an extra -inf at the beginning so that it
										# starts off increasing

