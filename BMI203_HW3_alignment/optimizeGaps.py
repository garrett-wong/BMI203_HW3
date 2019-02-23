import sys
import numpy as np
from .sequences import *
from .matrices import blosum50
from .sw import smithWaterman

TPR = 70
def getFpr(scoringMatrix, gapStart, gapExtend):
	# First, find all the true positive alignment scores.
	# We'll need these to find the cutoff that sets the TP rate
	# at 0.7.
	posScores = [smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2] for a,b in posPairs]
	cutoff = np.percentile(posScores, 100-TPR)
	# Now, find the false positive rate.
	negScores = [smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2] for a,b in negPairs]
	# It's just the count above the cutoff divided by the number of true negatives.
	return sum(map(lambda x: x > cutoff, negScores))/len(negScores)

def optimizeGaps():
	with open("optimizeGapPenaltiesTest.txt", "w") as f:
		for gapStart in range(-1, -21, -1):
			for gapExtend in range(-1, -6, -1):
				print("%s, %s:" % (gapStart, gapExtend))
				fpr = getFpr(blosum50, gapStart, gapExtend)
				print(" ", fpr)
				f.write("\t".join(map(str, [gapStart, gapExtend, fpr])) + "\n")

# Then run the R script to plot.