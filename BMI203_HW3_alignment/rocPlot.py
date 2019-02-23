import sys
from .sequences import *
from .sw import smithWaterman
import matplotlib.pyplot as plt
import sklearn.metrics as skm

def makeRocPlot(scoringMatrix, gapStart, gapExtend, name):
	'''
	given a scoring matrix and gap/extend parameters, finds the ROC using our
	known TP and TN data, and saves a PDF plot. 
	'''
	# Score TPs and TNs
	posScores = [smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2] for a,b in posPairs]
	negScores = [smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2] for a,b in negPairs]
	# build vectors that sklearn's ROC stuff expects: a single vector of TP/TN identities,
	# and a single vector of scores
	allScores = posScores + negScores
	allClasses = [1]*len(posScores) + [0]*len(negScores)
	fpr, tpr, threshold = skm.roc_curve(allClasses, allScores)
	roc_auc = skm.auc(fpr, tpr)
	#plot. Adapted from https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
	fig, axes = plt.figure(), plt.axes()
	plt.title('Receiver Operating Characteristic for %s\nwith gap start %s and extend %s' % \
	 (name, gapStart, gapExtend))
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	axes.set_aspect('equal', 'box')
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	#plt.show()
	fig.savefig("ROC_%s.pdf" % name)

def makeRocPlot_normScores(scoringMatrix, gapStart, gapExtend, name):
	'''
	given a scoring matrix and gap/extend parameters, finds the ROC using our
	known TP and TN data, when scaling by smaller string length, and saves a PDF plot. 
	'''
	# Score TPs and TNs
	posScores = [smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2]/min(len(a), len(b)) for a,b in posPairs]
	negScores = [smithWaterman(a,b, scoringMatrix, gapStart, gapExtend)[2]/min(len(a), len(b)) for a,b in negPairs]
	# build vectors that sklearn's ROC stuff expects: a single vector of TP/TN identities,
	# and a single vector of scores
	allScores = posScores + negScores
	allClasses = [1]*len(posScores) + [0]*len(negScores)
	fpr, tpr, threshold = skm.roc_curve(allClasses, allScores)
	roc_auc = skm.auc(fpr, tpr)
	#plot. Adapted from https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
	fig, axes = plt.figure(), plt.axes()
	plt.title('Receiver Operating Characteristic for %s\nwith gap start %s and extend %s\nwith normalized scores' % \
	 (name, gapStart, gapExtend))
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	axes.set_aspect('equal', 'box')
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	#plt.show()
	fig.savefig("ROC_normScores_%s.pdf" % name)

def makeRocPlot_givenScores(posScores, negScores, gapStart, gapExtend, name):
	'''
	given a list of TP and TP scores builds a ROC and PDF plots.
	'''
	# build vectors that sklearn's ROC stuff expects: a single vector of TP/TN identities,
	# and a single vector of scores
	allScores = posScores + negScores
	allClasses = [1]*len(posScores) + [0]*len(negScores)
	fpr, tpr, threshold = skm.roc_curve(allClasses, allScores)
	roc_auc = skm.auc(fpr, tpr)
	#plot. Adapted from https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
	fig, axes = plt.figure(), plt.axes()
	plt.title('Receiver Operating Characteristic for %s\nwith gap start %s and extend %s' % \
	 (name, gapStart, gapExtend))
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	axes.set_aspect('equal', 'box')
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	#plt.show()
	fig.savefig("ROC_normScores_%s.pdf" % name)