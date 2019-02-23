# BMI203 HW3: alignment

[![Build
Status](https://travis-ci.com/garrett-wong/BMI203_HW3.svg?branch=master)](https://travis-ci.com/garrett-wong/BMI203_HW3)

Alignment.

## structure

__main__.py calls the rest of the code to work through each of the tasks for the homework:

  * sequences.py reads in the pairs of true positive and false positive sequences;
  * matrices.py reads in the substitution matrices;
  * sw.py contains my implementation of Smith-Waterman local alignment;
  * optimizeGaps.py contains code to find the false positive rate for a range of gap opening and extension costs;
  * rocPlot.py contains code to create ROC plots;
  * optimizeMatrix.py contains code to run a genetic algorithm to improve a substitution matrix's performance

allSequences.txt is a list of the relative filepaths for the sequences. 
The rest of the repository (scoring matrices, Negpairs.txt and Pospairs.txt, sequences/, etc. is unchanged.)

## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `BMI203_HW3/__main__.py`) can be run to perform
subtasks:

```
python -m BMI203_HW3 -A   # find optimal gap penalties
python -m BMI203_HW3 -B   # compare scoring matrices with ROC curves
python -m BMI203_HW3 -C   # try length-normalized scoring
python -m BMI203_HW3 -D   # optimize blosum50
python -m BMI203_HW3 -E   # optimize matio
```

## testing

Run

```
python -m pytest
```

from the root directory of this project.

