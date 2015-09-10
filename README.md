# DNAshapedTFBS

## Python module for TFFM/PSSM + DNA shape classifiers

This module allows for:

1. Training TFFM/PSSM + DNA shape classifiers on ChIP-seq data
2. Applying TFFM/PSSM + DNA shape classifiers on ChIP-seq data

Note that only the best hit per ChIP-sequence is considered in the current
version of the module.

## Dependencies

The module requires:

* python2.7 (and does not work with python3 in its current
version).
* the BioPython module www.biopython.org.
* the [TFFM package](http://cisreg.cmmt.ubc.ca/TFFM/doc/index.html) accessed from your
PYTHONPATH environment variable.
* the [scikit-learn module](http://scikit-learn.org/stable).
* access to bigWig files providing the values of the DNA shape features HelT,
MGW, ProT, and Roll from your genome interest along with the second order
computation these features. Please visit the
[GBshape website](rohsdb.cmb.usc.edu/GBshape).

## Tutorial

You can find some examples of how to run the DNAshapedTFBS.py tool in the script
test.sh provided in the test/ repository of this package.

## Project home page

For information on the source tree, examples, issues, and pull requests, see

    http://github.com/amathelier/DNAshapedTFBS

## Cite

If you use this module, please cite

* A. Mathelier, L. Yang, T.-P. Chiu, R.R. Rohs, and W.W. Wasserman (2015)
Transcription factor binding site prediction in vivo using DNA sequence and
shape features. Submitted.
* A. Mathelier and W.W. Wasserman (2013) The Next Generation of Transcription
Factor Binding Site Predictions. PLoS Comput Biol 9(9):e1003214.
* T.-P. Chiu, L. Yang, T. Zhou, B.J. Main, S.C. Parker, S.V. Nuzhdin, T.D.
Tullius, and R.R. Rohs (2015) GBshape: a genome browser database for DNA shape
annotations. Nucleic Acids Res Jan,43(Database issue):D103-9.
* Pedregosa et al (2011) Scikit-learn: Machine Learning in Python. JMLR 12, 
pp. 2825-2830.
