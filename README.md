# DNAshapedTFBS

## Python module for TFFM/PSSM/4-bits + DNA shape classifiers

This module allows for:

1. Training PSSM/TFFM/4-bits + DNA shape classifiers on ChIP-seq data
2. Applying PSSM/TFFM/4-bits + DNA shape classifiers on ChIP-seq data

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
* the [pandas module](http://pandas.pydata.org).
* access to bigWig files providing the values of the DNA shape features HelT,
MGW, ProT, and Roll from your genome interest along with the second order
computation these features. Please visit the
[GBshape website](rohsdb.cmb.usc.edu/GBshape).

## Tutorial

You can find some examples of how to run the DNAshapedTFBS.py tool in the script
test.sh provided in the test/ repository of this package.

The script feature_importance_heatmap.py plots the heatmap(s) of trained
classifier(s). Note that the current version only works for PSSM/TFFM + DNA
shape classifiers. You can get help on how to use it by typing

python2.7 feature_importance_heatmap.py -h

## Project home page

For information on the source tree, examples, issues, and pull requests, see

    http://github.com/amathelier/DNAshapedTFBS

## Cite

If you use this module, please cite

* A. Mathelier, B. Xin, T.-P. Chiu, L. Yang, R.R. Rohs, and W.W. Wasserman (2016)
DNA shape features improve transcription factor binding site predictions in
vivo. Under revision.
* A. Mathelier and W.W. Wasserman (2013) The Next Generation of Transcription
Factor Binding Site Predictions. PLoS Comput Biol 9(9):e1003214.
* T.-P. Chiu, L. Yang, T. Zhou, B.J. Main, S.C. Parker, S.V. Nuzhdin, T.D.
Tullius, and R.R. Rohs (2015) GBshape: a genome browser database for DNA shape
annotations. Nucleic Acids Res Jan,43(Database issue):D103-9.
* T. Zhou, L. Yang, Y. Lu, I. Dror, AC. Dantas Machado, T. Ghane, R. Di Felice,
and R. Rohs (2013) DNAshape: a method for the high-throughput prediction of DNA
structural features on a genomic scale. Nucleic Acids Research 41(Web Server
issue):W56-62.
* Pedregosa et al. (2011) Scikit-learn: Machine Learning in Python. JMLR 12, 
pp. 2825-2830.
