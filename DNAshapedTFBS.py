#!/usr/bin/python
#*-* coding: utf-8 *-*

""" Train and apply TFFM/PSSM + DNAshape classifiers. """

import os
PATH = os.path.dirname(os.path.realpath(__file__))
import sys
# Local environment
sys.path.append('{0}/../TFFM/'.format(PATH))
from hit_module import HIT
from sklearn.externals import joblib
from argparsing import *
from constants import BWTOOL, DNASHAPEINTER
from shapes import *
from utils import *


def find_pssm_hits(pssm, seq_file):
    """ Predict hits in sequences using a PSSM. """
    from operator import itemgetter
    import math
    import Bio.SeqIO
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA as unambiguousDNA
    hits = []
    for record in Bio.SeqIO.parse(seq_file, "fasta", generic_dna):
        record.seq.alphabet = unambiguousDNA()
        scores = [(pos, ((score - pssm.min) / (pssm.max - pssm.min)))
                  for pos, score in pssm.search(record.seq, pssm.min) if not
                  math.isnan(score)]
        pos_maxi, maxi = max(scores, key=itemgetter(1))
        strand = "+"
        if pos_maxi < 0:
            strand = "-"
            pos_maxi = pos_maxi + len(record.seq)
        hits.append(HIT(record, pos_maxi + 1, pos_maxi + pssm.length, strand,
                        maxi))
    return hits


def find_tffm_hits(xml, seq_file):
    """ Predict hits in sequences using a TFFM. """
    import sys
    sys.path.append("/raid6/amathelier/TFFM+DNAshape/bin/TFFM/")
    import tffm_module
    from constants import TFFM_KIND  # TFFM-framework
    from hit_module import HIT
    tffm = tffm_module.tffm_from_xml(xml, TFFM_KIND.FIRST_ORDER)
    return tffm.scan_sequences(seq_file, only_best=True)


def fit_classifier(fg_train_hits, fg_train_shapes, bg_train_hits,
                   bg_train_shapes, extension=0):
    """ Fit the classifier to the training data. """
    from sklearn.ensemble import GradientBoostingClassifier
    fg_train = combine_hits_shapes(fg_train_hits, fg_train_shapes, extension)
    bg_train = combine_hits_shapes(bg_train_hits, bg_train_shapes, extension)
    data, classification = construct_classifier_input(fg_train, bg_train)
    classifier = GradientBoostingClassifier()
    classifier.fit(data, classification)
    return classifier


def construct_classifier_input(foreground, background):
    """ Make list of classes for foreground and background. """
    classes = [1.0] * len(foreground) + [0.0] * len(background)
    return foreground + background, classes


def make_predictions(clf, tests, hits, thr):
    """ Predict hits from the tests using the classifier. """
    predictions = {'peak_id': [], 'start': [], 'end': [], 'strand': [],
                   'sequence': [], 'proba': []}
    for indx, proba in enumerate(clf.predict_proba(tests)):
        if proba[1] >= thr:
            hit = hits[indx]
            predictions['peak_id'].append(hit.seq_record.name)
            predictions['start'].append(hit.start)
            predictions['end'].append(hit.end)
            predictions['strand'].append(hit.strand)
            if hit.strand == '-':
                sequence = ''.join(
                    hit.seq_record.seq[
                        hit.start - 1:hit.end].reverse_complement())
            else:
                sequence = ''.join(hit.seq_record[hit.start - 1:hit.end])
            predictions['sequence'].append(sequence)
            predictions['proba'].append(proba[1])
    return predictions


def print_predictions(predictions, output):
    """ Print the predictions in the output file. """
    import pandas as pd
    pd_predictions = pd.DataFrame(predictions)
    pd.set_option('display.max_rows', len(pd_predictions))
    with open(output, 'w') as stream:
        stream.write('{0}\n'.format(pd_predictions.to_string(
            index=False, columns=['peak_id', 'start', 'end', 'strand',
                                  'sequence', 'proba'])))


def apply_classifier(hits, argu):
    """ Apply the DNAshape-based classifier. """
    # Two options, 1) doing sequence by sequence but it means doing a lot of
    # bwtool calls and I/O, 2) doing all sequences as one batch but means that
    # we need to associate back probas to hits. I choose 2) to reduce I/O.

    hits_shapes = get_shapes(pssm_hits, argu.in_bed, argu.helt, argu.mgw,
            argu.prot, argu.roll, argu.helt2, argu.mgw2, argu.prot2, argu.roll2,
            argu.extension, argu.scaled)
    assert len(hits) == len(hits_shapes[0])
    classifier = joblib.load(argu.classifier)
    tests = combine_hits_shapes(hits, hits_shapes, argu.extension)
    # Need to print the results by associating the probas to the hits
    predictions = make_predictions(classifier, tests, hits, argu.threshold)
    print_predictions(predictions, argu.output)


def tffm_apply_classifier(argu):
    """ Apply the TFFM + DNA shape classifier. """
    hits = find_tffm_hits(argu.tffm_file, argu.in_fasta)
    apply_classifier(hits, argu)


def pssm_apply_classifier(argu):
    """ Apply the TFFM + DNA shape classifier. """
    if argu.jasparid:
        pssm = get_jaspar_pssm(argu.jasparid)
    else:
        pssm = get_jaspar_pssm(argu.jasparfile, False)
    hits = find_pssm_hits(pssm, argu.in_fasta)
    apply_classifier(hits, argu)


def train_classifier(fg_hits, bg_hits, argu):
    """ Train the DNAshape-based classifier. """
    fg_shapes = get_shapes(fg_hits, argu.fg_bed, argu.helt, argu.mgw, argu.prot,
            argu.roll, argu.helt2, argu.mgw2, argu.prot2, argu.roll2,
            argu.argu.extension, argu.scaled)
    bg_shapes = get_shapes(bg_hits, argu.bg_bed, argu.helt, argu.mgw, argu.prot,
            argu.roll, argu.helt2, argu.mgw2, argu.prot2, argu.roll2,
            argu.extension, argu.scaled)
    classifier = fit_classifier(fg_hits, fg_shapes, bg_hits, bg_shapes,
                                argu.extension)
    joblib.dump(classifier, '{0}.pkl'.format(argu.output))


def tffm_train_classifier(argu):
    """ Train a TFFM + DNA shape classifier. """
    fg_hits = find_tffm_hits(argu.tffm_file, argu.fg_fasta)
    bg_hits = find_tffm_hits(argu.tffm_file, argu.bg_fasta)
    train_classifier(fg_hits, bg_hits, argu)


def pssm_train_classifier(argu):
    """ Train a PSSM + DNA shape classifier. """
    if argu.jasparid:
        pssm = get_jaspar_pssm(argu.jasparid)
    else:
        pssm = get_jaspar_pssm(argu.jasparfile, False)
    fg_hits = find_pssm_hits(pssm, argu.fg_fasta)
    bg_hits = find_pssm_hits(pssm, argu.bg_fasta)
    train_classifier(fg_hits, bg_hits, argu)


##############################################################################
#                               MAIN
##############################################################################
if __name__ == "__main__":
    arguments = arg_parsing()
    arguments.func(arguments)
