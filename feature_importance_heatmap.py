#!/usr/bin/python
#*-* coding: utf-8 *-*

###############################################################################
#
###############################################################################

import argparse
from sklearn.externals import joblib


def load_classifiers(files):
    """ Load the input classifiers. """
    return [joblib.load(infile) for infile in files]


def plot_heatmap(feat_imp, infile, secorder, min_val, max_val):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib
    if secorder:
        len_motif = (len(feat_imp) - 1) / 8
        dico = {'Importance': list([feat_imp[0]]) * len_motif +
                list(feat_imp[1:]), 'Feature': ['hit_score'] * len_motif +
                ['HelT'] * len_motif + ['ProT'] * len_motif + ['MGW'] *
                len_motif + ['Roll'] * len_motif + ['HelT2'] * len_motif +
                ['ProT2'] * len_motif + ['MGW2'] * len_motif + ['Roll2'] *
                len_motif, 'Position': range(len_motif) + range(len_motif) +
                range(len_motif) + range(len_motif) + range(len_motif) +
                range(len_motif) + range(len_motif) + range(len_motif) +
                range(len_motif)}
    else:
        len_motif = (len(feat_imp) - 1) / 4
        dico = {'Importance': list([feat_imp[0]]) * len_motif +
                list(feat_imp[1:]), 'Feature': ['hit_score'] * len_motif +
                ['HelT'] * len_motif + ['ProT'] * len_motif + ['MGW'] *
                len_motif + ['Roll'] * len_motif, 'Position': range(len_motif) +
                range(len_motif) + range(len_motif) + range(len_motif) +
                range(len_motif)}
    map = pd.DataFrame.from_dict(dico)
    map = map.pivot(index='Feature', columns='Position', values='Importance')
    sns.heatmap(map, linewidth=.5, vmin=min_val, vmax=max_val)
    plt.savefig('{0}.svg'.format(infile))
    plt.clf()


def plot_average_heatmap(classifiers, output, secorder, min_val, max_val):
    import pandas as pd
    feat_imp = {}
    for indx, clf in enumerate(classifiers):
        feat_imp[indx] = pd.Series(clf.feature_importances_)
    df = pd.DataFrame(feat_imp)
    plot_heatmap(list(df.mean(1)), output, secorder, min_val, max_val)


def create_heatmap(argu):
    import matplotlib
    matplotlib.use('svg')
    infiles = argu.classif_files
    classifiers = load_classifiers(infiles)
    for indx, clf in enumerate(classifiers):
        plot_heatmap(clf.feature_importances_, infiles[indx], argu.secorder,
                argu.min_val, argu.max_val)
    if argu.output and len(classifiers) > 1:
        plot_average_heatmap(classifiers, argu.output, argu.secorder,
                argu.min_val, argu.max_val)


def arg_parsing():
    """ Parse the arguments. """

    descr = '''
    Plot the heatmap corresponding to the feature importance associated to the
    classifier(s) provided.
    '''
    parser = argparse.ArgumentParser(description=descr, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--classif', required=True, nargs='+',
            dest='classif_files', action='store',
            help='Classifier(s) to be used (.pkl file)')
    help_str='Basename of the output averaged heatmap over multiple '
    help_str += 'classifiers (.svg will be added)'
    parser.add_argument('-a', '--average', required=False, type=str,
            dest='output', action='store', default=None, help=help_str)
    parser.add_argument('-2', '--second', required=False, dest='secorder',
            action='store_true', default=False,
            help='Plot classifier using 2nd order DNA shape features.')
    parser.add_argument('-m', '--min', required=False, dest='min_val',
            action='store', type=float, default=None,
            help='Minimal value for the heat map range')
    parser.add_argument('-M', '--max', required=False, dest='max_val',
            action='store', type=float, default=None,
            help='Maximal value for the heat map range')
    parser.set_defaults(func=create_heatmap)
    argu = parser.parse_args()
    return argu


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    arguments = arg_parsing()
    arguments.func(arguments)
