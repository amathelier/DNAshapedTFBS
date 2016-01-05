from DNAshapedTFBS import tffm_train_classifier
from DNAshapedTFBS import pssm_train_classifier
from DNAshapedTFBS import tffm_apply_classifier
from DNAshapedTFBS import pssm_apply_classifier


def tffm_train_arg_parsing(subparsers):
    """ Train the TFFM + DNA shape classifier. """
    help_str = "Train the TFFM + DNA shape classifier."
    parser_t = subparsers.add_parser('trainTFFM', help=help_str)
    parser_t.add_argument('-T', '--tffmfile', required=True, dest='tffm_file',
            action='store', type=str, help='TFFM XML file.')
    help_str = 'Input fasta file containing the foreground sequences.'
    parser_t.add_argument('-i', '--fg_fasta', required=True, type=str,
                          dest='fg_fasta', action='store', help=help_str)
    help_str = 'Input bed file w/ positions of foreground sequences.'
    parser_t.add_argument('-I', '--fg_bed', required=True, type=str,
                          dest='fg_bed', action='store', help=help_str)
    help_str = 'Input fasta containing the background sequences.'
    parser_t.add_argument('-b', '--bg_fasta', required=True, type=str,
                          dest='bg_fasta', action='store', help=help_str)
    help_str = 'Input bed file w/ positions of background sequences.'
    parser_t.add_argument('-B', '--bg_bed', required=True, type=str,
                          dest='bg_bed', action='store', help=help_str)
    parser_t.add_argument('-o', '--outfile', required=True, type=str,
                          dest='output', action='store',
                          help='Output base name (.pkl will be added).')
    help_str = 'List of 1st-order bigWig files as\n'
    help_str += '-1 HelT ProT MGW Roll'
    help_str += ' (note that we assume such an order)'
    parser_t.add_argument('-1', '--1storder', required=True, type=str,
                          dest='first_shape', action='store', nargs=4,
                          help=help_str)
    help_str = 'List of 2nd-order bigWig files as\n'
    help_str += '-2 HelT2 ProT2 MGW2 Roll2'
    help_str += ' (note that we assume such an order)'
    parser_t.add_argument('-2', '--2ndorder', required=False, type=str,
                          dest='second_shape', action='store', nargs=4,
                          default=[], help=help_str)
    parser_t.add_argument('-n', '--scaled', required=False, dest='scaled',
                          action='store_true', default=False,
                          help='Scale DNAshape values in [0, 1]')
    help_str = 'Extension to be considered around TFBSs with DNAshapes'
    help_str += ' (default:0).'
    parser_t.add_argument('-e', '--extension', required=False, type=int,
                          dest='extension', action='store', default=0,
                          help=help_str)
    parser_t.set_defaults(func=tffm_train_classifier)


def tffm_apply_arg_parsing(subparsers):
    """ Apply the TFFM + DNA shape classifier. """
    help_str = 'Apply the TFFM + DNA shape classifier.'
    parser_a = subparsers.add_parser('applyTFFM', help=help_str)
    parser_a.add_argument('-T', '--tffmfile', required=True, dest='tffm_file',
            action='store', type=str, help='TFFM XML file.')
    help_str = 'Input fasta file containing the sequences.'
    parser_a.add_argument('-i', '--input_fasta', required=True, type=str,
                          dest='in_fasta', action='store', help=help_str)
    help_str = 'Input bed file w/ positions of the sequences.'
    parser_a.add_argument('-I', '--input_bed', required=True, type=str,
                          dest='in_bed', action='store', help=help_str)
    parser_a.add_argument('-c', '--classifier', required=True, type=str,
                          dest='classifier', action='store',
                          help='Classifier (.pkl file).')
    parser_a.add_argument('-o', '--outfile', required=True, type=str,
                          dest='output', action='store',
                          help='Output result file')
    help_str = 'List of 1st-order bigWig files as\n'
    help_str += '-1 HelT ProT MGW Roll'
    help_str += ' (note that we assume such an order)'
    parser_a.add_argument('-1', '--1storder', required=True, type=str,
                          dest='first_shape', action='store', nargs=4,
                          help=help_str)
    help_str = 'List of 2nd-order bigWig files as: '
    help_str += '-2 HelT2 ProT2 MGW2 Roll2'
    help_str += ' (note that we assume such an order)'
    parser_a.add_argument('-2', '--2ndorder', required=False, type=str,
                          dest='second_shape', action='store', nargs=4,
                          default=[], help=help_str)
    parser_a.add_argument('-n', '--scaled', required=False, dest='scaled',
                          action='store_true', default=False,
                          help='Scale DNAshape values in [0, 1]')
    help_str = 'Extension to be considered around TFBSs with DNAshapes'
    help_str += ' (default:0).'
    parser_a.add_argument('-e', '--extension', required=False, type=int,
                          dest='extension', action='store', default=0,
                          help=help_str)
    help_str = 'Probability threshold to predict a hit (default=0.5).'
    parser_a.add_argument('-v', '--threshold', required=False,
                          dest='threshold', action='store', default=0.5,
                          type=float, help=help_str)
    parser_a.set_defaults(func=tffm_apply_classifier)


def pssm_train_arg_parsing(subparsers):
    """ Train the PSSM + DNA shape classifier. """
    help_str = "Train the PSSM + DNA shape classifier."
    parser_t = subparsers.add_parser('trainPSSM', help=help_str)
    jaspar_grp = parser_t.add_mutually_exclusive_group(required=True)
    help_str = 'JASPAR ID corresponding to the TF '
    help_str += 'binding profile to be used.'
    jaspar_grp.add_argument('-j', '--jaspar', type=str, dest='jasparid',
            action='store', help=help_str)
    help_str = 'JASPAR file containing the TF binding profile in the '
    help_str += 'JASPAR format.'
    jaspar_grp.add_argument('-f', '--jasparfile', type=str, dest='jasparfile',
            action='store', help=help_str)
    help_str = 'Input fasta file containing the foreground sequences.'
    parser_t.add_argument('-i', '--fg_fasta', required=True, type=str,
                          dest='fg_fasta', action='store', help=help_str)
    help_str = 'Input bed file w/ positions of foreground sequences.'
    parser_t.add_argument('-I', '--fg_bed', required=True, type=str,
                          dest='fg_bed', action='store', help=help_str)
    help_str = 'Input fasta containing the background sequences.'
    parser_t.add_argument('-b', '--bg_fasta', required=True, type=str,
                          dest='bg_fasta', action='store', help=help_str)
    help_str = 'Input bed file w/ positions of background sequences.'
    parser_t.add_argument('-B', '--bg_bed', required=True, type=str,
                          dest='bg_bed', action='store', help=help_str)
    parser_t.add_argument('-o', '--outfile', required=True, type=str,
                          dest='output', action='store',
                          help='Output base name (.pkl will be added).')
    help_str = 'List of 1st-order bigWig files as\n'
    help_str += '-1 HelT ProT MGW Roll'
    help_str += ' (note that we assume such an order)'
    parser_t.add_argument('-1', '--1storder', required=True, type=str,
                          dest='first_shape', action='store', nargs=4,
                          help=help_str)
    help_str = 'List of 2nd-order bigWig files as\n'
    help_str += '-2 HelT2 ProT2 MGW2 Roll2'
    help_str += ' (note that we assume such an order)'
    parser_t.add_argument('-2', '--2ndorder', required=False, type=str,
                          dest='second_shape', action='store', nargs=4,
                          default=[], help=help_str)
    parser_t.add_argument('-n', '--scaled', required=False, dest='scaled',
                          action='store_true', default=False,
                          help='Scale DNAshape values in [0, 1]')
    help_str = 'Extension to be considered around TFBSs with DNAshapes'
    help_str += ' (default:0).'
    parser_t.add_argument('-e', '--extension', required=False, type=int,
                          dest='extension', action='store', default=0,
                          help=help_str)
    parser_t.set_defaults(func=pssm_train_classifier)


def pssm_apply_arg_parsing(subparsers):
    """ Apply the PSSM + DNA shape classifier. """
    help_str = 'Apply the PSSM + DNA shape classifier.'
    parser_a = subparsers.add_parser('applyPSSM', help=help_str)
    jaspar_grp = parser_a.add_mutually_exclusive_group(required=True)
    help_str = 'JASPAR ID corresponding to the TF '
    help_str += 'binding profile to be used.'
    jaspar_grp.add_argument('-j', '--jaspar', type=str, dest='jasparid',
            action='store', help=help_str)
    help_str = 'JASPAR file containing the TF binding profile in the '
    help_str += 'JASPAR format.'
    jaspar_grp.add_argument('-f', '--jasparfile', type=str, dest='jasparfile',
            action='store', help=help_str)
    help_str = 'Input fasta file containing the sequences.'
    parser_a.add_argument('-i', '--input_fasta', required=True, type=str,
                          dest='in_fasta', action='store', help=help_str)
    help_str = 'Input bed file w/ positions of the sequences.'
    parser_a.add_argument('-I', '--input_bed', required=True, type=str,
                          dest='in_bed', action='store', help=help_str)
    parser_a.add_argument('-c', '--classifier', required=True, type=str,
                          dest='classifier', action='store',
                          help='Classifier (.pkl file).')
    parser_a.add_argument('-o', '--outfile', required=True, type=str,
                          dest='output', action='store',
                          help='Output result file')
    help_str = 'List of 1st-order bigWig files as\n'
    help_str += '-1 HelT ProT MGW Roll'
    help_str += ' (note that we assume such an order)'
    parser_a.add_argument('-1', '--1storder', required=True, type=str,
                          dest='first_shape', action='store', nargs=4,
                          help=help_str)
    help_str = 'List of 2nd-order bigWig files as\n'
    help_str += '-2 HelT2 ProT2 MGW2 Roll2'
    help_str += ' (note that we assume such an order)'
    parser_a.add_argument('-2', '--2ndorder', required=False, type=str,
                          dest='second_shape', action='store', nargs=4,
                          default=[], help=help_str)
    parser_a.add_argument('-n', '--scaled', required=False, dest='scaled',
                          action='store_true', default=False,
                          help='Scale DNAshape values in [0, 1]')
    help_str = 'Extension to be considered around TFBSs with DNAshapes'
    help_str += ' (default:0).'
    parser_a.add_argument('-e', '--extension', required=False, type=int,
                          dest='extension', action='store', default=0,
                          help=help_str)
    help_str = 'Probability threshold to predict a hit (default=0.5).'
    parser_a.add_argument('-v', '--threshold', required=False,
                          dest='threshold', action='store', default=0.5,
                          type=float, help=help_str)
    parser_a.set_defaults(func=pssm_apply_classifier)


def arg_parsing():
    """ Parse the subcommand along with its arguments. """

    descr = '''
    Train or apply the DNAshape-based classifiers to a set of fasta sequences.
    '''
    import argparse
    parser = argparse.ArgumentParser(
        description=descr,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(
        help='Train or apply a TFFM/PSSM + DNA shape classifier',
        title='Subcommands', description='Valid subcommands')
    tffm_train_arg_parsing(subparsers)
    tffm_apply_arg_parsing(subparsers)
    pssm_train_arg_parsing(subparsers)
    pssm_apply_arg_parsing(subparsers)
    argu = parser.parse_args()
    return argu
