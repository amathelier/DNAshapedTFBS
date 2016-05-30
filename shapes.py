from utilities import *
from the_constants import BWTOOL, DNASHAPEINTER


def get_scores(in_file, shape=None, scaled=False):
    """ Get DNAshape values on single lines. """
    with open(in_file) as stream:
        scores = []
        for line in stream:
            values = [item for item in line.rstrip().split()[7].split(',')
                      if not_na(item)]
            values = [eval(value) for value in values]
            if scaled:
                mini, maxi = DNASHAPEINTER[shape]
                values, _, _ = scale01(values, mini, maxi)
            scores.append(values)
        return scores


def combine_hits_shapes(hits, shapes, extension=0, binary_encoding=False):
    """ Combine DNA sequence and shape features.
    
    The hit scores (PSSM or TFFM) or 4-bits encoding are combined with DNAshape
    features in vectors for classif.
    
    """
    comb = []
    index = -1
    for hit in hits:
        if hit:
            index += 1
            if shapes:
                hit_shapes = []
                for indx in xrange(len(shapes)):
                    hit_shapes += shapes[indx][index]
                if (not binary_encoding and
                        (len(hit_shapes) ==
                            len(shapes) * (len(hit) + 2 * extension)
                        )
                   ):
                    comb.append([hit.score] + hit_shapes)
                elif (binary_encoding and
                        (len(hit_shapes) == 
                            len(shapes) * (len(hit) / 4 + 2 * extension)
                        )
                     ):
                    comb.append(hit + hit_shapes)
            elif binary_encoding:
                comb.append(hit)
            else:
                comb.append([hit.score])
    return comb


def get_shapes(hits, bed_file, first_shape, second_shape, extension=0,
        scaled=False):
    """ Retrieve DNAshape feature values for the hits. """
    bigwigs = first_shape + second_shape
    import subprocess
    import os
    peaks_pos = get_positions_from_bed(bed_file)
    with open(os.devnull, 'w') as devnull:
        tmp_file = print_extended_hits(hits, peaks_pos, extension)
        shapes = ['HelT', 'ProT', 'MGW', 'Roll', 'HelT2', 'ProT2', 'MGW2',
                'Roll2']
        hits_shapes = []
        for indx, bigwig in enumerate(bigwigs):
            if bigwig:
                out_file = '{0}.{1}'.format(tmp_file, shapes[indx])
                subprocess.call([BWTOOL, 'ex', 'bed', tmp_file, bigwig, out_file],
                                stdout=devnull)
                if indx < 4:
                    hits_shapes.append(get_scores(out_file, shapes[indx], scaled))
                else:
                    hits_shapes.append(get_scores(out_file, shapes[indx]))
        subprocess.call(['rm', '-f', '{0}.HelT'.format(tmp_file),
            '{0}.MGW'.format(tmp_file), '{0}.ProT'.format(tmp_file),
            '{0}.Roll'.format(tmp_file), '{0}.HelT2'.format(tmp_file),
            '{0}.MGW2'.format(tmp_file), '{0}.ProT2'.format(tmp_file),
            '{0}.Roll2'.format(tmp_file),tmp_file])
        return hits_shapes
