def scale01(values, mini=None, maxi=None, tol=1e-6):
    """ Scale the values in [0, 1]. """
    from numpy import amax, amin
    if not mini:
        mini = amin(values)
    if not maxi:
        maxi = amax(values)
    scaled_values = [(val - mini) / (maxi - mini + tol) for val in values]
    return scaled_values, mini, maxi


def not_na(item):
    """ Remove NAs and empty values. """
    return not(item == "NA" or item == "")


def extended_hit_pos(hit, peak_chrom, peak_start, extension=0):
    """ Extend the hit by 'extension' nt to compute DNAshape features. """
    start = peak_start + hit.start - extension - 2  # BED
    end = peak_start + hit.end + extension - 1  # BED
    return peak_chrom, start, end


def print_extended_hits(hits, positions, extension=0):
    """
    Print the extended hits to a temporary bed file.

    :returns: The name of the temporary file.
    :rtype: str

    """

    import tempfile
    import os
    fdescr, tmp_file = tempfile.mkstemp()
    os.close(fdescr)
    with open(tmp_file, 'w') as stream:
        for hit in hits:
            if hit:
                identifier = hit.seq_record.id
                peak_chrom, peak_start, _ = positions[identifier]
                chrom, start, end = extended_hit_pos(hit, peak_chrom,
                                                     peak_start, extension)
                if not chrom.startswith("chr"):
                    chrom = "chr{0}".format(chrom)
                if hit.score >= 0.0 and hit.score <= 1.0:
                    stream.write("{0}\t{1:d}\t{2:d}\t{3}\t{4:d}\t{5}\n".format(
                        chrom, start, end, identifier, int(hit.score * 100),
                        hit.strand))
                else:
                    stream.write("{0}\t{1:d}\t{2:d}\t{3}\t{4:d}\t{5}\n".format(
                        chrom, start, end, identifier, 0,
                        hit.strand))
    return tmp_file


def get_positions_from_bed(bed_file):
    """ Get the positions of the sequences described in the bed file. """
    with open(bed_file) as stream:
        positions = {}
        for line in stream:
            spl = line.split()
            positions[spl[3]] = (spl[0], eval(spl[1]) + 1, eval(spl[2]))
    return positions


def get_jaspar_pssm(jaspar, bool_id=True):
    """ 
    
    Construct the PSSM from the JASPAR ID or JASPAR formatted file.

    We assume that we are using profiles from the CORE JASPAR db when providing
    a JASPAR ID. Hence the JASPAR ID should starts with 'MA'.
    If a filename is provided, we assume that the TF binding profile is using
    the JASPAR format as documented in the Bio.motifs.jaspar BioPython module.
    
    """
    import Bio.motifs
    if bool_id:
        from Bio.motifs.jaspar.db import JASPAR5
        jaspar_db_host = "vm5.cmmt.ubc.ca"
        jaspar_db_name = "JASPAR_2014"
        jaspar_db_user = "jaspar_r"
        jaspar_db_pass = ""
        jdb = JASPAR5(host=jaspar_db_host, name=jaspar_db_name,
                      user=jaspar_db_user, password=jaspar_db_pass)
        motif = jdb.fetch_motif_by_id(jaspar)
        motif.pseudocounts = Bio.motifs.jaspar.calculate_pseudocounts(motif)
    else:
        with open(jaspar) as handle:
            motif = Bio.motifs.read(handle, 'jaspar')
    return motif.pssm
