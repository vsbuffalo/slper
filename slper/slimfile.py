## slim.py -- parse slim results in a tabular format
import numpy as np
import glob
import os
import re
import sys
import random
from collections import namedtuple, OrderedDict, defaultdict
from itertools import product
import pandas as pd
from scipy.sparse import coo_matrix

SlimFreqs = namedtuple('SlimFreqs', ('params', 'positions', 'samples',
                                     'freqs', 'times'))
SlimStats = namedtuple('SlimStats', ('params', 'stats'))

def split_keyval(keyval):
    key, val = keyval.split('=')
    if key.isalpha():
        return key, val
    return key, float(val)

def parse_params(param_str):
    "Parse the parameter string."
    return dict([split_keyval(keyval) for keyval in param_str.split(';')])


def parse_slim_stats(filename, delim='\t', verbose=False):
    """
    Parse a tabular representation of population statistics.

    The return value is a (params, stats) tuple.
    """
    with open(filename) as fp:
        line = next(fp)  # graph first line
        assert(line.startswith('#'))
        params = parse_params(line[1:])
        # reset the read position
        fp.seek(0)
        statsdf = pd.read_csv(fp, comment='#', header=0, delimiter=delim)
    return SlimStats(params, statsdf)


def parse_slim_ragged_freqs(filename, delim='\t'):
    """
    Parse a ragged array of
      'gen', 'id;pos;freq' x number of polymorphic mutations

    Assumes generations are in order. Uses sparse matrices to efficiently
    load these data into a matrix.
    """
    muts = dict()
    row, col, data = [], [], []
    first_gen = None
    gens = []
    # mock genomic intervals
    positions = []
    # we keep our own internal IDs here
    with open(filename) as fp:
        line = next(fp)
        params = parse_params(line[1:])
        for line in fp:
            fields = line.strip().split(delim)
            gen = int(fields[0])
            gens.append(gen)
            if first_gen is None:
                first_gen = gen
            for mut in fields[1:]:
                mfs = mut.split(';')
                mid, pos, freq = int(mfs[0]), int(mfs[1]), float(mfs[2])
                muts[mid] = pos
                positions.append(pos)
                row.append(gen-first_gen)
                col.append(mid)
                data.append(freq)
    # remap the timepoints
    timepoints = {t: i for i, t in enumerate(set(row))}
    row = [timepoints[t] for t in row]
    # now map their IDs to our IDs
    key_map = dict((k, i) for i, k in enumerate(muts.keys()))
    new_col = [key_map[mid] for mid in col]
    assert(len(new_col) == len(row) == len(data))
    loci = OrderedDict((muts[mid], i) for i, mid in enumerate(col))
    # Even though SLIM doesn't really have concenpt of multiple chromosomes,
    # we use the defaultdict(OrderedDict) approach as in handling SyncFiles.
    # The chromosome name is None
    loci_dict = defaultdict(OrderedDict)
    loci_dict[None] = loci
    return SlimFreqs(params, positions, gens,
                     coo_matrix((data, (row, new_col)),
                                shape=(len(gens), len(new_col))).toarray(),
                     timepoints)


def parse_slim_freqs(filename, delim='\t', min_prop_samples=0,
                     missing_val=np.nan, verbose=False):
    """
    Parse a tabular representation of output frequencies.

    In my SLiM results, I preallocate a ngens x nloci matrix, and fill it with
    results as generations complete. Empty values are assigned a value of -1,
    which is replaced with `missing_val` (default: np.nan).

    The first row is a parameter string beginning with #, with 'key=val'
    parameter values.

    The return value is a (params, generations, frequency matrix) tuple.
    """
    with open(filename) as fp:
        next(fp)  # skip parameters
        header = next(fp).strip().split(delim)
        loci = np.array(header[1:], dtype='u4')  # grab loci
        # get the number of columns
        width = len(fp.readline().strip().split(delim))
        # reset the read position
        fp.seek(0)
        # parse the params
        line = next(fp).strip()
        if not line.startswith('#'):
            msg = ("error: SLiM results file does not begin with #-prefixed "
                   "string containing parameters.")
            raise ValueError(msg)
        params = parse_params(line[1:])
        fp.seek(0)
        next(fp); next(fp) # skip parameters, skip header
        # now, read the entire matrix
        #dtypes = ', '.join(['i4'] + ['f4'] * (width - 1))
        dtypes = 'float64'
        mat = np.loadtxt(fp, delimiter=delim, dtype=dtypes)
    # extract out the generation (samples here) and convert types
    samples = mat[:,0].astype('i4')
    # drop generation column, and convert -1 to nan
    mat = np.delete(mat, 0, 1)
    # convert -1s for non-poymorphic site to 0s
    mat[mat < 0] = missing_val
    # remove loci that are never polymorphic
    min_nsamples = int(min_prop_samples*mat.shape[0])
    keep_cols = np.logical_not(np.isnan(mat)).sum(0) > min_nsamples
    if verbose:
        msg = (f"pruning matrix from {mat.shape[1]} to "
               f"{keep_cols.sum()} columns (threshold: >{min_nsamples} samples)")
        print(msg)
    mat = mat[:, keep_cols]
    loci = loci[keep_cols]
    return SlimFreqs(params, loci, samples, mat)


def output_filename(input_filename, suffix, ext='.tsv'):
    _, inext = os.path.splitext(input_filename)
    return input_filename.replace(inext, f'-{suffix}{ext}')

if __name__ == '__main__':
    import sys
    import tempdata

    file = sys.argv[2]
    cmd = sys.argv[1]
    simdata = parse_slim_tsv(file, min_prop_samples=0.)
    d = tempdata.TemporalFreqs(simdata.freqs, simdata.samples, simdata.loci)

    if cmd == 'cov':
        print("calculating covariances...")
        d.calc_covs()
        print("writing covariances...")
        d.write_covs(output_filename(file, 'cov'), long=True)
    elif cmd == 'freq':
        print("writing frequencies...")
        d.write_freqs(output_filename(file, 'freq'))
    else:
        print("usage: [cov,freq] input_file.txt")
        sys.exit(1)
