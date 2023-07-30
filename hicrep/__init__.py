import os
import numpy as np
import math
import sys
import warnings
import hictkpy
from hicrep.utils import (
    getSubCoo,
    trimDiags, meanFilterSparse, varVstran,
    resample
    )
from hicrep.hicrep import (
    sccOfDiag, hicrepSCC
    )

def main(*args):
    import argparse
    import subprocess
    import re

    np.random.seed(10)

    parser = argparse.ArgumentParser()
    parser.add_argument("file1", type=str,
                        help="First .hic or .mcool contact file")
    parser.add_argument("file2", type=str,
                        help="Second .hic or .mcool contact file")
    parser.add_argument("fout", type=str,
                        help="Output results to this file. Output format would be\
                        one column of scc scores for each chromosome")
    parser.add_argument("--binSize", type=int, default=0,
                        help="Use this to select the bin size from the input mcool\
                        file. Default to 0, meaning that the inputs are treated as\
                        single-binsize .cool files")
    parser.add_argument("--h", type=int, required=True,
                        help="Smooth the input contact matrices using a 2d mean\
                        filter with window size of 1 + 2 * value. This should\
                        be set according to the bin size. For example, you can try the\
                        following settings: --binSize=10000 --h=20,\
                        --binSize=25000 --h=10, --binSize=40000 --h5. Beware that\
                        these examples might not work in all cases and the user\
                        should adjust them according to the specific application")
    parser.add_argument("--dBPMax", type=int, required=True,
                        help="Only consider contacts at most this number of bp away\
                        from the diagonal. For human genome, the value of\
                        5000000 was used in the original HiCRep paper.")
    parser.add_argument("--bDownSample", action='store_true', default=False,
                        help="Down sample the input with more contact counts to\
                        the the same number of counts as the other input with less\
                        contact counts. If turned off, the input matrices will be\
                        normalized by dividing the counts by their respective total\
                        number of contacts.")
    parser.add_argument("--chrNames", type=str, nargs='*', default=[],
                        help="Only compute the SCC scores on this subset of\
                        chromosomes whose names are provided. The output SCC\
                        scores will be ordered as the input chromosome names\
                        here")
    parser.add_argument("--excludeChr", type=str, nargs='*', default=['chrM', 'M'],
                        help="Exclude chromosomes from the SCC score calculations.\
                        Mitochondrial chromosomes named \"M\" are excluded by\
                        default. The output SCC scores will be ordered as the\
                        chromosomes in the input Cooler files by removing those\
                        chromosomes provided here")

    args = parser.parse_args()

    assert not (args.excludeChr != ['chrM', 'M'] and len(args.chrNames) > 0), f"""
        Please use --chrNames OR --excludeChr arguments but not both.
        """

    header = "#"+" ".join(sys.argv)+"\n"

    # Check if current script is under revision control
    gitls = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__)) +
                             ' && git ls-files --error-unmatch ' + __file__,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, encoding='utf-8')
    if not gitls.stderr.read():
        gitrev = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__))
                                  + ' && git rev-parse HEAD --abbrev-ref HEAD',
                                  shell=True, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, encoding='utf-8')
        if not gitrev.stderr.read():
            p = re.compile("\n")
            header += "# @rev " + p.sub("\n# @branch ",
                                        gitrev.stdout.read().strip()) + "\n"

    fout = args.fout
    binSize = args.binSize
    h = args.h
    dBPMax = args.dBPMax
    bDownSample = args.bDownSample
    chrNames = args.chrNames
    excludeChr = set(args.excludeChr)

    if len(excludeChr) != len(args.excludeChr):
        warnings.warn(f"""
            Duplicate excludeChr found. Please remove them in --excludeChr""")

    file1 = hictkpy.File(args.file1, binSize)
    file2 = hictkpy.File(args.file2, binSize)

    scc = hicrepSCC(file1, file2, h, dBPMax, bDownSample,
                    chrNames if len(chrNames) > 0 else None,
                    excludeChr if len(excludeChr) > 0 else None)

    np.savetxt(fout, scc, "%30.15e", header=header)
