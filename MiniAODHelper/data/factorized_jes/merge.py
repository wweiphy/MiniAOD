#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to merge CSV SF histograms into a single file.
Type "./merge.py --help" to see all arguments.

Example:
> ls
csv_rwt_fit_hf_v2_final_2017_6_7_JESAbsoluteMPFBias.root
csv_rwt_fit_hf_v2_final_2017_6_7_JESAbsoluteScale.root
csv_rwt_fit_hf_v2_final_2017_6_7_JESFlavorQCD.root
csv_rwt_fit_hf_v2_final_2017_6_7_Nominal.root

> ./merge.py -file-prefix csv_rwt_fit_hf_v2_final_2017_6_7_
output file   : csv_rwt_fit_hf_v2_final_2017_6_7_all.root
nominal file  : csv_rwt_fit_hf_v2_final_2017_6_7_Nominal.root
jes variations: JESAbsoluteMPFBias,JESAbsoluteScale,JESFlavorQCD
start merging ...
merge variation JESAbsoluteMPFBias
merge variation JESAbsoluteScale
merge variation JESFlavorQCD
done

> ls
csv_rwt_fit_hf_v2_final_2017_6_7_JESAbsoluteMPFBias.root
csv_rwt_fit_hf_v2_final_2017_6_7_JESAbsoluteScale.root
csv_rwt_fit_hf_v2_final_2017_6_7_JESFlavorQCD.root
csv_rwt_fit_hf_v2_final_2017_6_7_Nominal.root
csv_rwt_fit_hf_v2_final_2017_6_7_all.root # this is the merged one
"""


import os
import re
import sys
import glob
import shutil
import contextlib


def main(directory, output_file, file_prefix, nominal_name, jes_name):
    # check and update some args
    if not output_file:
        output_file = file_prefix + "all.root"
    output_file = os.path.normpath(output_file)
    directory = os.path.normpath(directory)

    # determine the nominal file and perform some checks
    nominal_file = os.path.join(directory, file_prefix + nominal_name + ".root")
    if not os.path.exists(nominal_file):
        abort("nominal file '%s' not existing" % nominal_file)
    elif output_file == nominal_file:
        abort("output file must not be the nominal file")
    elif os.path.exists(output_file):
        abort("output file '%s' already exists" % output_file)
    print("output file   : " + output_file)
    print("nominal file  : " + nominal_file)

    # determine jes variations
    jes_variations = []
    for file_path in glob.glob(os.path.join(directory, file_prefix + "*.root")):
        if file_path == nominal_file:
            continue
        variation = os.path.basename(file_path)[len(file_prefix):-5]
        jes_variations.append((variation, file_path))
    if not jes_variations:
        abort("no jes variations found")
    else:
        print("jes variations: " + ",".join(v for v, p in jes_variations))

    # copy the nominal file to the output file
    shutil.copy(nominal_file, output_file)

    # start merging
    print("start merging ...")
    n = 0
    with tfile_guard(output_file, mode="UPDATE") as f_out:
        # loop through all variations and collect histograms
        for variation, file_path in jes_variations:
            print("merge variation " + variation)
            with tfile_guard(file_path) as f_var:
                # loop trough all histograms and select the ones that end with <jes_name>(Up|Down)
                for tkey in f_var.GetListOfKeys():
                    m = re.match("^(.+)%s(Up|Down)$" % jes_name, tkey.GetName())
                    if not m:
                        continue

                    prefix = m.group(1)
                    direction = m.group(2)

                    f_out.cd()
                    clone = f_var.Get(tkey.GetName()).Clone(prefix + variation + direction)
                    clone.Write()
                    n += 1

                f_out.Write()

    print("done, copied %i histograms" % n)


def abort(msg=None, code=1):
    if msg:
        sys.stderr.write(msg + "\n")
        sys.stderr.flush()
    sys.exit(code)


@contextlib.contextmanager
def tfile_guard(*tfiles, **kwargs):
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch()

    mode = kwargs.get("mode", "READ")

    tfiles = list(tfiles)
    for i, f in enumerate(tfiles):
        if isinstance(f, basestring):
            tfiles[i] = ROOT.TFile.Open(f, mode)
    try:
        yield tfiles[0] if len(tfiles) == 1 else tfiles
    finally:
        for tfile in tfiles:
            if tfile and tfile.IsOpen():
                tfile.Close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    # setup argument parsing
    parser = ArgumentParser(description="CSV SF merger")
    parser.add_argument("--directory", "-d", default=os.getcwd(), help="location of the files to "
            "merge, default: pwd")
    parser.add_argument("--output-file", "-o", default=None, help="output file, default: "
            "<file-prefix>all.root")
    parser.add_argument("--file-prefix", "-p", default="csv_rwt_fit_hf_v2_final_2017_6_7_",
            help="input file prefix, default: csv_rwt_fit_hf_v2_final_2017_6_7_")
    parser.add_argument("--nominal-name", "-n", default="Nominal", help="name of the nominal file, "
            "default: Nominal")
    parser.add_argument("--jes-name", "-j", default="JES", help="name of the JES variation in "
            "histogram names, default: JES")
    args = parser.parse_args()

    # merge
    main(args.directory, args.output_file, args.file_prefix, args.nominal_name, args.jes_name)
