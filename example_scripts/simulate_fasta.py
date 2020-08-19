#!/usr/bin/env python
import os
import sys
import time

import warnings
import deeplc
import pandas as pd
import ursgal
from loguru import logger
from numpy.random import choice, normal, random, seed, uniform
from pyqms.chemical_composition import ChemicalComposition
from tqdm import tqdm
import tensorflow as tf
import smiter
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('darkgrid')

# seed(42)
warnings.simplefilter("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.get_logger().setLevel('ERROR')


def main(fasta_file, evidences=None):
    peak_props = {}
    cc = ChemicalComposition()

    if evidences is not None:
        ev = pd.read_csv(evidences)
        train_df = pd.DataFrame()
        train_df['seq'] = ev['Sequence']
        train_df['modifications'] = ev['Modifications'].fillna('')
        train_df['tr'] = ev['Retention Time (s)'] / 60
        dlc = deeplc.DeepLC()
        logger.info(f'Train DeepLC on {evidences}')
        dlc.calibrate_preds(seq_df=train_df)

    start = time.time()
    logger.info("Parsing fasta and calculating peak Properties")
    no_peptides = 0
    with open(fasta_file) as fin:
        for _id, seq in tqdm(
            ursgal.ucore.parse_fasta(fin),
            total=6691,
            # total=75004,
            desc="Calculating peak properties",
        ):
            peptides = ursgal.ucore.digest(seq, ("KR", "C"), no_missed_cleavages=True)

            r = uniform(1e4, 1e6)
            for p in peptides:

                if len(p) > 30 or len(p) < 10:
                    continue
                r_number = random()
                if r_number < 0.6:
                    continue
                no_peptides += 1
                peak_function = choice(["gauss", "gauss_tail"], p=[0, 1.0])
                if peak_function == "gauss":
                    peak_params = {"sigma": 2}
                elif peak_function == "gauss_tail":
                    peak_params = {"sigma": 2}
                try:
                    cc.use(p)
                except:
                    continue
                # logger.info(f'Add {p}')
                # r = random()
                scale_factor = uniform(r*1, r*2)
                charge_state = choice([2, 3, 4], 1, p=[0.6, 0.3, 0.1])[0]
                peak_width = uniform(0.3, 0.5)
                formula = cc.hill_notation_unimod()
                predicted_start = uniform(0, 120 * 60)
                # TODO run RT prediction before

                peak_props[p] = {
                    "chemical_formula": f"+{formula}",
                    "trivial_name": p,
                    "charge": charge_state,
                    "scan_start_time": predicted_start,
                    "peak_width": peak_width,
                    "peak_function": peak_function,
                    "peak_params": peak_params,
                    "peak_scaling_factor": normal(loc=scale_factor, scale=1),
                }

    logger.info(f"Predict RT for {no_peptides} simulated molecules")
    pred_df = pd.DataFrame()
    pred_df['seq'] = list(peak_props.keys())
    pred_df['modifications'] = ''
    rts = dlc.make_preds(seq_df=pred_df)
    pred_df['tr'] = rts

    # logger.info("Updating predicted RTs in peak properties")
    for i, line in pred_df.iterrows():
        peak_props[line['seq']]['scan_start_time'] = line['tr']
    duration = time.time() - start


    # d = {}
    # from pprint import pprint
    # for p in peak_props.keys():
    #     if peak_props[p]['peak_function'] == 'gauss_tail':
    #         pprint(peak_props[p])
    #         _t = round(peak_props[p]['scan_start_time'] + peak_props[p]['peak_width'])
    #         if _t not in d:
    #             d[_t] = 1
    #         else:
    #             d[_t] += 1
    # t = []
    # count = []
    # for rt in sorted(d.keys()):
    #     t.append(rt)
    #     count.append(d[rt])
    # plt.plot(t, count); plt.show()
    mzml_params = {
        "gradient_length": 120,  # in minutes
        "ms_rt_diff": 0.001,  # in minutes
        "min_intensity": 0
    }

    # _t = []
    # t = []
    # i = []
    # tree = smiter.synthetic_mzml.genereate_interval_tree(peak_props)
    # for x in range(0, mzml_params['gradient_length']):
    #     # x /= 2
    #     l = len(tree.at(x))
    #     t.append(x)
    #     i.append(l)
    #     # _t.append(x)
    # plt.plot(t, i); plt.show()
    # sns.distplot(i, bins=100);plt.show()
    # breakpoint()

    logger.info(f"Simulated {len(peak_props)} peptides in {duration} seconds")

    smiter.lib.peak_properties_to_csv(
        peak_props, "/mnt/Data/Proteomics/testing/test.csv"
    )

    fragmentor = smiter.fragmentation_functions.PeptideFragmentor()
    fragmentor = smiter.fragmentation_functions.PeptideFragmentorPyteomics()
    noise_injector = smiter.noise_functions.JamssNoiseInjector()
    file = "/mnt/Data/Proteomics/testing/test.mzML"

    start = time.time()
    mzml_path = smiter.synthetic_mzml.write_mzml(
        file, peak_props, fragmentor, noise_injector, mzml_params
    )
    logger.info(f'Wrote mzML  to {file}')
    # logger.info(f"Wrote mzML in {time.time() - start} seconds")


if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main(sys.argv[1], sys.argv[2])
