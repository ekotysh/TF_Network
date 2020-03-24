#!/usr/bin/env python3
"""
__author__  = "Eduard Kotysh"
__email__   = "kotysh@wustl.edu"

This script takes a set of elite enhancers from GeneHancer and finds TF
binding sites that are present in each enhancer. Then visualizes a network
of TF -> Enh Region (BS1, BS2) -> Target Gene. See README.md for more info.

Usage: python3 main.py
"""

import concurrent
import os
import statistics
import sys
import re
import logging.config
import configparser
from concurrent.futures import ALL_COMPLETED, as_completed
from concurrent.futures.thread import ThreadPoolExecutor
from pathlib import Path
from tqdm import tqdm
from enhancers import extract_enhancers, populate_enhancer_sequences
from motifs import get_motif, find_bs_motifs

# LOGGER
logging.config.fileConfig('logging.conf')
logger = logging.getLogger(__name__)

# LOAD CONFIG
config = configparser.ConfigParser()
config.read("config_dev.ini")
file_conf = config['FILES']


###############################################################
# Classes
###############################################################

class Gene:
    ensembl_gene_id = ""
    external_gene_name = ""
    gene_biotype = ""
    description = ""
    sequence = ""

    def __init__(self, column_list):
        self.ensembl_gene_id = column_list[0]
        self.external_gene_name = column_list[1]

        if (len(column_list) > 2):
            self.gene_biotype = column_list[2]

        if (len(column_list) > 3):
            self.description = column_list[3]


###############################################################
# Functions
###############################################################

# genes.ENSG.tbl is in this format:
# ensembl_gene_id	external_gene_name	gene_biotype	description
def preload_gene_names(filepath):
    logger.debug("Preloading gene names")
    gene_names_data = open(filepath, "r")

    # read and ignore the header line
    gene_names_data.readline()

    gene_dict = {}
    for line in gene_names_data:
        line_list = re.split(r'\t+', line.strip())  # split line on tabs into a list
        ensemble_name = line_list[0]  # use ensemble id as key

        # create a gene object out of values in this line
        gene_obj = Gene(line_list)

        # add the gene object as value in the dict, key'd by ensemble name.
        gene_dict[ensemble_name] = gene_obj

    gene_names_data.close()
    return gene_dict


###############################################################
# MAIN
###############################################################

# check that the correct number of arguments was given (1)
if (len(sys.argv) != 1):
    sys.exit(__doc__)

gene_names_dict = preload_gene_names(file_conf['gene_names_file'])
enhancers_dict = extract_enhancers(file_conf['gene_hancer_file'], gene_names_dict)
populate_enhancer_sequences(enhancers_dict, file_conf['gene_hancer_seq_file'])

# Read the PFMs one file at a time, and generate motifs
binding_motifs = []
binding_length_list = []
logger.debug("Generating motifs")

for filename in os.listdir(file_conf['hocomoco_pcm_mono_dir']):
    if filename.endswith(".pcm"):
        file_path = Path(file_conf['hocomoco_pcm_mono_dir']) / filename
        motif = get_motif(file_path)
        if (motif is not None):
            binding_motifs.append(motif)
            binding_length_list.append(len(str(motif.consensus)))

# print how many enhancers and binding motifs we loaded
logger.info("enhancers: " + str(len(enhancers_dict.keys())))
logger.info("binding_motifs: " + str(len(binding_motifs)))

# now go through each enhancer in a greedy sliding window fashion
# and generate binding scores
enhancer_bs_matches = {}

count = 0
num_matches_list = []
enhancer_length_list = []

# create a thread pool
executor = ThreadPoolExecutor(max_workers=(int)(config['DEFAULT']['threads']))
futures = []

logger.debug("Searching for motifs...")
for enh_key in enhancers_dict.keys():
    # get the enhancer
    enhancer = enhancers_dict.get(enh_key)
    enhancer_length_list.append(len(enhancer.sequence))

    # spin a thread for each enhancer region (using the pool of X thread workers)
    # each thread will find binding motifs within one enhancer region
    f = executor.submit(find_bs_motifs, enhancer, binding_motifs)
    futures.append(f)

# log progress as tasks complete
for f in tqdm(as_completed(futures), total=len(futures)):
    pass

# wait on the threads to finish
concurrent.futures.wait(futures, timeout=None, return_when=ALL_COMPLETED)

# collect the results
logger.debug("Collecting results")
for future in futures:
    hits = future.result()
    if (len(hits) > 0):
        enhancer_bs_matches[enh_key] = hits
        num_matches_list.append(len(hits))

avg_enh_length = statistics.mean(enhancer_length_list)
avg_matches_per_enh = statistics.mean(num_matches_list)
avg_bs_length = statistics.mean(binding_length_list)

logger.info("avg enhancer length: " + str(avg_enh_length) + " bp")
logger.info("avg bs length: " + str(avg_bs_length) + " bp")
logger.info("avg matches per enhancer: " + str(avg_matches_per_enh) + " binding sites")
# TODO: how many enhancers did each TF bind to?

