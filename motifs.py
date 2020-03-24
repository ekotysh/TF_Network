import logging.config

from Bio import motifs

from error_msg import ERROR_MSG_IO

# LOGGER
logging.config.fileConfig('logging.conf')
logger = logging.getLogger(__name__)

# CONSTS
HOCOMOCO_PARSE_FORMAT = "pfm-four-columns"


# Parses the pcm file into a Biopython motif
def get_motif(pcm_file):
    try:
        with open(pcm_file) as handle:
            m = motifs.parse(handle, HOCOMOCO_PARSE_FORMAT)
            if (m is not None and len(m) > 0):
                return m[0]
            else:
                return None
    except IOError as e:
        logger.error(ERROR_MSG_IO.format(pcm_file, e.strerror))


# Find exact matches of the consensus sequences of motifs in the enhancer
def find_bs_motifs(enhancer, bs_motifs):
    if enhancer is None or bs_motifs is None or len(bs_motifs) == 0:
        return []

    hits = []
    for m in bs_motifs:
        # search for exact matches of the consensus sequence of the motif
        hits.extend(enhancer.find_sequence(m.consensus))

    return hits


# finds the longest width of the motif sequences given
def get_max_motif_width(sequences):
    max_width = 0
    for seq in sequences:
        if (len(seq) > max_width):
            max_width = len(seq)

    return max_width
