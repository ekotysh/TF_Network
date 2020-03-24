from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import logging.config

# LOGGER
logging.config.fileConfig('logging.conf')
logger = logging.getLogger(__name__)


# CLASSES
class Enhancer:
    chr = ""
    start = -1
    stop = -1
    elite = False
    target_gene_ensemble_id = ""
    target_gene = None
    sequence = ""

    def __init__(self, chrom, start, stop, ensemble_id, elite):
        self.chr = chrom
        self.start = start
        self.stop = stop
        self.target_gene_ensemble_id = ensemble_id
        self.elite = elite

    def find_sequence(self, sequence_to_find):
        sequence_obj = Seq(self.sequence, generic_dna)
        reverse_comp = sequence_obj.reverse_complement()

        matches = []
        # search for exact matches of the consensus sequence
        pos = sequence_obj.find(sequence_to_find)
        while (pos >= 0):
            # got a hit, record position and matching sequence
            matches.append((pos, sequence_to_find, '+'))
            pos = sequence_obj.find(sequence_to_find, pos+1)

        # also, look for matches on the complementary strand
        pos = reverse_comp.find(sequence_to_find)
        while (pos >= 0):
            # got a hit, record position and sequence
            matches.append((pos, sequence_to_find, '-'))
            pos = reverse_comp.find(sequence_to_find, pos + 1)

        return matches

    @staticmethod
    def parse_from_bed_file(column_list, elite):
        chrom = column_list[0]
        start = column_list[1]
        stop = column_list[2]
        target_gene_ensemble_id = column_list[3]
        elite = elite
        return Enhancer(chrom, (int)(start), (int)(stop),
                        target_gene_ensemble_id, elite)

    @staticmethod
    def make_unique_key(enh_chr, enh_start, enh_stop):
        return ">" + enh_chr + ":" + str(enh_start) + "-" + str(enh_stop)


###############################################################
# Functions
###############################################################

def extract_enhancers(bed_file, gene_names_dict):
    # open the enhancer data file for reading
    gene_hancer_data = open(bed_file, "r")

    # initialize enhancer dictionary
    enhancer_dict = {}

    # Read GeneHancer bed file, line by line
    logger.debug("Loading enhancers")
    for enhancer_row in gene_hancer_data:
        row_values = enhancer_row.strip().split()  # split line into list at whitespace
        enhancer = Enhancer.parse_from_bed_file(row_values, True)  # create an elite enhancer obj from row

        # Look up the Hugo name of the target gene and its biotype/description in genes.ENSG.tbl
        gene = gene_names_dict.get(enhancer.target_gene_ensemble_id)
        enhancer.target_gene = gene

        # add enhancer to the map
        key = Enhancer.make_unique_key(enhancer.chr, enhancer.start, enhancer.stop)
        enhancer_dict[key] = enhancer

    # close the enhancer file
    gene_hancer_data.close()
    return enhancer_dict


def populate_enhancer_sequences(enhancer_map, sequences_file):
    # Sequences file is in this format:
    # >chr1:11369-12510
    # sequence
    # >chr1:12369-13510
    # sequence

    logger.debug("Populating enhancer sequences")

    enh_sequences_data = open(sequences_file, "r")
    while True:
        header = enh_sequences_data.readline().strip()
        sequence = enh_sequences_data.readline().strip()
        if not sequence:
            break

        # find this enhancer in our map
        enh = enhancer_map.get(header)
        if (enh is None):
            logger.error("Failed to find an enhnacer for " + header)
        else:
            # assign the sequence to it
            enh.sequence = sequence

    enh_sequences_data.close()

