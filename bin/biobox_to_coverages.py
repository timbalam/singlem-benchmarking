#!/usr/bin/env python3
from cami_opal.utils import load_data
import sys
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.condense import CondensedCommunityProfile

def read_biobox(input_biobox_path):
    levels = ['kingdom','phylum','class','order','family','genus','species']
    total_coverages = [0.0]*len(levels)

    tables = []
    sample_list = load_data(input_biobox_path, False)
    for sample_id, header, profile in sample_list:

        sample_summary_root_node = WordNode(None, "Root")
        for prediction in profile:
            sample_summary_root_node.add_words(prediction.taxpath.split('|'), prediction.percentage)
            
        condensed_table = CondensedCommunityProfile(sample_id, sample_summary_root_node)
        # Pass 2 - calculate percents
        for wordnode in condensed_table.breadth_first_iter():
            # Percent includes unassigned
            unassigned_cov = wordnode.coverage - wordnode.get_full_coverage()
            wordnode.coverage = unassigned_cov


        tables.append(condensed_table)

    return tables