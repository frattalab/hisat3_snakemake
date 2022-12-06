import os
from pathlib import Path

INPUT_DIR="/SAN/vyplab/alb_projects/data/4su_better_design_ward_ipsc/new_analysis/HISAT3N/"
OUTPUT_DIR="/SAN/vyplab/alb_projects/data/4su_better_design_ward_ipsc/new_analysis/spliced_nomask/"
ANNOTATED_JUNCTIONS="/SAN/vyplab/alb_projects/data/4su_full_ward_tdp_kd_ipsc/controlHumphreyCorticalNeuron-TDP43KDHumphreyCorticalNeuron_annotated_junctionscryptic_clusters.bed"
bam_suffix = '.sorted.bam'

basenameBed = Path(ANNOTATED_JUNCTIONS).stem

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(INPUT_DIR) if f.endswith(bam_suffix)]

rule all:
    input:
        expand(OUTPUT_DIR + "{sample}" + "_" + basenameBed + "_spliced_counts.csv", sample = SAMPLES)


rule calculate_splice_stability:
    input:
        bamfile = INPUT_DIR + "{sample}" + bam_suffix
    output:
        outputfile = OUTPUT_DIR + "{sample}" + "_" + basenameBed + "_spliced_counts.csv"
    shell:
        """
        echo {input.bamfile}
        python3 scripts/calculate_splice_stability.py -b {input.bamfile} -r {ANNOTATED_JUNCTIONS} -o {output.outputfile}
        """


# expected TDP43kd_1_8h_controlHumphreyCorticalNeuron-TDP43KDHumphreyCorticalNeuron_annotated_junctionscryptic_clusters_spliced_counts.csv
# output   TDP43kd_1_8h.snpmasked_controlHumphreyCorticalNeuron-TDP43KDHumphreyCorticalNeuron_annotated_junctionscryptic_clusters_spliced_counts.csv