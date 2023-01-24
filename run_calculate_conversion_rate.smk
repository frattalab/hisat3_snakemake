import os
from pathlib import Path

INPUT_DIR="/SAN/vyplab/alb_projects/data/4su_full_ward_tdp_kd_ipsc/HISAT3N/"
OUTPUT_DIR="/SAN/vyplab/alb_projects/data/4su_full_ward_tdp_kd_ipsc/HISAT3N_perbase_conv/"
INPUT_BED="/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40_3utr_unique.bed"
conversion_suffix = '.conversion.fake.bed.gz'

basenameBed = Path(INPUT_BED).stem

SAMPLES = [f.replace(conversion_suffix, "") for f in os.listdir(INPUT_DIR) if f.endswith(conversion_suffix)]

rule all:
    input:
        expand(OUTPUT_DIR + "{sample}" + "_" + basenameBed + "_perbase_cov.tsv", sample = SAMPLES)


rule calculate_splice_stability:
    input:
        conversion_file = INPUT_DIR + "{sample}" + conversion_suffix
    output:
        outputfile = OUTPUT_DIR + "{sample}" + "_" + basenameBed + "_perbase_cov.tsv"
    params:
        bed = INPUT_BED
    shell:
        """
        tabix {input.conversion_file}\
        -R {params.bed} | mawk -F '\t' '{{sum6+=$6; sum8+=$8}END{print sum6,sum8}}' > {output.outputfile}
        """
