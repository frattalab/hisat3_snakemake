import os
import subprocess
#mkdir -p {params.filter_output_dir}
# include: "helpers.py"
configfile: "config/config.yaml"

project_folder =  "/SAN/vyplab/alb_projects/data/4su_better_design_ward_ipsc/new_analysis/"
stability_output = "/SAN/vyplab/alb_projects/data/4su_better_design_ward_ipsc/new_analysis/spliced_read_stability_HISAT3N_snpmasked/"
bam_suffix = '.snpmasked.bam.bai'

bam_dir = '/SAN/vyplab/alb_projects/data/4su_better_design_ward_ipsc/new_analysis/'

SAMPLE_NAMES = [f.replace(bam_suffix, "") for f in os.listdir(bam_dir) if f.endswith(bam_suffix)]

bed_file = "/SAN/vyplab/alb_projects/pipelines/hisat3_snakemake/test_stability_input.bed"
bed_name  = os.path.splitext(os.path.basename(bed_file))[0]

print(SAMPLE_NAMES)

rule all:
    input:
        expand(stability_output + ".snpmasked_test_stability_input_spliced_counts.csv", sample = SAMPLE_NAMES)


rule stability_splice_count:
    input:
        input_bam = bam_dir + "{sample}.snpmasked.bam"
        #input_vcf = os.path.join(slamdunk_output + "snp/") + "{sample}_mapped_filtered_snp.vcf"
    output:
        stability_output + "{sample}.snpmasked_test_stability_input_spliced_counts.csv"
    params:
        bed = bed_file
    threads:
        8
    shell:
        """
        python3 /SAN/vyplab/alb_projects/pipelines/hisat3_snakemake/scripts/calculate_splice_stability.py \
        -b {input} \
        -r {params.bed} \
        -o /SAN/vyplab/alb_projects/pipelines/hisat3_snakemake -t 8
        """
