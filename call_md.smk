import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

#make sure the output folder exists before running anything
hisat_outdir = get_output_dir(config["project_top_level"], config['histat3n_output_folder'])
mdcalleddir = get_output_dir(config["project_top_level"], config['mdcalleddir'])
os.system("mkdir -p {0}".format(mdcalleddir))


SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

GENOME_DIR = config['GENOME_DIR']
GENOME_FA = config['fasta']



rule all_call:
    input:
        expand(mdcalleddir + "{name}.mdcalled.bam", name = SAMPLE_NAMES),
        expand(mdcalleddir + "{name}.mdcalled.bam.bai", name = SAMPLE_NAMES)


rule call_md_tag:
    wildcard_constraints:
        name="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sorted.bam"
    output:
        mdcalleddir + "{name}.mdcalled.bam"
    threads:
        4
    shell:
        """
        echo "Processing {input} file..."
        samtools calmd -@ 4 -b {input} {GENOME_FA} > {output}
    """

rule index_mdcalled:
    wildcard_constraints:
        name="|".join(SAMPLE_NAMES),
    input:
        mdcalleddir + "{name}.mdcalled.bam"
    output:
        mdcalleddir + "{name}.mdcalled.bam.bai"
    shell:
        """
        samtools index {input}
        """
