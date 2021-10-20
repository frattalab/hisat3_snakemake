import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

# RULE ORDER DIRECTIVE
# if paired end, use the paired end rule to run, if single end use the single end rule to run
# if config['end_type'] == "pe":
#     ruleorder: run_star_pe > run_star_se
# else:
#     ruleorder: run_star_se > run_star_pe


#make sure the output folder for STAR exists before running anything
hisat_outdir = get_output_dir(config["project_top_level"], config['histat3n_output_folder'])
os.system("mkdir -p {0}".format(star_outdir))

merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


GENOME_DIR = "/SAN/vyplab/vyplab_reference_genomes/hisat-3n/human/raw"

rule all_hisat3n:
    input:
        expand(hisat_outdir + "{name}.sam",name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.sorted.sam", name = SAMPLE_NAMES)

rule run_histat3n_pe:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        generated_index = GENOME_DIR,
        one = lambda wildcards: get_processed_fastq(wildcards.name, pair=1),
        two = lambda wildcards: get_processed_fastq(wildcards.name, pair=2)
    output:
        expand(hisat_outdir + "{name}.sam",name = SAMPLE_NAMES)
    params:
        genomeDir = GENOME_DIR,
        outTmpDir = os.path.join(star_outdir + "{name}_tmpdir"),
        outputPrefix = os.path.join(hisat_outdir + "{name}.sam"),
        baseChange = "T,C"
    threads:
        4
    shell:
        """
        /SAN/vyplab/alb_projects/tools/hisat-3n/hisat-3n -x \
        {params.genomeDir} \
        -1 {input.one} \
        -2 {input.two} \
        -q \
        -S {params.outputPrefix}
        --base-change {params.baseChange}\
        --threads {threads}
        """