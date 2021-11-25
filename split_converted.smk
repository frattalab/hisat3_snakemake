import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"
# RULE ORDER DIRECTIVE
# if paired end, use the paired end rule to run, if single end use the single end rule to run
if config['end_type'] == "pe":
	ruleorder: run_hisat3_pe > run_hisat3_se
else:
	ruleorder: run_hisat3_se > run_hisat3_pe
    
#make sure the output folder for STAR exists before running anything
hisat_outdir = get_output_dir(config["project_top_level"], config['histat3n_output_folder'])
os.system("mkdir -p {0}".format(hisat_outdir))

merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


GENOME_DIR = "/SAN/vyplab/vyplab_reference_genomes/hisat-3n/mouse/raw"
GENOME_FA = config['fasta']

bedGraph = '/SAN/vyplab/alb_projects/tools/bedGraphToBigWig'

rule all_split:
    input:
        expand(hisat_outdir + "split_conversions/" +  "{name}.sorted.convertedreads.bam", name = SAMPLE_NAMES),
        expand(hisat_outdir + "split_conversions/" +  "{name}.sorted.UNconvertedreads.bam", name = SAMPLE_NAMES)


rule split_converted:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sorted.bam"
    output:
        outCon = expand(hisat_outdir + "split_conversions/" +  "{name}.sorted.convertedreads.bam", name = SAMPLE_NAMES),
        outUNCon = expand(hisat_outdir + "split_conversions/" +  "{name}.sorted.UNconvertedreads.bam", name = SAMPLE_NAMES)
    params:
        outputPrefixCon = os.path.join(hisat_outdir + "{name}.sorted.convertedreads.bam")
        outputPrefixUNCon = os.path.join(hisat_outdir + "{name}.sorted.UNconvertedreads.bam")
    threads:
        4
    shell:
        """
        echo "Processing {input} file..."
        python3 /SAN/vyplab/alb_projects/pipelines/hisat3_snakemake/split_conversions.py -b {input} -m 2
        mv {params.outputPrefixCon} {output.outCon}
        mv {params.outputPrefixUNCon} {output.outUNCon}
        """
