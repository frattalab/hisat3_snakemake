import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

    
#make sure the output folder exists before running anything
hisat_outdir = get_output_dir(config["project_top_level"], config['histat3n_output_folder'])
os.system("mkdir -p {0}".format(hisat_outdir))

merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


GENOME_DIR = "/SAN/vyplab/vyplab_reference_genomes/hisat-3n/mouse/raw"
GENOME_FA = config['fasta']


rule all_split:
    input:
        expand(hisat_outdir + "split_conversions/" +  "{name}.snpmasked.convertedreads.bam", name = SAMPLE_NAMES),
        expand(hisat_outdir + "split_conversions/" +  "{name}.snpmasked.UNconvertedreads.bam", name = SAMPLE_NAMES)


rule split_converted:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.snpmasked.bam"
    output:
        outCon = hisat_outdir + "split_conversions/" +  "{name}.snpmasked.convertedreads.bam",
        outUNCon = hisat_outdir + "split_conversions/" +  "{name}.snpmasked.UNconvertedreads.bam"
    params:
        outputPrefix = os.path.join(hisat_outdir + "split_conversions/"),
        minimum_conversions = 2
    threads:
        4
    shell:
        """
        echo "Processing {input} file..."
        python3 /SAN/vyplab/alb_projects/pipelines/hisat3_snakemake/scripts/split_conversions.py -b {input} \
        -m {params.minimum_conversions} \
        -o {params.outputPrefix}
        """

rule index_split:
    wildcard_constraints:
        name="|".join(SAMPLE_NAMES),
    input:
        outCon = hisat_outdir + "split_conversions/" +  "{name}.snpmasked.convertedreads.bam",
        outUNCon = hisat_outdir + "split_conversions/" +  "{name}.snpmasked.UNconvertedreads.bam"
    output:
        hisat_outdir + "split_conversions/" +  "{name}.snpmasked.convertedreads.bam.bai",
        hisat_outdir + "split_conversions/" +  "{name}.snpmasked.UNconvertedreads.bam.bai"
    shell:
        """
        samtools index {input.outCon}
        samtools index {input.outUNCon}
        """