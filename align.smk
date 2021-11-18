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

rule all_hisat3n:
    input:
        expand(hisat_outdir + "{name}.conversion.tsv", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.sorted.bam.bai", name = SAMPLE_NAMES),
        # expand(hisat_outdir + "{name}.sorted.tagged.bam", name = SAMPLE_NAMES),
        # expand(hisat_outdir + "{name}.sorted.tagged.bam.bai", name = SAMPLE_NAMES)




rule run_hisat3_pe:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        generated_index = GENOME_DIR + ".3n.CT.1.ht2",
        one = lambda wildcards: get_processed_fastq(wildcards.name, pair=1),
        two = lambda wildcards: get_processed_fastq(wildcards.name, pair=2)
    output:
        temp(hisat_outdir + "{name}.sam")
    params:
        genomeDir = GENOME_DIR,
        outputPrefix = os.path.join(hisat_outdir + "{name}.sam"),
        strandedness = config['strandedness'],
        baseChange = "T,C"
    threads:
        4
    shell:
        """
        /SAN/vyplab/alb_projects/tools/hisat-3n/hisat-3n \
        -x {params.genomeDir} \
        -1 {input.one} \
        -2 {input.two} \
        -q \
        -S {params.outputPrefix} \
        --base-change {params.baseChange} \
        --threads {threads} \
        --rna-strandness {params.strandedness}
        """

rule run_hisat3_se:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        generated_index = GENOME_DIR + ".3n.CT.1.ht2",
        one = lambda wildcards: get_processed_fastq(wildcards.name, pair=1)
    output:
        temp(hisat_outdir + "{name}.sam")
    params:
        genomeDir = GENOME_DIR,
        outputPrefix = os.path.join(hisat_outdir + "{name}.sam"),
        strandedness = config['strandedness'],
        baseChange = "T,C"
    threads:
        4
    shell:
        """
        echo "This is our memory amount"
        free -mh
        echo "And the nproc"
        nproc
        /SAN/vyplab/alb_projects/tools/hisat-3n/hisat-3n \
        -x {params.genomeDir} \
        -U {input.one} \
        -q \
        -S {params.outputPrefix} \
        --base-change {params.baseChange} \
        --threads {threads} \
        --rna-strandness {params.strandedness}
    """

rule sort_histat:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sam"
    output:
        temp(hisat_outdir + "{name}.sorted.sam")
    threads:
        4
    shell:
        """
        t=/scratch0/$USER/$RANDOM
        mkdir -p $t
        samtools sort {input} -o {output} -T $t
        """
rule conversion_table:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sorted.sam"
    output:
        hisat_outdir + "{name}.conversion.tsv"
    params:
        genomeFA = GENOME_FA,
        outputPrefix = os.path.join(hisat_outdir + "{name}.conversion.tsv"),
        baseChange = "T,C"
    threads:
        4
    shell:
        """
        /SAN/vyplab/alb_projects/tools/hisat-3n/hisat-3n-table \
        --alignments {input} \
        --ref {params.genomeFA} \
        --output-name {params.outputPrefix} \
        --base-change {params.baseChange} \
        --threads {threads}
        """
rule sort_bams:
    input:
        hisat_outdir + "{name}.sorted.sam"
    output:
        hisat_outdir + "{name}.sorted.bam"
    shell:
        """
        samtools view -S --threads 4 -b {input} > {output} 
        """

rule index_bams:
    input:
        hisat_outdir + "{name}.sorted.bam"
    output:
        hisat_outdir + "{name}.sorted.bam.bai"
    threads:
        4
    shell:
        """
        samtools index {input} 
        """

# rule tag_bams:
#     input:
#         bam = hisat_outdir + "{name}.sorted.bam",
#         bai = hisat_outdir + "{name}.sorted.bam.bai"
#     output:
#         hisat_outdir + "{name}.sorted.tagged.bam"
#     params:
#         pickled = GENOME_FA + '.pickle'
#     shell:
#         """
#         python3 scripts/slamdunk_taggers.py -b {input.bam} -p {params.pickled}
#         """
# rule index_tagged_bams:
#     input:
#         hisat_outdir + "{name}.sorted.tagged.bam"
#     output:
#         hisat_outdir + "{name}.sorted.tagged.bam.bai"
#     threads:
#         4
#     shell:
#         """
#         samtools index {input} 
#         """