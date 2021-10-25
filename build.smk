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
os.system("mkdir -p {0}".format(hisat_outdir))

merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()


GENOME_DIR = "/SAN/vyplab/vyplab_reference_genomes/hisat-3n/human/raw"
GENOME_FA = config['fasta']
CHRMSIZES = config['fasta_sizes']
bedGraph = '/SAN/vyplab/alb_projects/tools/bedGraphToBigWig'

rule all_hisat3n:
    input:
        expand(hisat_outdir + "{name}.sorted.sam", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.conversion.tsv", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.count.+.bw"),
        expand(hisat_outdir + "{name}.count.-.bw"),
        expand(hisat_outdir + "{name}.rate.+.bw"),
        expand(hisat_outdir + "{name}.rate.-.bw")

rule run_histat3n_pe:
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
        --base-change {params.baseChange}\
        --threads {threads}
        """
rule sort_histat:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sam"
    output:
        hisat_outdir + "{name}.sorted.sam"
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
rule split_conversion:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.conversion.tsv"
    output:
        temp(hisat_outdir + "{name}.+.txt"),
        temp(hisat_outdir + "{name}.-.txt")
    params:
        outputPrefix = os.path.join(hisat_outdir + "{name}."),
    shell:
        """
        awk -F'\t' '{ print > "{params.name}" $3 ".txt" }' {input}
        """

rule split_toBedGraph:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plus = hisat_outdir + "{name}.+.txt",
        minus = hisat_outdir + "{name}.-.txt"
    output:
        plusC = temp(hisat_outdir + "{name}.+.count.bedgraph"),
        minusC = temp(hisat_outdir + "{name}.-.count.bedgraph"),
        plusR = temp(hisat_outdir + "{name}.+.rate.bedgraph"),
        minusR = temp(hisat_outdir + "{name}.-.rate.bedgraph")
    shell:
        """
        awk  -F "\t" '{print $1 "\t" $2 "\t" $2 + 1 "\t" $5}' {input.plus} > {output.plusC}
        awk  -F "\t" '{print $1 "\t" $2 "\t" $2 + 1"\t" $5}' {input.minus} > {output.minusC}
        awk  -F "\t" '{if($7+$5 ==0) $4 = 0; else $4 = ($5)/($7+$5)} {print $1 "\t" $2 "\t" $2 + 1"\t" $4}' {input.plus} > {output.plusR}
        awk  -F "\t" '{if($7+$5 ==0) $4 = 0; else $4 = ($5)/($7+$5)} {print $1 "\t" $2 "\t" $2 + 1"\t" $4}' {input.minus} > {output.minusR}
        """

rule sortBedGraph:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusC = hisat_outdir + "{name}.+.count.bedgraph",
        minusC = hisat_outdir + "{name}.-.count.bedgraph",
        plusR = hisat_outdir + "{name}.+.rate.bedgraph",
        minusR = hisat_outdir + "{name}.-.rate.bedgraph"
    output:
        plusC = temp(hisat_outdir + "{name}.+.count.sorted.bedgraph"),
        minusC = temp(hisat_outdir + "{name}.-.count.sorted.bedgraph"),
        plusR = temp(hisat_outdir + "{name}.+.rate.sorted.bedgraph"),
        minusR = temp(hisat_outdir + "{name}.-.rate.sorted.bedgraph")
    shell:
        """
        LC_COLLATE=C sort -k1,1 -k2,2n {input.plusC} > {output.plusC}
        LC_COLLATE=C sort -k1,1 -k2,2n {input.minusC} > {output.minusC}
        LC_COLLATE=C sort -k1,1 -k2,2n {input.plusR} > {output.plusR}
        LC_COLLATE=C sort -k1,1 -k2,2n {input.minusR} > {output.minusR}
        """

rule bedGraphtoBW:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusC = hisat_outdir + "{name}.+.count.sorted.bedgraph",
        minusC = hisat_outdir + "{name}.-.count.sorted.bedgraph",
        plusR = hisat_outdir + "{name}.+.rate.sorted.bedgraph",
        minusR = hisat_outdir + "{name}.-.rate.sorted.bedgraph"
    output:
        plusC = hisat_outdir + "{name}.count.+.bw",
        minusC = hisat_outdir + "{name}.count.-.bw",
        plusR = hisat_outdir + "{name}.rate.+.bw",
        minusR = hisat_outdir + "{name}.rate.-.bw"
    shell:
        """
        /SAN/vyplab/alb_projects/tools/bedGraphToBigWig {input.plusC} \
        {CHRMSIZES} {output.plusC}
        /SAN/vyplab/alb_projects/tools/bedGraphToBigWig {input.minusC} \
        {CHRMSIZES} {output.minusC}
        /SAN/vyplab/alb_projects/tools/bedGraphToBigWig {input.plusR} \
        {CHRMSIZES} {output.plusR}
        /SAN/vyplab/alb_projects/tools/bedGraphToBigWig {input.minusR} \
        {CHRMSIZES} {output.minusR}
        """