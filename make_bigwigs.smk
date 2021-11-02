
import os
configfile: "config/config.yaml"
cluster_config: "config/cluster.yaml"
include: "helpers.py"

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

rule all_makeBW:
    input:
        expand(hisat_outdir + "{name}.count.bw", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.rate.bw", name = SAMPLE_NAMES)

rule filter_conversion:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.conversion.tsv"
    output:
        temp(hisat_outdir + "{name}.filtered.conversion.tsv")
    params:
        outputPrefix = os.path.join(hisat_outdir + "{name}."),
        grep_pattern = "chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY|chrM"
    shell:
        """
        grep -E '{params.grep_pattern}' {input} > {output}
        """

rule split_toBedGraphRate:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plus = hisat_outdir + "{name}.filtered.conversion.tsv",
    output:
        plusR = temp(hisat_outdir + "{name}.rate.bedgraph")
    shell:
        """
        source scripts/rateAwk.sh {input.plus} {output.plusR}
        """

rule split_toBedGraphCount:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plus = hisat_outdir + "{name}.filtered.conversion.tsv"
    output:
        plusC = temp(hisat_outdir + "{name}.count.bedgraph")
    shell:
        """
        source scripts/countAwk.sh {input.plus} {output.plusC}
        """


rule sortBedGraphPlusCount:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusC = hisat_outdir + "{name}.count.bedgraph"
    output:
        plusC = temp(hisat_outdir + "{name}.count.sorted.bedgraph")
    threads:
        2
    shell:
        """
        LC_COLLATE=C sort -k1,1 -k2,2n {input.plusC} > {output.plusC}
        """

rule sortBedGraphPlusRate:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusR = hisat_outdir + "{name}.rate.bedgraph"
    output:
        plusR = temp(hisat_outdir + "{name}.rate.sorted.bedgraph")
    threads:
        2
    shell:
        """
        LC_COLLATE=C sort --parallel 2 -k1,1 -k2,2n {input.plusR} > {output.plusR}
        """

rule bedGraphtoBWPlusCount:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusC = hisat_outdir + "{name}.count.sorted.bedgraph"
    output:
        plusC = hisat_outdir + "{name}.count.bw"
    params:
        bedGraphToolPath = config['bedGraphToolPath']
    shell:
        """
        {params.bedGraphToolPath} {input.plusC} \
        {CHRMSIZES} {output.plusC}
        """

rule bedGraphtoBWPlusRate:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusR = hisat_outdir + "{name}.rate.sorted.bedgraph"
    output:
        plusR = hisat_outdir + "{name}.rate.bw"
    shell:
        """
        /SAN/vyplab/alb_projects/tools/bedGraphToBigWig {input.plusR} \
        {CHRMSIZES} {output.plusR}
        """

