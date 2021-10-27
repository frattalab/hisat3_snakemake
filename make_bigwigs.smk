
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
        expand(hisat_outdir + "{name}.count.plus.bw", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.count.minus.bw", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.rate.plus.bw", name = SAMPLE_NAMES),
        expand(hisat_outdir + "{name}.rate.minus.bw",name = SAMPLE_NAMES)

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
        awk -F"\\t" "{{ print > "{params.outputPrefix}" $3 ".txt" }}" {input}
        """


rule split_toBedGraph:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plus = hisat_outdir + "{name}.+.txt",
        minus = hisat_outdir + "{name}.-.txt"
    output:
        plusC = temp(hisat_outdir + "{name}.plus.count.bedgraph"),
        minusC = temp(hisat_outdir + "{name}.minus.count.bedgraph"),
        plusR = temp(hisat_outdir + "{name}.plus.rate.bedgraph"),
        minusR = temp(hisat_outdir + "{name}.minus.rate.bedgraph")
    shell:
        """
        countAwk.sh {input.plus} {output.plusC}
        countAwk.sh {input.minus} {output.minusC}
        rateAwk.sh {input.plus} {output.plusR}
        rateAwk.sh {input.minus} {output.minusR}
        """

rule sortBedGraph:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        plusC = hisat_outdir + "{name}.plus.count.bedgraph",
        minusC = hisat_outdir + "{name}.minus.count.bedgraph",
        plusR = hisat_outdir + "{name}.plus.rate.bedgraph",
        minusR = hisat_outdir + "{name}.minus.rate.bedgraph"
    output:
        plusC = temp(hisat_outdir + "{name}.plus.count.sorted.bedgraph"),
        minusC = temp(hisat_outdir + "{name}.minus.count.sorted.bedgraph"),
        plusR = temp(hisat_outdir + "{name}.plus.rate.sorted.bedgraph"),
        minusR = temp(hisat_outdir + "{name}.minus.rate.sorted.bedgraph")
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
        plusC = hisat_outdir + "{name}.plus.count.sorted.bedgraph",
        minusC = hisat_outdir + "{name}.minus.count.sorted.bedgraph",
        plusR = hisat_outdir + "{name}.plus.rate.sorted.bedgraph",
        minusR = hisat_outdir + "{name}.minus.rate.sorted.bedgraph"
    output:
        plusC = hisat_outdir + "{name}.count.plus.bw",
        minusC = hisat_outdir + "{name}.count.minus.bw",
        plusR = hisat_outdir + "{name}.rate.plus.bw",
        minusR = hisat_outdir + "{name}.rate.minus.bw"
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