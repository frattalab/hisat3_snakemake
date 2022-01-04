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
    
#make sure the output folder exists before running anything
hisat_outdir = get_output_dir(config["project_top_level"], config['histat3n_output_folder'])
os.system("mkdir -p {0}".format(hisat_outdir))

merged_outdir = get_output_dir(config['project_top_level'], config['merged_fastq_folder'])

SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

GENOME_DIR = config['GENOME_DIR']
GENOME_FA = config['fasta']

bedGraph = '/SAN/vyplab/alb_projects/tools/bedGraphToBigWig'

rule all_call:
    input:
        expand(hisat_outdir + "{name}.vcf", name = SAMPLE_NAMES)

MapQ = 20

rule call_samtools_mpileup:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sorted.bam"
    output:
        hisat_outdir + "{name}.vcf"
    params:
        referenceFile = GENOME_FA
    shell:
    """
    samtools mpileup -B -A -f {params.referenceFile} {input} | varscan pileup2snp  --strand-which filter 0 \
    --output-vcf \
    --min-var-freq {params.minVarFreq} \
    --min-coverage 20 \
    --variants 1
    """


# rule snp_mask_bams:
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