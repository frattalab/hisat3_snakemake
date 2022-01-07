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


SAMPLES = pd.read_csv(config["sampleCSVpath"], sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

GENOME_DIR = config['GENOME_DIR']
GENOME_FA = config['fasta']

chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 
'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

rule all_call:
    input:
        expand(hisat_outdir + "{name}.snpmasked.bam", name = SAMPLE_NAMES)

rule call_snp_mask:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        hisat_outdir + "{name}.sorted.bam",
    output:
        dataout = temp(["{{hisat_outdir}}{{name}}.sorted{}.snpmasked.bam".format(i) for i in chr_list]))
    params:
        vcf = config['vcf_path'] #TODO this could be altered to take a VCF per sample in the future
    threads:
        4
    shell:
        """
        echo "Processing {input} file..."
        python3 scripts/snp_mask.py --bam {input} --vcf {params.vcf} --cpu {threads}
    """

rule merge_chromosomes:
    wildcard_constraints:
        sample="|".join(SAMPLE_NAMES)
    input:
        expand(hisat_outdir + "{name}.sorted{chr}.snpmasked.bam", name = SAMPLE_NAMES,chromosome=chr_list)
    output:
        hisat_outdir + "{name}.snpmasked.bam"
    shell:
        """
        samtools merge -o {output} {input}
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