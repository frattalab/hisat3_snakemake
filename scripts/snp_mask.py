import os
import pysam
import pickle

bampath = "testing.bam"
vcfpath = "/SAN/vyplab/vyplab_reference_genomes/wtc11_vcf/wtc11.vcf.gz"


vcf = pysam.VariantFile(vcfpath,"r")
records = vcf.fetch()

infile = pysam.AlignmentFile(bampath, "rb")

outfile_path = os.path.splitext(bampath)[0] + ".snpmasked.bam"
outfile = pysam.AlignmentFile(outfile_path, "wb", template=infile)


reads_parsed = []
#go through all the SNPs in the vcf file
for r in records:
    #first check that the SNP is relevant for this - e.g. either a A > G or a T > C
    if (s.alleles[0] == 'T' and s.alleles[1] == 'C') | (s.alleles[0] == 'A' and s.alleles[1] == 'G') :
        coverage = infile.count_coverage(contig = r.chrom, start = r.start,stop = r.stop)
        #now check if there's any coverage at that location
        if(sum([x[0] for x in coverage]) > 0 ):
            #if there's coverage - get all the reads that overlap to that location with fetch
            overlapping_reads = infile.fetch(r.chrom, r.start, r.stop)
            for o in overlapping_reads:
                #now for the reads that actually overlap that position (e.g. not just spliced reads)
                if o.get_overlap(r.start, r.stop) > 0:
                    #fix the Yf tag - which contains the number of characteristic mismatches
                    new_yf = o.get_tag('Yf') - 1
                    outfile.write(o)
                    reads_parsed.append(o.query_name)


bampath = "testing.bam"
vcfpath = "/SAN/vyplab/vyplab_reference_genomes/wtc11_vcf/wtc11.vcf.gz"


vcf = pysam.VariantFile(vcfpath,"r")
records = vcf.fetch()

infile = pysam.AlignmentFile(bampath, "rb")
for read in infile:
    #SOLVE first for the condition that the library is stranded (e.g. read_1 matches the direction of the RNA)

    #if they're aligned to the plus strand we care about T > C mutations
    if read.get_tag('YZ') == "+":
        

    elifread.get_tag('YZ') == "-":



    
    outfile.write(o)

outfile.close()







bampath = "testing.bam"
vcfpath = "/SAN/vyplab/vyplab_reference_genomes/wtc11_vcf/wtc11.vcf.gz"


vcf = pysam.VariantFile(vcfpath,"r")
snps = vcf.fetch()


def build_dictionaries(vcfpath):
    vcf = pysam.VariantFile(vcfpath,"r")
    a_g_pickle = os.path.splitext(vcfpath)[0] + ".a_g_snps.pickle"
    t_c_pickle = os.path.splitext(vcfpath)[0] + ".t_c_snps.pickle"
    
    snps = vcf.fetch()

    if (os.path.isfile(a_g_pickle) and os.path.isfile(t_c_pickle)):
        with open(a_g_pickle, 'rb') as handle:
            a_g = pickle.load(handle)
        with open(t_c_pickle, 'rb') as handle:
            t_c = pickle.load(handle)

        return a_g, t_c

    else:
        a_g = {}
        t_c = {}
        for s in snps:
            #first check that the SNP is relevant for this - e.g. either a A > G or a T > C
            if(s.alleles[0] == 'A' and s.alleles[1] == 'G'):
                break
                #if the chromosome is already in the dictionary
                if s.chrom in a_g.keys():
                    a_g[s.chrom].append(s.pos)
                else:
                    a_g[s.chrom] = [s.pos]

            elif(s.alleles[0] == 'T' and s.alleles[1] == 'C'):
                #if the chromosome is already in the dictionary
                if s.chrom in t_c.keys():
                    t_c[s.chrom].append(s.pos)
                else:
                    t_c[s.chrom] = [s.pos]
            else:
                continue

        with open(a_g_pickle, 'wb') as handle:
            pickle.dump(a_g, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(t_c_pickle, 'wb') as handle:
            pickle.dump(t_c, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return a_g, t_c


