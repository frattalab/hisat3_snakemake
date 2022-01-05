__author__ = '@aleighbrown'
import os
import pysam
import argparse
import pickle


def snpmask_bams(bampath, vcfpath, outfile_path):
    print("Getting T>C and A>G SNPs from VCF")
    a_g, t_c = get_snp_dictionaries(vcfpath)
    print("Done getting T>C and A>G SNPs from VCF")

    infile = pysam.AlignmentFile(bampath, "rb")
    outfile = pysam.AlignmentFile(outfile_path, "wb", template=infile)

    print("Starting to parse BAM file")
    for i, read in enumerate(infile):
        if i > 100 and i % 1_000_000 == 0:
            print (f'{i} reads read')
        #first we're only interested in reads that have a conversion on them
        if read.has_tag('Yf'):
            conversion_tag = read.get_tag("Yf")
            if conversion_tag >= 1:
                #if they're aligned to the plus strand we care about T > C mutations
                #Note that if a read has no conversions, HISAT3N will assign it a value of YZ = +
                #This isn't a problem for this because we're only trying to adjust reads with "conversion" on it 
                if read.get_tag('YZ') == "+":
                    snp_dictionary = t_c
                elif read.get_tag('YZ') == "-":
                    snp_dictionary = a_g
                #now we're going to use the snp dictionary to see if the read positions have anything mapping
                #in the snp dictionary
                current_chrom = set(snp_dictionary[read.reference_name])
                read_positions = read.get_aligned_pairs(matches_only = True,with_seq = True)
                read_mismatches = set([x[1] + 1 for x in read_positions if x[2].islower()]) #convert to 1 based
                if current_chrom.intersection(read_mismatches):
                    new_yf = conversion_tag
                    for i in range(len(current_chrom.intersection(read_mismatches))):
                        new_yf = new_yf - 1
                    read.set_tag('Yf',new_yf)
                    outfile.write(read)         
            else:
                outfile.write(read)

    outfile.close()
    print("Success!")
    return 0


def build_dictionaries(vcfpath):
    vcf = pysam.VariantFile(vcfpath,"r")
    snps = vcf.fetch()
    a_g_pickle = os.path.splitext(vcfpath)[0] + ".a_g_snps.pickle"
    t_c_pickle = os.path.splitext(vcfpath)[0] + ".t_c_snps.pickle"
    snps = vcf.fetch()
    a_g = {}
    t_c = {}
    snps = vcf.fetch()
    for s in snps:
        #first check that the SNP is relevant for this - e.g. either a A > G or a T > C
        if(s.alleles[0] == 'A' and s.alleles[1] == 'G'):
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
    print("All done writing VCF pickles!!")
    return 0


def get_snp_dictionaries(vcfpath):
    #first check if the pickles exist
    a_g_pickle = os.path.splitext(vcfpath)[0] + ".a_g_snps.pickle"
    t_c_pickle = os.path.splitext(vcfpath)[0] + ".t_c_snps.pickle"
    if os.path.isfile(a_g_pickle) and os.path.getsize(a_g_pickle) > 300 and os.path.isfile(t_c_pickle) and os.path.getsize(t_c_pickle) > 300:
        print('VCF pickles exist and are above 300 Bytes - reading them')
        with open(a_g_pickle, 'rb') as handle:
            a_g = pickle.load(handle)
        with open(t_c_pickle, 'rb') as handle:
            t_c = pickle.load(handle)
    else:
        print("no pickles found - pickleing the relevant SNPs from the VCF")
        build_dictionaries(vcfpath)
        a_g, t_c = get_snp_dictionaries(vcfpath)
    return(a_g, t_c)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-v", "--vcf")
    args = parser.parse_args()

    bampath = args.bam
    vcfpath = args.vcf


    infile = pysam.AlignmentFile(bampath, "rb")

    outfile_path = os.path.splitext(bampath)[0] + ".snpmasked.bam"
    outfile = pysam.AlignmentFile(outfile_path, "wb", template=infile)
    
    result = snpmask_bams(bampath, vcfpath, outfile_path)

    return result


if __name__ == "__main__":
    main()

# reads_parsed = []
# #go through all the SNPs in the vcf file
# for r in records:
#     #first check that the SNP is relevant for this - e.g. either a A > G or a T > C
#     if (s.alleles[0] == 'T' and s.alleles[1] == 'C') | (s.alleles[0] == 'A' and s.alleles[1] == 'G') :
#         coverage = infile.count_coverage(contig = r.chrom, start = r.start,stop = r.stop)
#         #now check if there's any coverage at that location
#         if(sum([x[0] for x in coverage]) > 0 ):
#             #if there's coverage - get all the reads that overlap to that location with fetch
#             overlapping_reads = infile.fetch(r.chrom, r.start, r.stop)
#             for o in overlapping_reads:
#                 #now for the reads that actually overlap that position (e.g. not just spliced reads)
#                 if o.get_overlap(r.start, r.stop) > 0:
#                     #fix the Yf tag - which contains the number of characteristic mismatches
#                     new_yf = o.get_tag('Yf') - 1
#                     outfile.write(o)
#                     reads_parsed.append(o.query_name)


bampath = "testing.bam"
vcfpath = "/SAN/vyplab/vyplab_reference_genomes/wtc11_vcf/wtc11.vcf.gz"


