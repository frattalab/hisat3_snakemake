__author__ = '@aleighbrown'
import os
import pysam
import argparse
import pickle
from multiprocessing import Pool

def collect_stats(infile,annotated_sjintervals,cryptic_sjintervals):


    for read in infile:
    
        if read.is_unmapped:
            next
        total_reads += 1
        split_read = False
        intronic_region = []
        
        if read.is_reverse:
            conversions_on = 'A'
        else:
            conversions_on = 'T'
        

        if 'N' in read.cigarstring:
            split_read = True
            spliced_reads += 1

        #get all the "splice junctions" in a read
        if split_read:

            last_read_pos = False
            index = 0
            for read_loc, genome_loc in read.get_aligned_pairs():

                # get_aligned_pairs returns None for the read location when there is a gap(splice)
                #so this sets the start of the splice junction to the genomic location if the read loc is none
                if read_loc is None and last_read_pos:
                    start = genome_loc

                # this is to deal with cases where the SJ overhang is only 1 bp and its on the first base 
                if index == 1 and read_loc is None:

                    start = genome_loc

                elif read_loc and last_read_pos is None:

                    stop = genome_loc  # we are right exclusive ,so this is correct

                    intronic_region.append((start, stop))

                    del start
                    del stop

                last_read_pos = read_loc
                index +=1

            if set(intronic_region).intersection(annotated_sjintervals):
                spliced_annotated += 1
                t_coverage_annotated += read.get_reference_sequence().upper().count(conversions_on)

                conversion_tag = read.get_tag("Yf")

                if conversion_tag > 0:
                    spliced_annotated_coverted += 1


            elif set(intronic_region).intersection(cryptic_sjintervals):
                spliced_cryptic += 1
                t_coverage_cryptic += read.get_reference_sequence().upper().count(conversions_on)

                conversion_tag = read.get_tag("Yf")

                if conversion_tag > 0:
                    spliced_cryptic_coverted += 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-r", "--regions")


    args = parser.parse_args()

    bampath = args.bam
    vcfpath = args.vcf
    cpu = int(args.cpu)



    return 0


if __name__ == "__main__":
    main()

