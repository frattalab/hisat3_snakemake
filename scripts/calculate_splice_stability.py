__author__ = '@aleighbrown'

import pyranges as pr
import pandas as pd
from pathlib import Path
import pysam
import sys
import argparse
import os
def add_key(dictionary, key, value):
    if key not in dictionary.keys():
        dictionary[key] = [value]
    else:
        dictionary[key].append(value)
    return dictionary

def bed_to_odd_dict(bedpath):
    infile = open(bedpath, 'r')

    #create an empty dictionary
    bed_dict = {}

    #loop over every line in the file
    for line in infile:
        #split the line into its parts
        fields = line.strip().split()
        #extract the chromosome, start and end values
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        #add the chromosome as the key, and the start and end values as the values
        add_key(bed_dict,chrom, (start, end))

    #close the file
    infile.close()
    nested_dict = {key:{value:[0,0,0,0] for value in bed_dict[key]} for key in bed_dict} 
    
    return nested_dict
def splice_junctions(read):
    """
    This function takes in a read (SAM/BAM format) and returns a list of 
    tuples containing the starting and ending positions of all of the 
    splice junctions in the read.
    """
    current_read_pos = read.reference_start 
    splice_starts = list()
    splice_ends = list()
    if read.cigartuples:
        for i in range(len(read.cigartuples)):
            if read.cigartuples[i][0] == 0:
                current_read_pos = current_read_pos + read.cigartuples[i][1]
            if read.cigartuples[i][0] == 3:
                splice_starts.append(current_read_pos)
                current_read_pos = current_read_pos + read.cigartuples[i][1]
                splice_ends.append(current_read_pos)
        splice_starts = [s - 1 for s in splice_starts]
        introns = list(zip(splice_starts,splice_ends))

        return introns

def bed_to_dict(bedpath):
    infile = open(bedpath, 'r')

    #create an empty dictionary
    bed_dict = {}

    #loop over every line in the file
    for line in infile:
        #split the line into its parts
        fields = line.strip().split()
        #extract the chromosome, start and end values
        chrom = fields[0]
        start = fields[1]
        end = fields[2]
        #add the chromosome as the key, and the start and end values as the values
        bed_dict[chrom] = (start, end)

    #close the file
    infile.close()
    

def get_tag(read):
    spliced_coverted = 0
    if read.is_reverse:
        conversions_on = 'A'
    else:
        conversions_on = 'T'

    t_coverage = read.get_reference_sequence().upper().count(conversions_on)

    conversion_tag = read.get_tag("Yf")

    if conversion_tag > 0:
        spliced_coverted = 1
    else:
        spliced_converted = 0


    return spliced_coverted, conversion_tag,t_coverage
    
def process_bams(bam_file_path,bed_file_path):
    
    bam = pysam.AlignmentFile(bam_file_path, "rb")
    the_bed = pr.read_bed(bed_file_path)
    
    search_range = the_bed.merge()
    storage_dict = bed_to_odd_dict(bed_file_path)
    
    for row in search_range.df.itertuples():
        set1 = set(storage_dict[row.Chromosome].keys())

        for read in  bam.fetch(row.Chromosome,row.Start,row.End):
            spliced_converted = 0
            conversion_tag = 0
            t_coverage = 0 
            introns = splice_junctions(read)
            if introns:
                set2 = set(introns)
                overlap = set1.intersection(set2)
                if len(overlap) > 0:
                    spliced_converted, conversion_tag, t_coverage = get_tag(read)
                    for i in introns:
                        if i in storage_dict[row.Chromosome].keys():
                            storage_dict[row.Chromosome][i][0] += 1
                            storage_dict[row.Chromosome][i][1] += spliced_converted
                            storage_dict[row.Chromosome][i][2] += conversion_tag
                            storage_dict[row.Chromosome][i][3] += t_coverage

    data = [(k1, k2[0], k2[1], v[0], v[1], v[2], v[3]) for k1 in storage_dict for k2, v in storage_dict[k1].items()]

    df = pd.DataFrame(data, columns=['Chromosome', 'Start', 'End', 'n_spliced_reads', 'n_converted_spliced', 'n_conversions', 'n_possible_conversions'])
    merged_df = pd.merge(the_bed.df, df, on=['Chromosome','Start','End'], how='left')
                
    return merged_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-r", "--regions")
    parser.add_argument("-o", "--outputfolder")


    args = parser.parse_args()

    bampath = args.bam
    bedpath = args.regions
    outfolder = args.outputfolder



    basenameBam = Path(bampath).stem
    basenameBed = Path(bedpath).stem
    
    counts = process_bams(bampath,bedpath)

    outputfile = os.path.join(outfolder, basenameBam + "_" + basenameBed + "_" + "spliced_counts.csv")
    

    counts.to_csv(outputfile,index = False)



    return 0

if __name__ == "__main__":
    main()

