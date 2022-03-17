__author__ = '@aleighbrown'
import os
import pysam
import argparse
from pathlib import Path
import pyranges
import pandas as pd


def collect_stats(infile,sjintervals):
    spliced = 0
    t_coverage = 0
    spliced_coverted = 0
    
    

    for read in infile.fetch(sj_info[0]):
        
        split_read = False
        intronic_region = []
        
        if read.is_unmapped:
            next
        else:
            
            if read.is_reverse:
                conversions_on = 'A'
            else:
                conversions_on = 'T'
            
            try:
                if 'N' in read.cigarstring:
                    split_read = True
            except:
                print("somehow I'm here")
                if 'N' in read.cigarstring:
                    split_read = True
                    print('now i did that')

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


                if set(intronic_region).intersection(sjintervals):

                    spliced += 1

                    t_coverage += read.get_reference_sequence().upper().count(conversions_on)

                    conversion_tag = read.get_tag("Yf")

                    if conversion_tag > 0:
                        spliced_coverted += 1

    rows = [spliced,t_coverage,spliced_coverted]

    return rows

def process_bam(infile,regions,chr):
    rows = []
    junctions = pyranges.readers.read_bed(regions, as_df=True, nrows=None)
    
    print("Hi I'm about to start reading in the file...")
    samfile = pysam.AlignmentFile(infile, "rb")
    print('Processing Junctions:',flush=True)
    for index, row in junctions.iterrows():
        print(row['Name'],flush=True)
        jnc = [row['Chromosome'], row['Start'],row['End']]
        sjintervals = set([(jnc[1],jnc[2])])
        conversion_stats = collect_stats(samfile,sjintervals)
        junc_conv = conversion_stats + jnc + [row['Name']]
        rows.append(junc_conv)

    df = pd.DataFrame(rows, columns=["spliced","t_coverage","spliced_coverted","chr","start","end","name"])

    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-r", "--regions")
    parser.add_argument("-o", "--outputfolder")
    parser.add_argument("-c", "--cpu")


    args = parser.parse_args()

    bampath = args.bam
    bedpath = args.regions
    outfolder = args.outputfolder

    basenameBam = Path(bampath).stem
    basenameBed = Path(bedpath).stem
    
    junctions = pyranges.readers.read_bed(bedpath, as_df=True, nrows=None)
    chroms = junctions['Chromosome']
    # chroms = set(a_g_dict.keys()).union(set(t_c_dict.keys()))
    #only normal chromosomes
    chroms = ["chr" + str(x+1) for x in range(22)]
    chroms.append("chrX")
    chroms.append("chrY")
    chroms.append("chrM")

    pool = Pool(processes=cpu)

    for x in range(len(chroms)):
        pool.apply_async(snpmask_bams,(bampath,a_g_dict, t_c_dict,chroms[x]))
    #countReads(chroms[x],bamfh)
    pool.close()
    pool.join()
    outputfile = os.path.join(outfolder, basenameBam + "_" + basenameBed + "_" + "spliced_counts.csv")
    
    counts = process_bam(bampath,bedpath)


    counts.to_csv(outputfile,index = False)



    return 0

if __name__ == "__main__":
    main()

