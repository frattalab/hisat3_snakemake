__author__ = '@aleighbrown'
import os
import pysam
import argparse
from pathlib import Path
import pyranges
import pandas as pd
from multiprocessing import Pool, cpu_count
from timeit import default_timer as timer


def collect_stats(infile,sj_info):
    spliced = 0
    t_coverage = 0
    spliced_coverted = 0

    sjintervals = set([(sj_info[1],sj_info[2])])

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

def process_bam(junctions,infile):
    rows = []
    

    samfile = pysam.AlignmentFile(infile, "rb")

    print('Processing Junctions:',flush=True)
    
    for index, row in junctions.iterrows():
        print(row['Name'],flush=True)
        jnc = [row['Chromosome'], row['Start'],row['End']]

        conversion_stats = collect_stats(samfile,jnc)
        
        junc_conv = conversion_stats + jnc + [row['Name']]
        rows.append(junc_conv)

    df = pd.DataFrame(rows, columns=["spliced","t_coverage","spliced_coverted","chr","start","end","name"])

    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-r", "--regions")
    parser.add_argument("-o", "--outputfolder")
    parser.add_argument("-t", "--threads")


    args = parser.parse_args()

    bampath = args.bam
    bedpath = args.regions
    outfolder = args.outputfolder
    threads = int(args.threads)



    basenameBam = Path(bampath).stem
    basenameBed = Path(bedpath).stem
    
    junctions = pyranges.readers.read_bed(bedpath, as_df=True, nrows=None)
    n = junctions.shape[0] // threads #chunk row size
    list_df = [junctions[i:i+n] for i in range(0,junctions.shape[0],n)]
    input_file = [bampath] * len(list_df) 
    junctions_bamtuple = tuple(zip(list_df,input_file))

    start = timer()

    print(f'starting computations on {threads} cores')

    

    with Pool() as pool:
        res = pool.starmap(process_bam, junctions_bamtuple)
        counts = pd.concat(res)


    end = timer()
    print(f'elapsed time: {end - start}')
    
    outputfile = os.path.join(outfolder, basenameBam + "_" + basenameBed + "_" + "spliced_counts.csv")
    

    counts.to_csv(outputfile,index = False)



    return 0

if __name__ == "__main__":
    main()

