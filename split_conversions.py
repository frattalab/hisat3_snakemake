import os
import pysam
import argparse
import pickle




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-m", "--mimimum",  default=2)
    skipped = []

    infile = pysam.AlignmentFile(args.bam, "rb")

    outfileConverted_path = os.path.splitext(args.bam)[0] + ".convertedreads.bam"
    outfileConverted = pysam.AlignmentFile(outfileConverted_path, "wb", template=infile)

    outfileUnConverted_path = os.path.splitext(args.bam)[0] + ".UNconvertedreads.bam"
    outfileUnConverted = pysam.AlignmentFile(outfileConverted_path, "wb", template=infile)


    for i, read in enumerate(infile):

        if read.has_tag('Yf'):
            conversion_tag = read.get_tag("Yf")
            
            if conversion_tag >= args.mimimum:
                outfileConverted.write(record)
            else:
                outfileUnConverted.write(record)


            if i % 1_000_000 == 0:
                print (i)
                print(seq)
                print(output)
                
    infile.close()
    outfileConverted.close()
    outfileUnConverted.close()



if __name__ == "__main__":
    main()
