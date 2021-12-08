import os
import pysam
import argparse




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-m", "--mimimum",  default=2)
    
    args = parser.parse_args()


    infile = pysam.AlignmentFile(args.bam, "rb")

    outfileConverted_path = os.path.splitext(args.bam)[0] + ".convertedreads.bam"
    outfileConverted = pysam.AlignmentFile(outfileConverted_path, "wb", template=infile)

    outfileUnConverted_path = os.path.splitext(args.bam)[0] + ".UNconvertedreads.bam"
    outfileUnConverted = pysam.AlignmentFile(outfileUnConverted_path, "wb", template=infile)


    for i, read in enumerate(infile):

        if read.has_tag('Yf'):
            conversion_tag = read.get_tag("Yf")
            
            if conversion_tag >= int(args.mimimum):
                outfileConverted.write(read)
            else:
                outfileUnConverted.write(read)


            if i % 1_000_000 == 0:
                print (i)

                
    infile.close()
    outfileConverted.close()
    outfileUnConverted.close()



if __name__ == "__main__":
    main()
