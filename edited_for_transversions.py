import os
import pysam
import argparse
import pickle


TDP43kd_1_8h.sorted.bam real	123m36.579s
def read_fasta(filename):
    """
    This is a simple function written in base python that
    returns a dictionary made from a fasta file
    """
    # Read in the transcript fasta
    fasta = {}


    with open(filename, 'r') as file:
        for i, line in enumerate(file):

            if line.rstrip()[0:1] == ">":
                this_tx_name = line.rstrip().replace(">", "").split(" ")[0]
                print("reading " + this_tx_name)
            else:
                try:
                    fasta[this_tx_name].append(line.rstrip())
                except KeyError:
                    fasta[this_tx_name] = [line.rstrip()]

    for key in fasta.keys():
        fasta[key] = ''.join(fasta[key]).upper()
    return fasta


def gen_dict_matrix():
    """
    Generate a dictionary that returns the index for a given transversion
    should be used like: d[reference_base+read_base]
    """

    nts = ["A", "C", "G", "T", "N"]

    d = {}

    i = 0
    for ref in nts:
        for read in nts:
            d[ref+read] = i
            i += 1

    return d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam")
    parser.add_argument("-f", "--fasta")
    parser.add_argument("-p", "--picklefasta")
    args = parser.parse_args()

    matrix_dict = gen_dict_matrix()

    # Reading fasta
    
    if(args.fasta is  None):
        print("Reading pickled")
        with open(args.picklefasta, 'rb') as f:
            fasta = pickle.load(f)
    if(args.fasta is not None):
        print("Reading fasta")
        fasta = read_fasta(args.fasta)
    

    skipped = []

    infile = pysam.AlignmentFile(args.bam, "rb")

    outfile_path = os.path.splitext(args.bam)[0] + ".tagged.bam"
    outfile = pysam.AlignmentFile(outfile_path, "wb", template=infile)


    for i, record in enumerate(infile):

        rname = record.reference_name

        if rname is None:
            continue

        if rname not in fasta.keys():
            if rname not in skipped:
                skipped.append(rname)
                print("Skipping read aligned to unknown sequence " + rname)
            continue

        matrix = [0]*25

        positions = record.positions
        seq = record.query_alignment_sequence

        for j, p in enumerate(positions):
            matrix[matrix_dict[fasta[rname][p]+seq[j]]] += 1

        output = ','.join([str(a) for a in matrix])
        record.tags = record.tags + [('RA', output)]
        outfile.write(record)
    
    infile.close()
    outfile.close()



if __name__ == "__main__":
    main()
