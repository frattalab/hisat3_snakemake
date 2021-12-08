__author__ = '@oscarwilkins'
import os
import pysam
import argparse
import pickle




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
    parser.add_argument("-r", "--read1", default=1, type=int, help="Which read correponds to the forward strand of the RNA.\
         Default = 1, alternative = 2 (if read 2 is the RNA sense strand)")
    args = parser.parse_args()

    assert args.read1 in [1, 2], "Must be 1 or 2"

    rev_c_d = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

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
        mismatch_pos = []

        positions = record.get_reference_positions(full_length=True)   # fix 1
        seq = record.query_sequence # fix 3

        # Check what genomic strand the original RNA is from
        if record.is_reverse:  #
            if record.is_read1:
                strand = "-"
            else:
                strand = "+"
        else:
            if record.is_read1:
                strand = "+"
            else:
                strand = "-"

        # Flip if read 2 is the sense strand
        if args.read1 == 2:
            if strand == "+":
                strand = "-"
            else:
                strand = "+"            


        for j, p in enumerate(positions):
            if p is None:   # fix 2
                continue
            try:

                if strand == "-":
                    ref_RNA_base = rev_c_d[fasta[rname][p]]
                    read_RNA_base = rev_c_d[seq[j]]
                else:
                    ref_RNA_base = fasta[rname][p]
                    read_RNA_base = seq[j]                    

                matrix[matrix_dict[ref_RNA_base + read_RNA_base]] += 1

                if ref_RNA_base != read_RNA_base:
                    mmtype = matrix_dict[ref_RNA_base+read_RNA_base]
                    #MP:Z:2:38:38,16:70:70,10:80:80,7:110:110,10:189:189,16:202:202
                    ###TODO calculate the relative position from left most on the read
                    relative_position = p - record.query_alignment_start
                    mismatch_pos.append(f'{mmtype}:{j}:{relative_position}')
            except:
                print(f'failure: {rname}')
                print(positions)
                print(f'p is: {p}')
                print(record)
                print(seq)
                break
        
        ra_tag = ','.join([str(a) for a in matrix])
        record.tags = record.tags + [('RA', ra_tag)]
        
        if len(mismatch_pos) > 0:
            ma_tag = ','.join([str(a) for a in mismatch_pos])
            record.tags = record.tags + [('MP', ma_tag)]


        outfile.write(record)

        if i % 1_000_000 == 0:
            print (f'{i} reads read')

            
    infile.close()
    outfile.close()
    print("I finished tagging without any problems")



if __name__ == "__main__":
    main()
