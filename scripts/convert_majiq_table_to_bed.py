import pyranges as pr
import pandas as pd
import argparse
from pathlib import Path

def convert_to_bed(path_to_csv,outfolder):
    new_bed_name = Path(path_to_csv).stem + 'cryptic_clusters.bed'
    new_bed_name = outfolder + new_bed_name
    df = pd.read_csv(path_to_csv)
    ori_colnames = list(df.columns)

    ori_colnames[11] = "baseline_PSI"
    ori_colnames[12] = "contrast_PSI"
    df.columns = ori_colnames
    mini_df = df[["seqnames", "start", "end", "gene_name", "paste_into_igv_junction", "junc_cat", "baseline_PSI","contrast_PSI", "strand"]]
    mini_df['Name'] = mini_df['gene_name'] + '|' + mini_df['junc_cat']
    mini_df['start'] = mini_df['start'] -1
    mini_df['end'] = mini_df['end'] -1

    mini_df = mini_df.rename(columns = {'seqnames': 'Chromosome',
                                        'start': 'Start',
                                        'end': 'End',
                                        'strand': 'Strand',
                                        'baseline_PSI': 'Baseline_PSI',
                                        'contrast_PSI': 'Contrast_PSI'
                                    })

    #Convert mini_df into a PyRanges object
    pyranges_df = pr.PyRanges(mini_df)

    #Cluster the ranges
    clustered_pyranges_df = pyranges_df.cluster()

    #Create a new dataframe from the clustered PyRanges object
    new_df = clustered_pyranges_df.df

    #Filter out rows where Baseline_PSI < 0.05 and Contrast_PSI > 0.1
    cluster_ids = new_df[(new_df['Baseline_PSI'] < 0.05) & (new_df['Contrast_PSI'] > 0.1)]['Cluster'].tolist()

    #Create a new dataframe from the clusters ids
    filtered_df = new_df[new_df['Cluster'].isin(cluster_ids)]
    #Get unique rows of filtered-df
    filtered_df = filtered_df[["Chromosome", "Start", "End", "Strand", "Name", "Cluster"]]
    filtered_df = filtered_df.drop_duplicates()

    pr.PyRanges(filtered_df).to_bed(new_bed_name)
    print(f'Ranges have been written to {new_bed_name}')
    
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--majiq")
    parser.add_argument("-o", "--outputfolder")


    args = parser.parse_args()

    csv_file = args.majiq
    outputfolder = args.outputfolder


    convert_to_bed(csv_file,outputfolder)




if __name__ == "__main__":
    main()
