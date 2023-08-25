import os
from glob import glob
import pandas as pd
import argparse


def readABCounts(file):
    """read `{sample_id}.AB.match.counts.txt` file as data frame
    :type file: str
    """
    df = pd.read_csv(file, sep='\t', header=None)
    df.rename(columns={0: 'sgID_AB', 1: 'count'}, inplace=True)

    return df


def makeABCountMatrix(countFiles):
    """read and merge `{sample_id}.AB.match.counts.txt` counts files
    :type countFiles: list
    """

    samples = [os.path.basename(f).split('.')[0] for f in countFiles]

    sgIDs = []
    for file in countFiles:
        for sg_ab in readABCounts(file)['sgID_AB'].tolist():
            sgIDs.append(sg_ab)

    counts = pd.DataFrame(index=set(sgIDs), columns=samples)

    for sample, file in zip(samples, countFiles):
        counts.loc[
            readABCounts(file)['sgID_AB'].tolist(), sample
        ] = readABCounts(file)['count'].tolist()

    counts = counts.fillna(0)

    counts.insert(0, 'target', counts.index.str.split('_').str[0])

    return counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Load `{sample_id}.AB.match.counts.txt` count files and merge them into one counts matrix.'
                    'This is the input to calculate the enrichment of each guide between screen arms.'
    )
    parser.add_argument('Count_Files_Path', help='Directory where module1 output files located.')
    parser.add_argument('Out_File_Path', help='Directory where output files should be written.')

    args = parser.parse_args()

    countsDirectory = args.Count_Files_Path
    if countsDirectory[-1] == '/':
        countsDirectory = countsDirectory[:-1]

    outputDirectory = args.Out_File_Path
    if outputDirectory[-1] == '/':
        outputDirectory = outputDirectory[:-1]

    if not os.path.exists(outputDirectory):
        # Create a new directory because it does not exist
        os.makedirs(outputDirectory)
        print("Output Directory is created!")

    countFiles = glob(f'{args.Count_Files_Path}/*.AB.match.counts.txt')
    counts = makeABCountMatrix(countFiles)

    # write raw counts table
    counts.to_csv(f'{outputDirectory}/counts.txt', sep='\t')
