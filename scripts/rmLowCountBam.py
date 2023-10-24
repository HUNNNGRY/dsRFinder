#!/usr/bin/env python2
import os, sys
import argparse
import pysam

def filter_bam(input_file, output_file, min_count):
    # Open the input BAM file
    bamfile = pysam.AlignmentFile(input_file, "rb")
    # if not os.path.exists(bamfile+'.bai'):
    #     sys.stderr.write('Indexing RNA-Seq BAM file.\n')
    #     pysam.index(bamfile)

    # Get the list of chromosomes in the BAM file
    chromosomes = bamfile.references

    # Create a new BAM file for the output
    outfile = pysam.AlignmentFile(output_file, "wb", template=bamfile)

    # Iterate over each chromosome in the BAM file
    for chrom in chromosomes:
        # Get the count of reads for the current chromosome
        count = bamfile.count(chrom)

        # If the count is greater than or equal to the minimum threshold, keep the chromosome
        if count >= min_count:
            # Iterate over each read in the current chromosome and write it to the output file
            for read in bamfile.fetch(chrom):
                outfile.write(read)

    # Close the input and output BAM files
    bamfile.close()
    outfile.close()

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Filter BAM file by chromosome count')
    parser.add_argument('--input', help='input BAM file')
    parser.add_argument('--output', help='output BAM file')
    parser.add_argument('--min-count', type=int, default=10, help='minimum count threshold')
    args = parser.parse_args()

    # Call the filter_bam function with the specified arguments
    filter_bam(args.input, args.output, args.min_count)