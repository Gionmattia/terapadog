# About the testing data

This document serves the purpose to illustrate where the testing datasets were obtained and how they were modified to better suit the testing purposes.

## Data Source:

The files ribo_counts, rna_counts, and sample_info were downloaded from:

<https://github.com/SGDDNB/translational_regulation/tree/5812738bcd4219d1163eda26ece346fd491babf6/sample_data>

Which is github repository linked to the original DeltaTE paper, see [Chotani, S.et al. Current Protocols in Molecular Biology (2019)](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpmb.108)

## Modifications:

The datasets where then manually modified as it follows:

-   The extension was changed to .tsv
-   Sample names were changed between matching RNA-Seq and Ribo-Seq samples, so to match TERAPADOG's required formatting. \
    See the schema below for info:

| Original Name                                                                   | New name   | SeqType | Condition |
|---------------------|-----------------|-----------------|-----------------|
| Fib2424hr5R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam          | Sample2424 | RIBO    | 2         |
| Fib2724hr10R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam         | Sample2724 | RIBO    | 2         |
| Fib3424hr15R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam         | Sample3424 | RIBO    | 2         |
| Fib4124hr20R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam         | Sample4124 | RIBO    | 2         |
| Fib24bsl1R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam           | Sample24   | RIBO    | 1         |
| Fib27bsl6R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam           | Sample27   | RIBO    | 1         |
| Fib34bsl11R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam          | Sample34   | RIBO    | 1         |
| Fib41bsl16R1001trimmedNoAbundant.fastq.gzAligned.sortedByCoord.out.bam          | Sample41   | RIBO    | 1         |
| Pao24RIBO24hpA5R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam  | Sample2424 | RNA     | 2         |
| Pao27RIBO24hpA10R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam | Sample2724 | RNA     | 2         |
| Pao34RIBO24hpA15R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam | Sample3424 | RNA     | 2         |
| Pao41RIBO24hpA20R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam | Sample4124 | RNA     | 2         |
| Pao24RIBOBSLpA1R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam  | Sample24   | RNA     | 1         |
| Pao27RIBOBSLpA6R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam  | Sample27   | RNA     | 1         |
| Pao34RIBOBSLpA11R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam | Sample34   | RNA     | 1         |
| Pao41RIBOBSLpA16R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam | Sample41   | RNA     | 1         |

-   The sample_info file was also stripped of metadata which was not required by TERAPADOG. The column "SampleID" was replaced by "SampleName" ("New name" in the table above) and the "Batch" column was removed.

These operations were done only to illustrate users how an input file should be formatted and to help ease the learning curve for TERAPADOG.\
See the original paper mentioned for additional information on the biological source of the data.
