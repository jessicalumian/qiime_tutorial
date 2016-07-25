First we will start by downloading the data. Copy the following command into your terminal.

```text
wget ftp://ftp.microbio.me/qiime/tutorial_files/moving_pictures_tutorial-1.9.0.tgz || curl -O ftp://ftp.microbio.me/qiime/tutorial_files/moving_pictures_tutorial-1.9.0.tgz
tar -xzf moving_pictures_tutorial-1.9.0.tgz
```

Now we will change directory into the correct folder.

```text
cd moving_pictures_tutorial-1.9.0/
cd illumina/
```
# Checking the mapping files for errors

The QIIME mapping file contains all of the per-sample metadata, including technical information such as primers and barcodes that were used for each sample, and information about the samples, including what body site they were taken from. In this data set we're looking at human microbiome samples from four sites on the bodies of two individuals at mutliple time points. The metadata in this case therefore includes a subject identifier, a timepoint, and a body site for each sample. You can review the map.tsv file at the link in the previous cell to see an example of the data (or view the published Google Spreadsheet version, which is more nicely formatted).

In this step, we run validate_mapping_file.py to ensure that our mapping file is compatible with QIIME.

```text
validate_mapping_file.py -o vmf-map/ -m map.tsv
```

In this case there were no errors, but if there were we would review the resulting HTML summary to find out what errors are present. You could then fix those in a spreadsheet program or text editor and rerun validate_mapping_file.py on the updated mapping file.

For the sake of illustrating what errors in a mapping file might look like, we've created a bad mapping file (map-bad.tsv). We'll next call validate_mapping_file.py on the file map-bad.tsv. Review the resulting HTML report. What are the issues with this mapping file?

```text
validate_mapping_file.py -o vmf-map-bad/ -m map-bad.tsv
```

# Demultiplexing and quality filtering sequences

We next need to demultiplex and quality filter our sequences (i.e. assigning barcoded reads to the samples they are derived from). In general, you should get separate fastq files for your sequence and barcode reads. Note that we pass these files while still gzipped. split_libraries_fastq.py can handle gzipped or unzipped fastq files. The default strategy in QIIME for quality filtering of Illumina data is described in Bokulich et al (2013).

```text
split_libraries_fastq.py -o slout/ -i forward_reads.fastq.gz -b barcodes.fastq.gz -m map.tsv
```

We can see how many sequences we ended up with using count_seqs.py.

```text
count_seqs.py -i slout/seqs.fna
```
