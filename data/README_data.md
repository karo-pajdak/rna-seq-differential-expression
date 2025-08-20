## Download FASTQ data from specific project (run in data/raw directory)
### Single-threaded
1. Download all .sra files first:
```
esearch -db sra -query PRJNA1224751 | efetch -format runinfo | cut -d',' -f1 | grep SRR | while read srr; do prefetch --output-directory . $srr; done
```
2. Then convert all .sra files to FASTQ
```
find . -name "*.sra" | while read sra_file; do fastq-dump --split-files --gzip "$sra_file"; done
```

Notes: 
- `esearch` finds the project in SRA 
- `efetch` retrieves run metadata (including SRR IDs) in csv format
- `cut` and `grep` take the first column (SRR IDs) and removes header
- `prefetch` downloads the .sra file and 'fastq-dump' converts to FASTQ

(See --help for each tool for more information)

### Parallel Version 
1. Download all .sra files first:
```
esearch -db sra -query PRJNA1224751 | efetch -format runinfo | cut -d',' -f1 | grep SRR | parallel -j 4 'prefetch --output-directory . {}'
```
2. Then convert all .sra files to FASTQ
```
find . -name "*.sra" | parallel -j 4 "fastq-dump --split-files --gzip {}"
```

Notes:
- `parallel -j N` runs N jobs simultaneously (if system allows it)

### Other option (less memory intensive)
To download FastQ files directly from Sequence Read Archive, use https://sra-explorer.info/#. This is allows you to skip the step of downloading .sra files first and then converting. You can search by studies, samples, experiments, or runs and then have the option to download with either:
- Raw FastQ Download URLs
- Bash script for downloading FastQ files