# demultiplex_fastq
Simple, fast and memory efficient demultiplexer for FASTQ sequencing files
written by Andy Hauser <andreas.hauser@LMU.de>

# Usage

> Usage: demultiplex_fastq [-p PREFIX] --r1 FASTQ1 [--r2 FASTQ2] --i1 FASTQ_INDEX [--i2 FASTQ_INDEX2] -b BARCODE1,BARCODE2[,...]

# Example

```
$ demultiplex_fastq --r1 lane8_R1.fastq.gz --i1 lane8_R2.fastq.gz -b ACACGC,GGTATA
```

