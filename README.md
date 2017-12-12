# 10x BAM to FASTQ converter

Tool for converting 10x BAMs produced by Cell Ranger or Long Ranger back to FASTQ files that can be used as inputs to re-run analysis.
The FASTQ files emitted by the tool should contain the same set of sequences that were input to the original pipeline run, although the order will
not be preserved.  The FASTQs will be emitted into a directory structure that is compatible with the directories created by the 'mkfastq' tool.

## Building
bamtofastq is standard Rust executable project, that works with stable Rust >=1.13.  Install Rust through the standard channels, then type `cargo build --release`.
The executable will appear at `target/release/bamtofastq`.  As usual it's important to use a release build to get good performance.

## Running

```
10x Genomics BAM to FASTQ converter.

Usage:
  bamtofastq [options] <bam> <output-path>
  bamtofastq (-h | --help)
  bamtofastq --version

Usage:
  bamtofastq [options] <bam> <output-path>
  bamtofastq (-h | --help)

Options:
  --locus=<locus>      Optional. Only include read pairs mapping to locus. Use chrom:start-end format.
  --reads-per-fastq=N  Number of reads per FASTQ chunk [default: 50000000]
  --gemcode            Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)
  --lr20               Convert a BAM produced by Longranger 2.0
  --cr11               Convert a BAM produced by Cell Ranger 1.0-1.1
  --bx-list=L          Only include BX values listed in text file L. Requires BX-sorted and index BAM file (see Long Ranger support for details).
  -h --help            Show this screen.
```  


## BAM file format support

10x BAMs produced by Long Ranger v2.1+ and Cell Ranger v1.2+ contain header fields that permit automatic conversion to the correct FASTQ sequences.
Older 10x pipelines require arguments to indicate which pipeline created the BAM.

### Sequence Layout

Special entries are inserted into @CO (comment) lines in the BAM header, indicating how to recover the original FASTQ sequences from the BAM record.

```
10x_bam_to_fastq:R1(RX:QX,TR:TQ,SEQ:QUAL)
10x_bam_to_fastq:R2(SEQ:QUAL)
10x_bam_to_fastq:I1(BC:QT)
```

The 'R1' line indicates that the R1 FASTQ sequence is composed of the RX tag, followed by the TR tag, followed by the sequence stored in the BAM record.
The QV sequence is composed of the QX tag, the TQ tag followed by the QVs stored in the record. The seqeunce and qv stored in the record need to be
reverse-complemented if the record has the 'reverse' flag set.

Some sequencing formats (in particular Single Cell 3') do not include generate alignment records for the R2 read, and the sequence of the R2 read 
are encoded in tags of the R1 alignment record. This mode is activated if the R2 specification does not contain a 'SEQ:QUAL' entry.


### Read groups

If the BAM file contains '@RG' headers and tags to indicate the source of each read, separate FASTQs will be created for each read group present
in the BAM file.


### Known Issues

* Multi-Gem Group BAM files created by Cell Ranger 1.2 and earlier do not carry Read Group tags, so reads from different GEM groups cannot be distinguished.
* 'Unaligned' BAM files created by the BASIC pipeline in Long Ranger versions prior to 2.1.3, due to missing R1/R2 flags.
