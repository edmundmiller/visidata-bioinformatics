# Visidata Plugins for Bioinformatics Formats

Formats:

- GFF2/3
- BED

## Installation and Usage

1. Save the code as `gff.py` in your VisiData plugins directory (`~/.visidata/plugins/`)
2. Launch VisiData with any GFF file: `vd your_file.gff`

Now your GFF attributes are neatly spread across columns, making it easy to sort, filter, and analyze your genomic features.

### Features

- Automatic handling of BED format variations (3-12 columns)
- Color-coding of strands (+ in red, - in blue)
- Chromosome-level statistics (counts, lengths, etc.)
- Region merging and filtering capabilities
- Distance calculations between adjacent regions
- Detailed region information view
- Fallback to TSV for non-standard BED files

# Further Reading

- [Visidata Docs on Plugins](https://www.visidata.org/docs/plugins/)
- [BED Format Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)


