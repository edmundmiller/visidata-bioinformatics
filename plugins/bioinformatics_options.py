"""Shared options for bioinformatics plugins."""

from visidata import options

# BED options
options.bed_skip_validation = False
options.bed_default_score = "0"
options.bed_default_name = "."
options.bed_default_strand = "."
options.bed_color_strands = True
options.bed_max_region_size = 1000000
options.bed_min_region_size = 0
options.bed_chrom_order = "natural"
options.bed_show_gc = False
options.bed_to_gff_type = "region"
options.bed_to_gff_source = "bed2gff"

# GFF options
options.gff_to_bed_name_attr = "Name"
options.gff_to_bed_score_attr = "score"
options.gff_default_score = "0"
