"""VisiData loader for BED (Browser Extensible Data) files.

Uses pybedlite to properly parse BED format files with the following fields:
- chrom      - chromosome name
- chromStart - starting position (0-based)
- chromEnd   - ending position (exclusive)
- name       - feature name (optional)
- score      - score between 0-1000 (optional) 
- strand     - + or - (optional)
- thickStart - starting position at which feature is drawn thickly (optional)
- thickEnd   - ending position at which feature is drawn thickly (optional)
- itemRgb    - RGB color value (optional)
- blockCount - number of blocks/exons (optional)
- blockSizes - comma-separated list of block sizes (optional)
- blockStarts - comma-separated list of block starts (optional)
"""

from visidata import *
import pybedlite as pybed
from pathlib import Path
import re

@VisiData.api
def guess_bed(vd, p):
    """Guess if file is a BED format based on content"""
    with p.open_text() as fp:
        # Skip any header lines
        for line in fp:
            if line.startswith(('#', 'track', 'browser')):
                continue
            if not line.strip():
                continue
            # Check if line matches BED format (tab-separated, at least 3 fields,
            # 2nd and 3rd fields are integers)
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                try:
                    int(fields[1])
                    int(fields[2])
                    return dict(filetype='bed', _likelihood=9)
                except ValueError:
                    pass
            break
    return None

@VisiData.api
def open_bed(vd, p):
    return BedSheet(p.name, source=p)

class BedSheet(Sheet):
    rowtype = 'regions'  # rowdef: BedRecord
    
    def iterload(self):
        with pybed.reader(path=Path(str(self.source))) as bed_reader:
            for record in Progress(bed_reader, 'loading'):
                yield record

    def reload(self):
        self.columns = []
        
        # Required fields
        self.addColumn(ColumnAttr('chrom', 'chrom'))
        self.addColumn(ColumnAttr('chromStart', 'start', type=int))
        self.addColumn(ColumnAttr('chromEnd', 'end', type=int))
        
        # Optional fields that are attributes of BedRecord
        optional_cols = [
            ('name', 'name', str),
            ('score', 'score', int),
            ('strand', 'strand', str),
            ('thickStart', 'thick_start', int),
            ('thickEnd', 'thick_end', int),
            ('itemRgb', 'item_rgb', str),
            ('blockCount', 'block_count', int),
            ('blockSizes', 'block_sizes', str),
            ('blockStarts', 'block_starts', str)
        ]
        
        # Add all possible columns - they'll be None if not present
        for name, attr, typ in optional_cols:
            self.addColumn(ColumnAttr(name, attr, type=typ))

        super().reload()
