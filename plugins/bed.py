"""VisiData loader for BED (Browser Extensible Data) files.

Handles the standard BED format with the following fields:
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

@VisiData.api
def open_bed(vd, p):
    return BedSheet(p.name, source=p)

class BedSheet(TsvSheet):
    rowtype = 'regions'  # rowdef: list of fields
    delimiter = '\t'
    
    def iterload(self):
        # Skip comment/header lines
        self.options.regex_skip = r'^(#|track|browser)'
        
        for row in super().iterload():
            # Ensure minimum required fields
            if len(row) < 3:
                continue
                
            # Convert positions to integers
            try:
                row[1] = int(row[1])  # chromStart
                row[2] = int(row[2])  # chromEnd
                
                if len(row) > 4:
                    row[4] = int(row[4])  # score
                if len(row) > 6:
                    row[6] = int(row[6])  # thickStart
                    row[7] = int(row[7])  # thickEnd
                if len(row) > 9:
                    row[9] = int(row[9])  # blockCount
            except ValueError as e:
                continue
                
            yield row

    def reload(self):
        self.columns = []
        
        # Required fields
        self.addColumn(ColumnItem('chrom', 0))
        self.addColumn(ColumnItem('chromStart', 1, type=int))
        self.addColumn(ColumnItem('chromEnd', 2, type=int))
        
        # Optional fields
        optional_cols = [
            ('name', 3, str),
            ('score', 4, int),
            ('strand', 5, str),
            ('thickStart', 6, int),
            ('thickEnd', 7, int),
            ('itemRgb', 8, str),
            ('blockCount', 9, int),
            ('blockSizes', 10, str),
            ('blockStarts', 11, str)
        ]
        
        # Add columns based on first data row
        for line in self.source.open_text():
            if self.options.regex_skip and re.match(self.options.regex_skip, line):
                continue
            if not line.strip():
                continue
                
            ncols = len(line.split(self.delimiter))
            for name, i, typ in optional_cols:
                if i < ncols:
                    self.addColumn(ColumnItem(name, i, type=typ))
            break

        super().reload()
