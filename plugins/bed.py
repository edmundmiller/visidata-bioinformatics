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

class BedSheet(Sheet):
    rowtype = 'regions'  # rowdef: list of fields
    
    def iterload(self):
        for line in Progress(self.source.open_text(), 'loading'):
            # Skip comment/header lines
            if line.startswith(('#', 'track', 'browser')):
                continue
            if not line.strip():
                continue
                
            fields = line.rstrip('\n\r').split('\t')
            
            # Ensure minimum required fields
            if len(fields) < 3:
                continue
                
            # Convert positions to integers
            try:
                fields[1] = int(fields[1])  # chromStart
                fields[2] = int(fields[2])  # chromEnd
                
                if len(fields) > 4:
                    fields[4] = int(fields[4])  # score
                if len(fields) > 6:
                    fields[6] = int(fields[6])  # thickStart
                    fields[7] = int(fields[7])  # thickEnd
                if len(fields) > 9:
                    fields[9] = int(fields[9])  # blockCount
            except ValueError as e:
                continue
                
            yield fields

    def reload(self):
        self.columns = []
        
        # Required fields
        self.addColumn(ColumnItem('chrom', 0))
        self.addColumn(ColumnItem('chromStart', 1, type=int))
        self.addColumn(ColumnItem('chromEnd', 2, type=int))
        
        # Optional fields based on number of columns in first data row
        for line in self.source.open_text():
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
                
            fields = line.split('\t')
            if len(fields) > 3:
                self.addColumn(ColumnItem('name', 3))
            if len(fields) > 4:
                self.addColumn(ColumnItem('score', 4, type=int))
            if len(fields) > 5:
                self.addColumn(ColumnItem('strand', 5))
            if len(fields) > 7:
                self.addColumn(ColumnItem('thickStart', 6, type=int))
                self.addColumn(ColumnItem('thickEnd', 7, type=int))
            if len(fields) > 8:
                self.addColumn(ColumnItem('itemRgb', 8))
            if len(fields) > 11:
                self.addColumn(ColumnItem('blockCount', 9, type=int))
                self.addColumn(ColumnItem('blockSizes', 10))
                self.addColumn(ColumnItem('blockStarts', 11))
            break

        super().reload()
