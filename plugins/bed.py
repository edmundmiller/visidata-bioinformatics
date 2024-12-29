# /// script
# requires-python = ">=3.8"
# dependencies = [
#   "pybedlite",
# ]
# ///

"""VisiData loader for BED (Browser Extensible Data) files."""

from visidata import *

@VisiData.api
def guess_bed(vd, p):
    """Guess if file is a BED format based on content"""
    with p.open_text() as fp:
        for line in fp:
            if line.startswith(("#", "track", "browser")):
                continue
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                try:
                    int(fields[1])
                    int(fields[2])
                    return dict(filetype="bed", _likelihood=9)
                except ValueError:
                    pass
            break
    return None

@VisiData.api
def open_bed(vd, p):
    return BedSheet(p.name, source=p)

class BedSheet(Sheet):
    """Sheet for displaying BED format data"""
    rowtype = "regions"  # rowdef: list of fields
    
    def iterload(self):
        with self.source.open_text() as fp:
            for line in Progress(fp, 'loading'):
                line = line.strip()
                if not line or line.startswith(('#', 'track', 'browser')):
                    continue
                yield line.split('\t')

    def reload(self):
        self.columns = []
        
        # Required BED fields
        self.addColumn(Column('chrom', 0))
        self.addColumn(Column('start', 1, type=int))
        self.addColumn(Column('end', 2, type=int))
        
        # Optional BED fields
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
        
        for name, i, typ in optional_cols:
            self.addColumn(Column(name, i, type=typ))

        super().reload()
