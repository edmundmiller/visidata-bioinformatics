"""VisiData loader for BED (Browser Extensible Data) files using pybedlite."""

import pybedlite as pybed
from visidata import VisiData, Sheet, Column, vd, asyncthread

@VisiData.api
def open_bed_pybedlite(vd, p):
    """Open a BED file using pybedlite."""
    return BedPyblSheet(p.name, source=p)

class BedPyblSheet(Sheet):
    """Sheet for displaying BED format data using pybedlite."""
    rowtype = "regions"  # rowdef: BedRecord objects
    
    def __init__(self, name, source=None, **kwargs):
        super().__init__(name, source=source, **kwargs)
        self.columns = []
        
        # Define columns based on BedRecord attributes
        self.addColumn(Column("chrom", getter=lambda col,row: row.chrom))
        self.addColumn(Column("start", type=int, getter=lambda col,row: row.start))
        self.addColumn(Column("end", type=int, getter=lambda col,row: row.end))
        self.addColumn(Column("name", getter=lambda col,row: row.name))
        self.addColumn(Column("score", type=float, getter=lambda col,row: row.score))
        self.addColumn(Column("strand", getter=lambda col,row: "+" if not row.strand else "-" if row.strand.negative else "+"))
        self.addColumn(Column("thickStart", type=int, getter=lambda col,row: row.thick_start))
        self.addColumn(Column("thickEnd", type=int, getter=lambda col,row: row.thick_end))
        self.addColumn(Column("itemRgb", getter=lambda col,row: ','.join(map(str, row.item_rgb)) if row.item_rgb else None))
        self.addColumn(Column("blockCount", type=int, getter=lambda col,row: row.block_count))
        self.addColumn(Column("blockSizes", getter=lambda col,row: ','.join(map(str, row.block_sizes)) if row.block_sizes else None))
        self.addColumn(Column("blockStarts", getter=lambda col,row: ','.join(map(str, row.block_starts)) if row.block_starts else None))

    @asyncthread
    def reload(self):
        """Load BED records from file."""
        self.rows = []
        try:
            with pybed.reader(self.source.resolve()) as reader:
                for record in reader:
                    self.addRow(record)
        except Exception as e:
            vd.warning(f"Error loading BED file: {str(e)}")

        # Register BED format detection
        vd.option("filetype", "bed", "BED", BedPyblSheet)
