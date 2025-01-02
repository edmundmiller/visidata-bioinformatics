"""VisiData loader for BED (Browser Extensible Data) files using pybedlite."""

import pybedlite as pybed
from pathlib import Path
from visidata import VisiData, Sheet, Column, vd, asyncthread, options


@VisiData.api
def open_bed(vd, p):
    """Try to open as BED, fall back to TSV if parsing fails"""
    try:
        sheet = BedPyblSheet(p.name, source=p)
        sheet.reload()
        if not sheet.rows:  # If no rows were successfully parsed
            vd.warning("No valid BED records found, falling back to TSV")
            return vd.openSource(p, filetype='tsv')
        return sheet
    except Exception as e:
        vd.warning(f"Failed to parse as BED ({str(e)}), falling back to TSV")
        return vd.openSource(p, filetype='tsv')


class BedPyblSheet(Sheet):
    """Sheet for displaying BED format data using pybedlite."""

    rowtype = "regions"  # rowdef: BedRecord objects

    def __init__(self, name, source=None, **kwargs):
        super().__init__(name, source=source, **kwargs)
        self.columns = []
        self.header_lines = []  # Store browser/track/comment lines

        # Define columns based on BedRecord attributes
        self.addColumn(Column("chrom", getter=lambda col, row: row.chrom))
        self.addColumn(Column("start", type=int, getter=lambda col, row: row.start))
        self.addColumn(Column("end", type=int, getter=lambda col, row: row.end))
        self.addColumn(Column("name", getter=lambda col, row: row.name))
        self.addColumn(Column("score", type=float, getter=lambda col, row: row.score))
        self.addColumn(Column("strand", getter=lambda col, row: row.strand))
        self.addColumn(Column("thickStart", type=int, getter=lambda col, row: row.thick_start))
        self.addColumn(Column("thickEnd", type=int, getter=lambda col, row: row.thick_end))
        self.addColumn(Column("itemRgb", getter=lambda col, row: ",".join(map(str, row.item_rgb)) if row.item_rgb else None))
        self.addColumn(Column("blockCount", type=int, getter=lambda col, row: row.block_count))
        self.addColumn(Column("blockSizes", getter=lambda col, row: ",".join(map(str, row.block_sizes)) if row.block_sizes else None))
        self.addColumn(Column("blockStarts", getter=lambda col, row: ",".join(map(str, row.block_starts)) if row.block_starts else None))

    @asyncthread
    def reload(self):
        """Load BED records from file."""
        self.rows = []
        
        vd.status('Starting BED file load...')
        
        # First pass to collect header lines
        header_count = 0
        with self.source.open_text() as fp:
            for line in fp:
                line = line.rstrip('\n')
                if line.startswith(('#', 'browser', 'track')):
                    self.header_lines.append(line)
                    header_count += 1
        
        vd.status(f'Found {header_count} header lines')
                    
        # Second pass to load records using pybedlite
        try:
            bed_path = Path(self.source.resolve())
            vd.status(f'Processing BED file: {bed_path}')
            
            # Debug: Print first few lines of file
            with bed_path.open() as preview:
                first_lines = [next(preview) for _ in range(5)]
                vd.status(f'First line preview: {first_lines[0].strip()}')
            
            with bed_path.open() as bed_file:
                vd.status('Creating pybedlite reader...')
                reader = pybed.reader(bed_file)
                
                count = 0
                for record in reader:
                    if count == 0:
                        vd.status(f'First record found: {record}')
                    self.addRow(record)
                    count += 1
                    if count % 1000 == 0:
                        vd.status(f'Loaded {count} records...')
                    
                vd.status(f'Completed loading {count} BED records')
                
        except Exception as e:
            vd.warning(f"Error reading BED file: {str(e)}")
            import traceback
            vd.debug(traceback.format_exc())

