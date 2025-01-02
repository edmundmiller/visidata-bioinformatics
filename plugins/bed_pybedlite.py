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
        # Wait for async reload to complete
        while sheet.loading:
            pass
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
        self.addColumn(Column("score", getter=lambda col, row: row.score))  # Keep as string to handle non-numeric values
        self.addColumn(Column("strand", getter=lambda col, row: row.strand))
        self.addColumn(Column("thickStart", getter=lambda col, row: self._safe_convert(row.thick_start, int)))
        self.addColumn(Column("thickEnd", getter=lambda col, row: self._safe_convert(row.thick_end, int)))
        self.addColumn(Column("itemRgb", getter=lambda col, row: ",".join(map(str, row.item_rgb)) if row.item_rgb else None))
        self.addColumn(Column("blockCount", getter=lambda col, row: self._safe_convert(row.block_count, int)))
        self.addColumn(Column("blockSizes", getter=lambda col, row: ",".join(map(str, row.block_sizes)) if row.block_sizes else None))
        self.addColumn(Column("blockStarts", getter=lambda col, row: ",".join(map(str, row.block_starts)) if row.block_starts else None))

    def _safe_convert(self, value, type_func):
        """Safely convert value to given type, return None if fails"""
        if value is None or value == '':
            return None
        try:
            return type_func(value)
        except (ValueError, TypeError):
            vd.debug(f'Could not convert {value} to {type_func.__name__}')
            return value  # Return original value instead of None

    def reload(self):
        """Load BED records from file."""
        self.rows = []
        self.loading = True
        
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
                # Skip header lines
                for _ in range(len(self.header_lines)):
                    next(bed_file)
                    
                count = 0
                for line in bed_file:
                    line = line.strip()
                    if not line or line.startswith(('#', 'browser', 'track')):
                        continue
                        
                    fields = line.split('\t')
                    if len(fields) < 3:  # Must have at least chrom, start, end
                        continue
                        
                    try:
                        record = pybed.BedRecord(
                            chrom=fields[0],
                            start=int(fields[1]),
                            end=int(fields[2]),
                            name=fields[3] if len(fields) > 3 else None,
                            score=fields[4] if len(fields) > 4 else None,  # Keep as string
                            strand=fields[5] if len(fields) > 5 else None,
                            thick_start=int(fields[6]) if len(fields) > 6 else None,
                            thick_end=int(fields[7]) if len(fields) > 7 else None,
                            item_rgb=fields[8].split(',') if len(fields) > 8 and fields[8] != '.' else None,
                            block_count=int(fields[9]) if len(fields) > 9 and fields[9] != '.' else None,
                            block_sizes=[int(x) for x in fields[10].rstrip(',').split(',')] if len(fields) > 10 and fields[10] != '.' else None,
                            block_starts=[int(x) for x in fields[11].rstrip(',').split(',')] if len(fields) > 11 and fields[11] != '.' else None
                        )
                        if count == 0:
                            vd.status(f'First record found: {record}')
                        self.addRow(record)
                        count += 1
                        if count % 1000 == 0:
                            vd.status(f'Loaded {count} records...')
                    except (ValueError, IndexError) as e:
                        vd.debug(f'Skipping malformed record: {line} ({str(e)})')
                        continue
                    
                vd.status(f'Completed loading {count} BED records')
                
        except Exception as e:
            vd.warning(f"Error reading BED file: {str(e)}")
            import traceback
            vd.debug(traceback.format_exc())
        finally:
            self.loading = False

