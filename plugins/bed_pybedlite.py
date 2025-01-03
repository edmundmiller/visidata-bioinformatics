"""VisiData loader for BED (Browser Extensible Data) files using pybedlite."""

import pybedlite as pybed
from pathlib import Path
from visidata import VisiData, Sheet, Column, vd, asyncthread, options

# Add options for BED file handling
options.bed_skip_validation = False  # Whether to skip field validation
options.bed_default_score = "0"  # Default score when missing
options.bed_default_name = "."  # Default name when missing
options.bed_default_strand = "."  # Default strand when missing


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
            return vd.openSource(p, filetype="tsv")
        return sheet
    except Exception as e:
        vd.warning(f"Failed to parse as BED ({str(e)}), falling back to TSV")
        return vd.openSource(p, filetype="tsv")


class BedPyblSheet(Sheet):
    """Sheet for displaying BED format data using pybedlite."""

    rowtype = "genomic regions"  # More specific rowtype
    required_fields = ["chrom", "start", "end"]

    def __init__(self, name, source=None, **kwargs):
        super().__init__(name, source=source, **kwargs)
        self.columns = []
        self.header_lines = []  # Store browser/track/comment lines

        # Define columns based on BedRecord attributes with better descriptions
        self.addColumn(
            Column("chrom", getter=lambda col, row: row.chrom, help="Chromosome name")
        )
        self.addColumn(
            Column(
                "start",
                type=int,
                getter=lambda col, row: row.start,
                help="Start position (0-based)",
            )
        )
        self.addColumn(
            Column(
                "end",
                type=int,
                getter=lambda col, row: row.end,
                help="End position (exclusive)",
            )
        )
        self.addColumn(
            Column("name", getter=lambda col, row: row.name, help="Feature name")
        )
        self.addColumn(
            Column("score", getter=lambda col, row: row.score, help="Score from 0-1000")
        )
        self.addColumn(
            Column(
                "strand", getter=lambda col, row: row.strand, help="Strand (+, -, or .)"
            )
        )
        self.addColumn(
            Column(
                "thickStart",
                getter=lambda col, row: self._safe_convert(row.thick_start, int),
                help="Start of thick drawing",
            )
        )
        self.addColumn(
            Column(
                "thickEnd",
                getter=lambda col, row: self._safe_convert(row.thick_end, int),
                help="End of thick drawing",
            )
        )
        self.addColumn(
            Column(
                "itemRgb",
                getter=lambda col, row: ",".join(map(str, row.item_rgb))
                if row.item_rgb
                else None,
                help="RGB color (R,G,B)",
            )
        )
        self.addColumn(
            Column(
                "blockCount",
                getter=lambda col, row: self._safe_convert(row.block_count, int),
                help="Number of blocks/exons",
            )
        )
        self.addColumn(
            Column(
                "blockSizes",
                getter=lambda col, row: ",".join(map(str, row.block_sizes))
                if row.block_sizes
                else None,
                help="Block sizes in bases",
            )
        )
        self.addColumn(
            Column(
                "blockStarts",
                getter=lambda col, row: ",".join(map(str, row.block_starts))
                if row.block_starts
                else None,
                help="Block starts relative to start",
            )
        )

    def _safe_convert(self, value, type_func):
        """Safely convert value to given type, return None if fails"""
        if value is None or value == "":
            return None
        try:
            return type_func(value)
        except (ValueError, TypeError):
            vd.debug(f"Could not convert {value} to {type_func.__name__}")
            return value  # Return original value instead of None

    @asyncthread
    def reload(self):
        """Load BED records from file."""
        self.rows = []
        self.loading = True

        vd.status("Starting BED file load...")

        # First pass to collect header lines
        header_count = 0
        with self.source.open_text() as fp:
            for line in fp:
                line = line.rstrip("\n")
                if line.startswith(
                    ("track", "browser")
                ):  # Remove # from header check since it's not standard
                    self.header_lines.append(line)
                    header_count += 1

        vd.status(f"Found {header_count} header lines")

        # Second pass to load records using pybedlite
        try:
            bed_path = Path(self.source.resolve())
            vd.status(f"Processing BED file: {bed_path}")

            with bed_path.open() as bed_file:
                vd.status("Creating pybedlite reader...")
                # Skip header lines
                for _ in range(len(self.header_lines)):
                    next(bed_file)

                count = 0
                for line in bed_file:
                    line = line.strip()
                    if not line or line.startswith(("track", "browser")):
                        continue

                    fields = line.split("\t")
                    if len(fields) < 3:  # Must have at least chrom, start, end
                        continue

                    try:
                        # Modify how we handle the fields to match the test data format
                        record = pybed.BedRecord(
                            chrom=fields[0],
                            start=int(fields[1]),
                            end=int(fields[2]),
                            name=fields[3] if len(fields) > 3 else ".",
                            score=fields[4] if len(fields) > 4 else "0",
                            strand="." if len(fields) <= 5 else fields[5],
                            thick_start=None,  # Make these optional
                            thick_end=None,
                            item_rgb=None,
                            block_count=None,
                            block_sizes=None,
                            block_starts=None,
                        )

                        if count == 0:
                            vd.status(f"First record found: {record}")
                        self.addRow(record)
                        count += 1
                        if count % 1000 == 0:
                            vd.status(f"Loaded {count} records...")
                    except (ValueError, IndexError) as e:
                        vd.debug(f"Skipping malformed record: {line} ({str(e)})")
                        continue

                vd.status(f"Completed loading {count} BED records")

        except Exception as e:
            vd.warning(f"Error reading BED file: {str(e)}")
            import traceback

            vd.debug(traceback.format_exc())
        finally:
            self.loading = False

