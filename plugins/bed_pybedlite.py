"""VisiData loader for BED (Browser Extensible Data) files using pybedlite."""

import pybedlite as pybed
from pathlib import Path
from visidata import (
    VisiData,
    Sheet,
    Column,
    vd,
    asyncthread,
    options,
    ENTER,
    TextSheet,
    IndexSheet,
)

# Add options for BED file handling and visualization
options.bed_skip_validation = False  # Whether to skip field validation
options.bed_default_score = "0"  # Default score when missing
options.bed_default_name = "."  # Default name when missing
options.bed_default_strand = "."  # Default strand when missing
options.bed_color_strands = True  # Color + and - strands differently
options.bed_max_region_size = 1000000  # Warning threshold for large regions

# Add more BED-specific options
options.bed_min_region_size = 0  # Minimum region size filter
options.bed_chrom_order = "natural"  # 'natural' or 'lexical' chromosome sorting
options.bed_show_gc = False  # Show GC content (requires reference genome)


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

        # Add commands specific to BED files
        self.bindkey("g#", "show-region-stats")  # Show statistics about genomic regions
        self.bindkey("zs", "select-by-strand")  # Select rows by strand
        self.bindkey("zl", "select-large-regions")  # Select unusually large regions
        self.bindkey(ENTER, "view-region-details")  # Show detailed view of region

        # Add more commands
        self.bindkey("zc", "summarize-by-chrom")  # Chromosome summary
        self.bindkey("zf", "filter-by-size")  # Filter by region size
        self.bindkey("zm", "merge-overlapping")  # Merge overlapping regions
        self.bindkey("gd", "distance-to-next")  # Calculate distances between regions

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

        # Add computed columns
        self.addColumn(
            Column(
                "length",
                getter=lambda col, row: self.get_region_length(row),
                type=int,
                help="Region length in bp",
            )
        )

        self.addColumn(
            Column(
                "distance_to_next",
                getter=lambda col, row: self._get_distance_to_next(row),
                type=int,
                help="Distance to next region",
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

    def colorize_strand(self, row):
        """Return color based on strand."""
        if not options.bed_color_strands:
            return None
        return "red" if row.strand == "+" else "blue" if row.strand == "-" else None

    def get_region_length(self, row):
        """Calculate length of genomic region."""
        return row.end - row.start

    def show_region_stats(self, vd):
        """Display statistics about the genomic regions."""
        total_regions = len(self.rows)
        total_bases = sum(self.get_region_length(row) for row in self.rows)
        strands = {}
        chroms = {}

        for row in self.rows:
            strands[row.strand] = strands.get(row.strand, 0) + 1
            chroms[row.chrom] = chroms.get(row.chrom, 0) + 1

        vd.status(
            f"Regions: {total_regions}, Total bases: {total_bases}, "
            f"Strands: {dict(strands)}, Chromosomes: {len(chroms)}"
        )

    def select_by_strand(self, vd):
        """Interactive strand selection."""
        strand = vd.input("Select strand (+/-/./other): ")
        for row in self.rows:
            if row.strand == strand:
                row.selected = True

    def select_large_regions(self, vd):
        """Select regions larger than threshold."""
        threshold = vd.input(
            "Minimum region size: ", value=str(options.bed_max_region_size)
        )
        try:
            threshold = int(threshold)
            for row in self.rows:
                if self.get_region_length(row) > threshold:
                    row.selected = True
        except ValueError:
            vd.warning("Invalid threshold value")

    def view_region_details(self, vd):
        """Show detailed information about the current region."""
        row = self.cursorRow
        length = self.get_region_length(row)
        details = f"""
        Region Details:
        ---------------
        Chromosome: {row.chrom}
        Start: {row.start}
        End: {row.end}
        Length: {length:,} bp
        Name: {row.name}
        Score: {row.score}
        Strand: {row.strand}
        """
        if length > options.bed_max_region_size:
            details += f"\nWARNING: Region size ({length:,}) exceeds threshold ({options.bed_max_region_size:,})"

        vd.push(TextSheet(f"details_{row.name}", source=details))

    def _get_distance_to_next(self, row):
        """Calculate distance to next region on same chromosome."""
        try:
            idx = self.rows.index(row)
            next_row = next(
                (r for r in self.rows[idx + 1 :] if r.chrom == row.chrom), None
            )
            if next_row:
                return next_row.start - row.end
        except (ValueError, IndexError):
            pass
        return None

    def summarize_by_chrom(self, vd):
        """Create a summary sheet with chromosome statistics."""
        chrom_stats = {}
        for row in self.rows:
            if row.chrom not in chrom_stats:
                chrom_stats[row.chrom] = {
                    "count": 0,
                    "total_length": 0,
                    "min_length": float("inf"),
                    "max_length": 0,
                }
            stats = chrom_stats[row.chrom]
            length = self.get_region_length(row)
            stats["count"] += 1
            stats["total_length"] += length
            stats["min_length"] = min(stats["min_length"], length)
            stats["max_length"] = max(stats["max_length"], length)

        summary_sheet = IndexSheet(f"{self.name}_chrom_summary", source=chrom_stats)
        summary_sheet.addColumn(Column("chromosome", getter=lambda c, r: r))
        summary_sheet.addColumn(
            Column("region_count", getter=lambda c, r: chrom_stats[r]["count"])
        )
        summary_sheet.addColumn(
            Column("total_bp", getter=lambda c, r: chrom_stats[r]["total_length"])
        )
        summary_sheet.addColumn(
            Column("min_length", getter=lambda c, r: chrom_stats[r]["min_length"])
        )
        summary_sheet.addColumn(
            Column("max_length", getter=lambda c, r: chrom_stats[r]["max_length"])
        )
        summary_sheet.addColumn(
            Column(
                "avg_length",
                getter=lambda c, r: chrom_stats[r]["total_length"]
                / chrom_stats[r]["count"],
                type=float,
            )
        )

        vd.push(summary_sheet)

    def filter_by_size(self, vd):
        """Filter regions by size range."""
        min_size = vd.input(
            "Minimum size (bp): ", value=str(options.bed_min_region_size)
        )
        max_size = vd.input(
            "Maximum size (bp): ", value=str(options.bed_max_region_size)
        )
        try:
            min_size = int(min_size)
            max_size = int(max_size)
            for row in self.rows:
                length = self.get_region_length(row)
                row.selected = min_size <= length <= max_size
        except ValueError:
            vd.warning("Invalid size value")

    def merge_overlapping(self, vd):
        """Merge overlapping regions on same chromosome."""
        merged = []
        current = None

        # Sort rows by chromosome and start position
        sorted_rows = sorted(self.rows, key=lambda r: (r.chrom, r.start))

        for row in sorted_rows:
            if not current:
                current = row
                continue

            if current.chrom == row.chrom and row.start <= current.end:
                # Merge overlapping regions
                current.end = max(current.end, row.end)
            else:
                merged.append(current)
                current = row

        if current:
            merged.append(current)

        self.rows = merged
        vd.status(f"Merged into {len(merged)} regions")
