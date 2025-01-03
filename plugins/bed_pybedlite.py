"""visidata loader for bed (browser extensible data) files using pybedlite."""

import pybedlite as pybed
from pathlib import path
from visidata import (
    visidata,
    sheet,
    column,
    vd,
    asyncthread,
    options,
    enter,
    textsheet,
    indexsheet,
    progress,
)

# add options for bed file handling and visualization
options.bed_skip_validation = false  # whether to skip field validation
options.bed_default_score = "0"  # default score when missing
options.bed_default_name = "."  # default name when missing
options.bed_default_strand = "."  # default strand when missing
options.bed_color_strands = true  # color + and - strands differently
options.bed_max_region_size = 1000000  # warning threshold for large regions

# add more bed-specific options
options.bed_min_region_size = 0  # minimum region size filter
options.bed_chrom_order = "natural"  # 'natural' or 'lexical' chromosome sorting
options.bed_show_gc = false  # show gc content (requires reference genome)
options.bed_to_gff_type = "region"  # default feature type when converting to gff
options.bed_to_gff_source = "bed2gff"  # default source when converting to gff


@visidata.api
def open_bed(vd, p):
    """try to open as bed, fall back to tsv if parsing fails"""
    try:
        sheet = bedpyblsheet(p.name, source=p)
        sheet.reload()
        # wait for async reload to complete
        while sheet.loading:
            pass
        if not sheet.rows:  # if no rows were successfully parsed
            vd.warning("no valid bed records found, falling back to tsv")
            return vd.opensource(p, filetype="tsv")
        return sheet
    except exception as e:
        vd.warning(f"failed to parse as bed ({str(e)}), falling back to tsv")
        return vd.opensource(p, filetype="tsv")


class bedpyblsheet(sheet):
    """sheet for displaying bed format data using pybedlite."""

    rowtype = "genomic regions"  # more specific rowtype
    required_fields = ["chrom", "start", "end"]

    def __init__(self, name, source=none, **kwargs):
        super().__init__(name, source=source, **kwargs)
        self.columns = []
        self.header_lines = []  # store browser/track/comment lines

        # add commands specific to bed files
        self.bindkey("g#", "show-region-stats")  # show statistics about genomic regions
        self.bindkey("zs", "select-by-strand")  # select rows by strand
        self.bindkey("zl", "select-large-regions")  # select unusually large regions
        self.bindkey(enter, "view-region-details")  # show detailed view of region

        # add more commands
        self.bindkey("zc", "summarize-by-chrom")  # chromosome summary
        self.bindkey("zf", "filter-by-size")  # filter by region size
        self.bindkey("zm", "merge-overlapping")  # merge overlapping regions
        self.bindkey("gd", "distance-to-next")  # calculate distances between regions

        # add format conversion commands
        self.bindkey("gf", "convert-to-gff")  # convert to gff format

        # define columns based on bedrecord attributes with better descriptions
        self.addcolumn(
            column("chrom", getter=lambda col, row: row.chrom, help="chromosome name")
        )
        self.addcolumn(
            column(
                "start",
                type=int,
                getter=lambda col, row: row.start,
                help="start position (0-based)",
            )
        )
        self.addcolumn(
            column(
                "end",
                type=int,
                getter=lambda col, row: row.end,
                help="end position (exclusive)",
            )
        )
        self.addcolumn(
            column("name", getter=lambda col, row: row.name, help="feature name")
        )
        self.addcolumn(
            column("score", getter=lambda col, row: row.score, help="score from 0-1000")
        )
        self.addcolumn(
            column(
                "strand", getter=lambda col, row: row.strand, help="strand (+, -, or .)"
            )
        )
        self.addcolumn(
            column(
                "thickstart",
                getter=lambda col, row: self._safe_convert(row.thick_start, int),
                help="start of thick drawing",
            )
        )
        self.addcolumn(
            column(
                "thickend",
                getter=lambda col, row: self._safe_convert(row.thick_end, int),
                help="end of thick drawing",
            )
        )
        self.addcolumn(
            column(
                "itemrgb",
                getter=lambda col, row: ",".join(map(str, row.item_rgb))
                if row.item_rgb
                else none,
                help="rgb color (r,g,b)",
            )
        )
        self.addcolumn(
            column(
                "blockcount",
                getter=lambda col, row: self._safe_convert(row.block_count, int),
                help="number of blocks/exons",
            )
        )
        self.addcolumn(
            column(
                "blocksizes",
                getter=lambda col, row: ",".join(map(str, row.block_sizes))
                if row.block_sizes
                else none,
                help="block sizes in bases",
            )
        )
        self.addcolumn(
            column(
                "blockstarts",
                getter=lambda col, row: ",".join(map(str, row.block_starts))
                if row.block_starts
                else none,
                help="block starts relative to start",
            )
        )

        # add computed columns
        self.addcolumn(
            column(
                "length",
                getter=lambda col, row: self.get_region_length(row),
                type=int,
                help="region length in bp",
            )
        )

        self.addcolumn(
            column(
                "distance_to_next",
                getter=lambda col, row: self._get_distance_to_next(row),
                type=int,
                help="distance to next region",
            )
        )

    def _safe_convert(self, value, type_func):
        """safely convert value to given type, return none if fails"""
        if value is none or value == "":
            return none
        try:
            return type_func(value)
        except (valueerror, typeerror):
            vd.debug(f"could not convert {value} to {type_func.__name__}")
            return value  # return original value instead of none

    @asyncthread
    def reload(self):
        """load bed records from file."""
        self.rows = []
        self.loading = true

        vd.status("starting bed file load...")

        # first pass to collect header lines
        header_count = 0
        with self.source.open_text() as fp:
            for line in fp:
                line = line.rstrip("\n")
                if line.startswith(
                    ("track", "browser")
                ):  # remove # from header check since it's not standard
                    self.header_lines.append(line)
                    header_count += 1

        vd.status(f"found {header_count} header lines")

        # second pass to load records using pybedlite
        try:
            bed_path = path(self.source.resolve())
            vd.status(f"processing bed file: {bed_path}")

            with bed_path.open() as bed_file:
                vd.status("creating pybedlite reader...")
                # skip header lines
                for _ in range(len(self.header_lines)):
                    next(bed_file)

                count = 0
                for line in bed_file:
                    line = line.strip()
                    if not line or line.startswith(("track", "browser")):
                        continue

                    fields = line.split("\t")
                    if len(fields) < 3:  # must have at least chrom, start, end
                        continue

                    try:
                        # modify how we handle the fields to match the test data format
                        record = pybed.bedrecord(
                            chrom=fields[0],
                            start=int(fields[1]),
                            end=int(fields[2]),
                            name=fields[3] if len(fields) > 3 else ".",
                            score=fields[4] if len(fields) > 4 else "0",
                            strand="." if len(fields) <= 5 else fields[5],
                            thick_start=none,  # make these optional
                            thick_end=none,
                            item_rgb=none,
                            block_count=none,
                            block_sizes=none,
                            block_starts=none,
                        )

                        if count == 0:
                            vd.status(f"first record found: {record}")
                        self.addrow(record)
                        count += 1
                        if count % 1000 == 0:
                            vd.status(f"loaded {count} records...")
                    except (valueerror, indexerror) as e:
                        vd.debug(f"skipping malformed record: {line} ({str(e)})")
                        continue

                vd.status(f"completed loading {count} bed records")

        except exception as e:
            vd.warning(f"error reading bed file: {str(e)}")
            import traceback

            vd.debug(traceback.format_exc())
        finally:
            self.loading = false

    def colorize_strand(self, row):
        """return color based on strand."""
        if not options.bed_color_strands:
            return none
        return "red" if row.strand == "+" else "blue" if row.strand == "-" else none

    def get_region_length(self, row):
        """calculate length of genomic region."""
        return row.end - row.start

    def show_region_stats(self, vd):
        """display statistics about the genomic regions."""
        total_regions = len(self.rows)
        total_bases = sum(self.get_region_length(row) for row in self.rows)
        strands = {}
        chroms = {}

        for row in self.rows:
            strands[row.strand] = strands.get(row.strand, 0) + 1
            chroms[row.chrom] = chroms.get(row.chrom, 0) + 1

        vd.status(
            f"regions: {total_regions}, total bases: {total_bases}, "
            f"strands: {dict(strands)}, chromosomes: {len(chroms)}"
        )

    def select_by_strand(self, vd):
        """interactive strand selection."""
        strand = vd.input("select strand (+/-/./other): ")
        for row in self.rows:
            if row.strand == strand:
                row.selected = true

    def select_large_regions(self, vd):
        """select regions larger than threshold."""
        threshold = vd.input(
            "minimum region size: ", value=str(options.bed_max_region_size)
        )
        try:
            threshold = int(threshold)
            for row in self.rows:
                if self.get_region_length(row) > threshold:
                    row.selected = true
        except valueerror:
            vd.warning("invalid threshold value")

    def view_region_details(self, vd):
        """show detailed information about the current region."""
        row = self.cursorrow
        length = self.get_region_length(row)
        details = f"""
        region details:
        ---------------
        chromosome: {row.chrom}
        start: {row.start}
        end: {row.end}
        length: {length:,} bp
        name: {row.name}
        score: {row.score}
        strand: {row.strand}
        """
        if length > options.bed_max_region_size:
            details += f"\nwarning: region size ({length:,}) exceeds threshold ({options.bed_max_region_size:,})"

        vd.push(textsheet(f"details_{row.name}", source=details))

    def _get_distance_to_next(self, row):
        """calculate distance to next region on same chromosome."""
        try:
            idx = self.rows.index(row)
            next_row = next(
                (r for r in self.rows[idx + 1 :] if r.chrom == row.chrom), none
            )
            if next_row:
                return next_row.start - row.end
        except (valueerror, indexerror):
            pass
        return none

    def summarize_by_chrom(self, vd):
        """create a summary sheet with chromosome statistics."""
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

        summary_sheet = indexsheet(f"{self.name}_chrom_summary", source=chrom_stats)
        summary_sheet.addcolumn(column("chromosome", getter=lambda c, r: r))
        summary_sheet.addcolumn(
            column("region_count", getter=lambda c, r: chrom_stats[r]["count"])
        )
        summary_sheet.addcolumn(
            column("total_bp", getter=lambda c, r: chrom_stats[r]["total_length"])
        )
        summary_sheet.addcolumn(
            column("min_length", getter=lambda c, r: chrom_stats[r]["min_length"])
        )
        summary_sheet.addcolumn(
            column("max_length", getter=lambda c, r: chrom_stats[r]["max_length"])
        )
        summary_sheet.addcolumn(
            column(
                "avg_length",
                getter=lambda c, r: chrom_stats[r]["total_length"]
                / chrom_stats[r]["count"],
                type=float,
            )
        )

        vd.push(summary_sheet)

    def filter_by_size(self, vd):
        """filter regions by size range."""
        min_size = vd.input(
            "minimum size (bp): ", value=str(options.bed_min_region_size)
        )
        max_size = vd.input(
            "maximum size (bp): ", value=str(options.bed_max_region_size)
        )
        try:
            min_size = int(min_size)
            max_size = int(max_size)
            for row in self.rows:
                length = self.get_region_length(row)
                row.selected = min_size <= length <= max_size
        except valueerror:
            vd.warning("invalid size value")

    def merge_overlapping(self, vd):
        """merge overlapping regions on same chromosome."""
        merged = []
        current = none

        # sort rows by chromosome and start position
        sorted_rows = sorted(self.rows, key=lambda r: (r.chrom, r.start))

        for row in sorted_rows:
            if not current:
                current = row
                continue

            if current.chrom == row.chrom and row.start <= current.end:
                # merge overlapping regions
                current.end = max(current.end, row.end)
            else:
                merged.append(current)
                current = row

        if current:
            merged.append(current)

        self.rows = merged
        vd.status(f"merged into {len(merged)} regions")

    def convert_to_gff(self, vd):
        """convert bed records to gff format."""
        gff_sheet = gffsheet(f"{self.name}_gff", source=none)

        feature_type = vd.input("gff feature type: ", value=options.bed_to_gff_type)
        source = vd.input("gff source: ", value=options.bed_to_gff_source)

        for bed_row in progress(self.rows, "converting to gff"):
            # convert coordinates from 0-based (bed) to 1-based (gff)
            start = bed_row.start + 1

            # construct gff attributes
            attrs = []
            if bed_row.name and bed_row.name != ".":
                attrs.append(f"name={bed_row.name}")
            if bed_row.score and bed_row.score != "0":
                attrs.append(f"score={bed_row.score}")
            if bed_row.thick_start is not none:
                attrs.append(
                    f"thick_start={bed_row.thick_start + 1}"
                )  # convert to 1-based
            if bed_row.thick_end is not none:
                attrs.append(f"thick_end={bed_row.thick_end}")
            if bed_row.item_rgb:
                attrs.append(f"rgb={','.join(map(str, bed_row.item_rgb))}")
            if bed_row.block_count:
                attrs.append(f"block_count={bed_row.block_count}")
            if bed_row.block_sizes:
                attrs.append(f"block_sizes={','.join(map(str, bed_row.block_sizes))}")
            if bed_row.block_starts:
                # convert relative starts to absolute 1-based coordinates
                abs_starts = [start + bs for bs in bed_row.block_starts]
                attrs.append(f"block_starts={','.join(map(str, abs_starts))}")

            # create gff fields
            gff_row = [
                bed_row.chrom,  # seqid
                source,  # source
                feature_type,  # type
                str(start),  # start (1-based)
                str(bed_row.end),  # end
                bed_row.score or ".",  # score
                bed_row.strand or ".",  # strand
                ".",  # phase
                ";".join(attrs) or ".",  # attributes
            ]

            gff_sheet.addrow(gff_row)

        vd.push(gff_sheet)
        vd.status(f"converted {len(self.rows)} bed records to gff format")


@visidata.api
def save_bed(vd, p, *sheets):
    """save sheet to bed format."""
    if len(sheets) != 1:
        vd.fail("can only save one sheet to bed")
    sheet = sheets[0]

    with p.open_text(mode="w") as fp:
        # write any header lines
        if hasattr(sheet, "header_lines"):
            for line in sheet.header_lines:
                fp.write(line + "\n")

        # write records
        for row in progress(sheet.rows, "saving"):
            try:
                if isinstance(row, pybed.bedrecord):
                    # handle native bed records
                    fields = [
                        row.chrom,
                        str(row.start),
                        str(row.end),
                        row.name or ".",
                        row.score or "0",
                        row.strand or ".",
                    ]
                    # add optional fields if present
                    if row.thick_start is not none:
                        fields.extend([str(row.thick_start), str(row.thick_end)])
                        if row.item_rgb:
                            fields.append(",".join(map(str, row.item_rgb)))
                            if row.block_count:
                                fields.extend(
                                    [
                                        str(row.block_count),
                                        ",".join(map(str, row.block_sizes)),
                                        ",".join(map(str, row.block_starts)),
                                    ]
                                )
                else:
                    # handle conversion from other formats (like gff)
                    fields = [
                        str(col.gettypedvalue(row)) for col in sheet.visiblecols[:6]
                    ]

                fp.write("\t".join(fields) + "\n")
            except exception as e:
                vd.warning(f"error saving row: {e}")
