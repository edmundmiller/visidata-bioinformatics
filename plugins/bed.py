"""VisiData loader for BED (Browser Extensible Data) files."""

from copy import copy

from visidata import (
    Sheet,
    TsvSheet,
    Column,
    options,
    vd,
    VisiData,
    asyncthread,
    Progress,
)


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


class BedSheet(TsvSheet):
    """Sheet for displaying BED format data"""

    rowtype = "regions"  # rowdef: list of fields

    def __init__(self, name, source=None, **kwargs):
        super().__init__(name, source=source, delimiter="\t", headerlines=0, **kwargs)

    @asyncthread
    def reload(self):
        self.columns = []
        self.rows = []

        def make_getter(idx, type_func=str):
            def getter(col, row):
                try:
                    return type_func(row[idx]) if row and len(row) > idx else None
                except (IndexError, ValueError, TypeError):
                    return None

            return getter

        # Required BED fields
        self.addColumn(Column(name="chrom", type=str, getter=make_getter(0)))
        self.addColumn(Column(name="start", type=int, getter=make_getter(1, int)))
        self.addColumn(Column(name="end", type=int, getter=make_getter(2, int)))

        # Optional BED fields with their types
        optional_cols = [
            ("name", 3, str),  # Name of region
            ("score", 4, float),  # Score from 0-1000
            ("strand", 5, str),  # + or - for strand
            ("thickStart", 6, int),  # Starting position at which feature is drawn thickly
            ("thickEnd", 7, int),  # Ending position at which feature is drawn thickly
            ("itemRgb", 8, str),  # RGB color value (e.g., 255,0,0)
            ("blockCount", 9, int),  # Number of blocks (exons)
            ("blockSizes", 10, str),  # Comma-separated list of block sizes
            ("blockStarts", 11, str),  # Comma-separated list of block starts
        ]

        for name, idx, type_func in optional_cols:
            self.addColumn(Column(name=name, getter=make_getter(idx, type_func)))

        # Load the data
        with self.source.open_text() as fp:
            lines = fp.readlines()
            for line in Progress(lines, 'loading BED file'):
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(("browser", "#")):
                    continue
                if line.startswith("track"):
                    # Parse track line but don't add as row
                    continue
                try:
                    fields = line.split(self.delimiter)
                    # Ensure minimum 3 fields
                    if len(fields) < 3:
                        vd.warning(f"skipping line with too few fields: {line[:50]}...")
                        continue
                    # Validate integer fields
                    fields[1] = int(fields[1])  # start
                    fields[2] = int(fields[2])  # end
                    if len(fields) > 4:
                        try:
                            fields[4] = float(fields[4])  # score
                        except (ValueError, TypeError):
                            fields[4] = 0.0
                    self.addRow(fields)
                except Exception as e:
                    vd.warning(f"error parsing line: {line[:50]}... {str(e)}")
