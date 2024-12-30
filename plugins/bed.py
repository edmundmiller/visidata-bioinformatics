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
        self.addColumn(Column(name="chrom", getter=make_getter(0)))
        self.addColumn(Column(name="start", getter=make_getter(1, int)))
        self.addColumn(Column(name="end", getter=make_getter(2, int)))

        # Optional BED fields with their types
        optional_cols = [
            ("name", 3, str),
            ("score", 4, int),
            ("strand", 5, str),
            ("thickStart", 6, int),
            ("thickEnd", 7, int),
            ("itemRgb", 8, str),
            ("blockCount", 9, int),
            ("blockSizes", 10, str),
            ("blockStarts", 11, str),
        ]

        for name, idx, type_func in optional_cols:
            self.addColumn(Column(name=name, getter=make_getter(idx, type_func)))

        # Load the data
        with self.source.open_text() as fp:
            for line in Progress(fp.readlines(), 'loading BED file'):
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                try:
                    fields = line.split(self.delimiter)
                    self.addRow(fields)
                except Exception as e:
                    vd.warning(f"error parsing line: {line[:50]}... {str(e)}")
                    self.addRow([line]) # Add the problematic line as a single-column row
