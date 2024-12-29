"""VisiData loader for BED (Browser Extensible Data) files."""

from copy import copy

from visidata import Sheet, TsvSheet, Column, options, vd, VisiData


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

    def reload(self):
        super().reload()
        
        # Clear existing columns and add our BED-specific ones
        self.columns = []
        
        def make_getter(idx, type_func=str):
            def getter(col, row):
                try:
                    return type_func(row[idx])
                except (IndexError, ValueError):
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

    def iterload(self):
        with self.source.open_text() as fp:
            for line in fp:
                line = line.rstrip('\n')
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split(self.delimiter)
                yield fields


@VisiData.api
def open_usv(vd, p):
    return TsvSheet(p.base_stem, source=p, delimiter="\u241f", row_delimiter="\u241e")


@VisiData.api
def save_usv(vd, p, vs):
    vd.save_tsv(p, vs, row_delimiter="\u241e", delimiter="\u241f")
