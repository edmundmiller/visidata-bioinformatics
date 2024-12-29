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
        super().__init__(name, source=source, delimiter="\t", **kwargs)

    def reload(self):
        self.columns = []
        # Required BED fields
        self.addColumn(Column("chrom", 0))
        self.addColumn(Column("start", 1, type=int))
        self.addColumn(Column("end", 2, type=int))

        # Optional BED fields
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

        for name, i, typ in optional_cols:
            self.addColumn(Column(name, i, type=typ))

        super().reload()

    def iterload(self):
        for line in super().iterload():
            if not line or (
                isinstance(line[0], str)
                and line[0].startswith(("#", "track", "browser"))
            ):
                continue
            yield line


@VisiData.api
def open_usv(vd, p):
    return TsvSheet(p.base_stem, source=p, delimiter="\u241f", row_delimiter="\u241e")


@VisiData.api
def save_usv(vd, p, vs):
    vd.save_tsv(p, vs, row_delimiter="\u241e", delimiter="\u241f")
