"""https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""
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

__version__ = "0.1"


class GffSheet(TsvSheet):
    """Sheet for displaying GFF (General Feature Format) data"""

    rowtype = "features"  # rowdef: list of values
    columns = [
        Column("seqid", 0, type=str),
        Column("source", 1, type=str),
        Column("type", 2, type=str),
        Column("start", 3, type=int),
        Column("end", 4, type=int),
        Column("score", 5, type=float),
        Column("strand", 6, type=str),
        Column("phase", 7, type=str),
        Column("attributes", 8, type=str),
    ]

    def iterload(self):
        for line in self.source:
            if line.startswith("#"):  # Skip comment/directive lines
                continue
            line = line.rstrip("\n")
            if not line:  # Skip empty lines
                continue
            row = line.split("\t")
            if len(row) >= 9:  # Only load valid rows with all 9 fields
                # Convert score to float or None if '.'
                if row[5] == ".":
                    row[5] = None
                yield row


@VisiData.api
def guess_gff(vd, p):
    """Guess if file is a GFF format based on content"""
    with p.open_text() as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:  # GFF requires exactly 9 fields
                try:
                    int(fields[3])  # start position
                    int(fields[4])  # end position
                    return dict(filetype='gff', _likelihood=9)
                except ValueError:
                    pass
            break
    return None


@VisiData.api
def open_gff(vd, p):
    return GffSheet(p.name, source=p)


# Add GFF format detection
VisiData.guess_handlers.append(guess_gff)
