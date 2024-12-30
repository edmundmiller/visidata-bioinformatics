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


# Add GFF format with .gff and .gff3 extensions
vd.save_filetype("gff", "tsv")
vd.save_filetype("gff3", "tsv")


def open_gff(p):
    return GffSheet(p.name, source=p)


# Register the GFF opener
vd.openfile_extensions.append("gff")
vd.openfile_extensions.append("gff3")
