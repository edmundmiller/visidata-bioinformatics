"""GFF3 file format loader for VisiData
See: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""
from visidata import (
    Sheet, 
    Column,
    Progress,
    options,
    vd,
    VisiData,
    asyncthread,
)

__version__ = "0.1"


class GffSheet(Sheet):
    """Sheet for displaying GFF (General Feature Format) data"""
    rowtype = "features"  # rowdef: list of [seqid, source, type, start, end, score, strand, phase, attrs]
    
    columns = [
        Column(name="seqid", getter=lambda c,r: r[0], type=str),
        Column(name="source", getter=lambda c,r: r[1], type=str), 
        Column(name="type", getter=lambda c,r: r[2], type=str),
        Column(name="start", getter=lambda c,r: int(r[3]) if r[3] != '.' else None, type=int),
        Column(name="end", getter=lambda c,r: int(r[4]) if r[4] != '.' else None, type=int),
        Column(name="score", getter=lambda c,r: float(r[5]) if r[5] != '.' else None, type=float),
        Column(name="strand", getter=lambda c,r: r[6] if r[6] != '.' else None, type=str),
        Column(name="phase", getter=lambda c,r: r[7] if r[7] != '.' else None, type=str),
        Column(name="attributes", getter=lambda c,r: r[8] if r[8] != '.' else None, type=str)
    ]

    def iterload(self):
        with self.source.open_text() as fp:
            lines = list(fp)
            for line in Progress(lines, 'loading'):
                try:
                    if line.startswith('#') or not line.strip():
                        continue
                        
                    fields = line.rstrip('\n').split('\t')
                    if not fields or len(fields) < 9:
                        vd.warning(f'skipping invalid line (expected 9 fields, got {len(fields)}): {line}')
                        continue
                    # Pad any missing fields with '.'
                    while len(fields) < 9:
                        fields.append('.')
                    yield fields
                except Exception as e:
                    vd.warning(f'error parsing line: {e}')


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
                    # Validate required integer fields
                    int(fields[3])  # start position
                    int(fields[4])  # end position
                    
                    # Check strand field
                    if fields[6] in ['+', '-', '.']:
                        return dict(filetype='gff', _likelihood=9)
                except ValueError:
                    pass
            break
    return None


@VisiData.api
def open_gff(vd, p):
    """Open a GFF file and return a GffSheet"""
    return GffSheet(p.name, source=p)

@VisiData.api
def save_gff(p, *sheets):
    """Save sheet to a GFF file"""
    if len(sheets) != 1:
        vd.fail("can only save one sheet to GFF")
    sheet = sheets[0]
    
    with p.open_text(mode='w') as fp:
        # Write format version
        fp.write("##gff-version 3\n")
        
        for row in Progress(sheet.rows, 'saving'):
            try:
                fields = [str(col.getTypedValue(row)) for col in sheet.visibleCols]
                if len(fields) != 9:
                    continue
                fp.write('\t'.join(fields) + '\n')
            except Exception as e:
                vd.warning(f'error saving row: {e}')
                
# Register GFF format detection
vd.filetype('gff', GffSheet)
vd.guess_handlers.append(guess_gff)
