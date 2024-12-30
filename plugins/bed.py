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
    AttrDict,
    ItemColumn,
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
    """Try to open as BED, fall back to TSV if parsing fails"""
    try:
        sheet = BedSheet(p.name, source=p)
        sheet.reload()
        if not sheet.rows:  # If no rows were successfully parsed
            return TsvSheet(p.name, source=p)
        return sheet
    except Exception as e:
        vd.warning(f"Failed to parse as BED ({str(e)}), falling back to TSV")
        return TsvSheet(p.name, source=p)


class TrackAttributesSheet(Sheet):
    """Sheet for displaying parsed track line attributes"""
    rowtype = "attributes"  # rowdef: AttrDict of key-value pairs
    
    def iterload(self):
        track_str = self.source
        if not track_str or not track_str.startswith('track'):
            return
            
        attrs = AttrDict()
        # Remove 'track' from the start
        parts = track_str[5:].strip().split()
        
        for part in parts:
            if '=' in part:
                key, value = part.split('=', 1)
                # Remove quotes if present
                value = value.strip('"\'')
                attrs[key.strip()] = value.strip()
            else:
                attrs[part.strip()] = ''
                
        yield attrs

    def addRow(self, row, index=None):
        super().addRow(row, index=index)
        
        # Add columns for any new keys
        for k in row:
            if not any(c.name == k for c in self.columns):
                self.addColumn(ItemColumn(k))


class BedSheet(TsvSheet):
    """Sheet for displaying BED format data"""

    rowtype = "regions"  # rowdef: list of fields
    
    def __init__(self, name, source=None, **kwargs):
        super().__init__(name, source=source, delimiter="\t", headerlines=0, **kwargs)
        self.track_lines = []  # Store track lines for later reference

    @asyncthread
    def reload(self):
        self.columns = []
        self.rows = []

        def make_getter(idx, type_func=str, validator=None):
            def getter(col, row):
                try:
                    val = type_func(row[idx]) if row and len(row) > idx else None
                    if validator and val is not None:
                        val = validator(val)
                    return val
                except (IndexError, ValueError, TypeError):
                    return None
            return getter

        def validate_score(score):
            """Validate and clamp score between 0-1000"""
            try:
                score = float(score)
                return min(max(score, 0), 1000)
            except (ValueError, TypeError):
                return 0

        def validate_strand(strand):
            """Validate strand is +, -, or ."""
            return strand if strand in ('+', '-', '.') else '.'

        def validate_rgb(rgb):
            """Validate RGB string format"""
            try:
                r,g,b = map(int, rgb.split(','))
                return f"{min(max(r,0),255)},{min(max(g,0),255)},{min(max(b,0),255)}"
            except:
                return "0,0,0"

        # Required BED fields with validation
        self.addColumn(Column(name="chrom", type=str, getter_type=make_getter(0)))
        self.addColumn(Column(name="start", type=int, getter_type=make_getter(1, int)))
        self.addColumn(Column(name="end", type=int, getter_type=make_getter(2, int)))

        # Optional BED fields with their types and validators
        optional_cols = [
            ("name", 3, str, None),  # Name of region
            ("score", 4, float, validate_score),  # Score from 0-1000
            ("strand", 5, str, validate_strand),  # + or - for strand
            ("thickStart", 6, int, None),  # Starting position at which feature is drawn thickly
            ("thickEnd", 7, int, None),  # Ending position at which feature is drawn thickly
            ("itemRgb", 8, str, validate_rgb),  # RGB color value (e.g., 255,0,0)
            ("blockCount", 9, int, None),  # Number of blocks (exons)
            ("blockSizes", 10, str, None),  # Comma-separated list of block sizes
            ("blockStarts", 11, str, None),  # Comma-separated list of block starts
        ]

        for name, idx, type_func, validator in optional_cols:
            self.addColumn(Column(name=name, getter=make_getter(idx, type_func, validator)))

        # Load the data
        with self.source.open_text() as fp:
            lines = fp.readlines()
            for line in Progress(lines, "loading BED file"):
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                if line.startswith("browser"):
                    # Store browser lines but don't add as rows
                    continue
                if line.startswith("track"):
                    self.track_lines.append(line)
                    continue
                try:
                    fields = line.split("\t")  # Explicitly split on tabs
                    # Ensure minimum 3 fields
                    if len(fields) < 3:
                        vd.warning(f"skipping line with too few fields: {line[:50]}...")
                        continue

                    # Validate required fields
                    if not fields[0].strip():
                        vd.warning(f"skipping line with empty chromosome: {line[:50]}...")
                        continue

                    # Convert coordinates - BED uses 0-based start and 1-based end
                    try:
                        start = int(fields[1])
                        end = int(fields[2])
                        if start < 0:
                            vd.warning(f"invalid negative start coordinate: {line[:50]}...")
                            continue
                        if end <= start:
                            vd.warning(f"end coordinate must be greater than start: {line[:50]}...")
                            continue
                        fields[1] = start
                        fields[2] = end
                    except ValueError as e:
                        vd.warning(f"invalid coordinates in line: {line[:50]}... {str(e)}")
                        continue

                    # Handle optional fields
                    if len(fields) >= 7:  # thickStart/End exist
                        try:
                            thick_start = int(fields[6])
                            thick_end = int(fields[7])
                            # Validate thick coordinates are within feature bounds
                            fields[6] = max(thick_start, start)
                            fields[7] = min(thick_end, end)
                        except (ValueError, TypeError):
                            fields[6] = start  # default to chromStart
                            fields[7] = end    # default to chromEnd

                    if len(fields) >= 10:  # blockCount/Sizes/Starts exist
                        try:
                            block_count = int(fields[9])
                            if block_count > 0:
                                # Validate block sizes and starts if present
                                if len(fields) >= 12:
                                    sizes = fields[10].rstrip(',').split(',')
                                    starts = fields[11].rstrip(',').split(',')
                                    if len(sizes) != block_count or len(starts) != block_count:
                                        vd.warning(f"mismatched block counts in line: {line[:50]}...")
                                        fields[9] = 0  # Reset block count if invalid
                        except (ValueError, TypeError):
                            fields[9] = 0

                    self.addRow(fields)
                except Exception as e:
                    vd.warning(f"error parsing line: {line[:50]}... {str(e)}")
