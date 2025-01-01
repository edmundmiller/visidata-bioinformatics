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
    """Sheet for displaying BED format data.
    
    Handles the Browser Extensible Data (BED) format which defines genomic regions.
    Supports required fields (chrom, start, end) and optional fields including:
    - name: Name/ID of the feature
    - score: Score from 0-1000
    - strand: + or - for strand
    - thickStart/End: Positions for thick drawing
    - itemRgb: Color in RGB or hex format
    - blockCount/Sizes/Starts: Exon/block structure
    """

    rowtype = "regions"  # rowdef: list of fields
    
    def __init__(self, name, source=None, **kwargs):
        super().__init__(name, source=source, delimiter="\t", headerlines=0, **kwargs)
        self.track_lines = []  # Store track lines for later reference
        
    def openRow(self, row):
        """Allow diving into track attributes when row is a track line"""
        if isinstance(row, str) and row.startswith('track'):
            return TrackAttributesSheet(name=f"track_attributes", source=row)
        return None

    def validate_blocks(self, start, end, block_count, block_sizes, block_starts):
        """Validate block coordinates are within feature bounds"""
        if not (block_sizes and block_starts):
            return 0, "", ""
            
        sizes = [int(x) for x in block_sizes.rstrip(',').split(',')]
        starts = [int(x) for x in block_starts.rstrip(',').split(',')]
        
        if len(sizes) != len(starts):
            return 0, "", ""
            
        # Validate each block
        valid_sizes = []
        valid_starts = []
        for size, rel_start in zip(sizes, starts):
            abs_start = start + rel_start
            abs_end = abs_start + size
            
            # Block must be within feature bounds
            if abs_start >= start and abs_end <= end:
                valid_sizes.append(str(size))
                valid_starts.append(str(rel_start))
                
        if not valid_sizes:
            return 0, "", ""
            
        return len(valid_sizes), ",".join(valid_sizes), ",".join(valid_starts)

    @asyncthread
    def reload(self):
        self.columns = []
        self.rows = []

        def make_getter(idx, type_func=str, validator=None):
            def getter(row):
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
            """Validate RGB string format (comma-separated or hex)"""
            try:
                if rgb.startswith('#'):
                    # Convert hex to RGB
                    rgb = rgb.lstrip('#')
                    r = int(rgb[0:2], 16)
                    g = int(rgb[2:4], 16)
                    b = int(rgb[4:6], 16)
                else:
                    r,g,b = map(int, rgb.split(','))
                return f"{min(max(r,0),255)},{min(max(g,0),255)},{min(max(b,0),255)}"
            except:
                return "0,0,0"

        # Required BED fields with validation
        self.addColumn(Column(name="chrom", type=str, getter=make_getter(0)))
        self.addColumn(Column(name="start", type=int, getter=make_getter(1, int)))
        self.addColumn(Column(name="end", type=int, getter=make_getter(2, int)))

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
            self.addColumn(Column(name=name, type=type_func, getter=make_getter(idx, type_func, validator)))

        # Load the data
        with self.source.open_text() as fp:
            lines = fp.readlines()
            for line in Progress(lines, "loading BED file"):
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "browser", "track")):
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
                            if block_count > 0 and len(fields) >= 12:
                                # Validate and clean block coordinates
                                valid_count, valid_sizes, valid_starts = self.validate_blocks(
                                    start, end, block_count,
                                    fields[10], fields[11]
                                )
                                fields[9] = valid_count
                                fields[10] = valid_sizes
                                fields[11] = valid_starts
                        except (ValueError, TypeError):
                            fields[9] = 0

                    self.addRow(fields)
                except Exception as e:
                    vd.warning(f"error parsing line: {line[:50]}... {str(e)}")

        # Register BED format detection
        vd.option('filetype', 'bed', 'BED', BedSheet)
