import pytest
from visidata import Path


def test_gff_loading(vd, load_plugin, test_data_dir):
    """Test basic GFF file loading"""
    load_plugin("gff")

    # Load test GFF file
    gff_path = test_data_dir / "gff" / "test.gff"
    sheet = vd.openSource(Path(gff_path))

    # Verify basic sheet properties
    assert sheet.__class__.__name__ == "GffSheet"
    assert len(sheet.columns) == 9
    assert [col.name for col in sheet.columns] == [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]


def test_gff_to_bed_conversion(vd, load_plugin, test_data_dir):
    """Test GFF to BED conversion"""
    load_plugin("gff")
    load_plugin("bed_pybedlite")

    gff_path = test_data_dir / "gff" / "test.gff"
    sheet = vd.openSource(Path(gff_path))

    # Convert to BED
    sheet.convert_to_bed(vd)

    # Verify conversion results
    bed_sheet = vd.sheets[1]  # Converted sheet should be pushed
    assert bed_sheet.__class__.__name__ == "BedPyblSheet"
    assert len(bed_sheet.columns) >= 6  # BED6 format minimum


def test_gff_attributes_parsing(vd, load_plugin, test_data_dir):
    """Test GFF attributes parsing"""
    load_plugin("gff")

    gff_path = test_data_dir / "gff" / "test.gff"
    sheet = vd.openSource(Path(gff_path))

    # Get a row with attributes
    row = next(r for r in sheet.rows if r[8] != ".")

    # Open attributes sheet
    attrs_sheet = sheet.openRow(row)
    assert attrs_sheet.__class__.__name__ == "AttributesSheet"
    assert len(attrs_sheet.columns) == 2
    assert [col.name for col in attrs_sheet.columns] == ["attribute", "value"]
