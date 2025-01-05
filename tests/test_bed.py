import pytest
from visidata import Path


def test_bed_loading(vd, load_plugin, test_data_dir):
    """Test basic BED file loading"""
    load_plugin("bed_pybedlite")

    bed_path = test_data_dir / "bed" / "ucsc_example1.bed"
    sheet = vd.openSource(Path(bed_path))

    assert sheet.__class__.__name__ == "BedPyblSheet"
    assert len(sheet.columns) >= 6
    assert [col.name for col in sheet.columns[:6]] == [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]


def test_bed_strand_colors(vd, load_plugin, test_data_dir):
    """Test BED strand coloring"""
    load_plugin("bed_pybedlite")

    bed_path = test_data_dir / "bed" / "ColorByStrandDemo.bed"
    sheet = vd.openSource(Path(bed_path))

    # Verify strand column exists
    strand_col = next(col for col in sheet.columns if col.name == "strand")
    assert strand_col is not None


def test_bed_to_gff_conversion(vd, load_plugin, test_data_dir):
    """Test BED to GFF conversion"""
    load_plugin("bed_pybedlite")
    load_plugin("gff")

    bed_path = test_data_dir / "bed" / "ucsc_example1.bed"
    sheet = vd.openSource(Path(bed_path))

    # Convert to GFF
    sheet.convert_to_gff(vd)

    # Verify conversion results
    gff_sheet = vd.sheets[1]
    assert gff_sheet.__class__.__name__ == "GffSheet"
    assert len(gff_sheet.columns) == 9
