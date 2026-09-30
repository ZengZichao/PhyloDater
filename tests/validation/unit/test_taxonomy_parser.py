"""Unit tests for TaxonomyParser."""

from phylodater.services import TaxonomyParser


class TestTaxonomyParser:
    def test_parse_embedded_format(self):
        parser = TaxonomyParser()
        result = parser.parse(
            "GB_GCA_001_d_Bacteria_p_Proteobacteria_c_Gammaproteobacteria"
        )
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Proteobacteria"
        assert result["class"] == "Gammaproteobacteria"

    def test_parse_table_format(self):
        parser = TaxonomyParser()
        result = parser.parse("d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia")
        assert result["domain"] == "Bacteria"
        assert result["phylum"] == "Cyanobacteriota"

    def test_load_from_file_two_column(self, taxonomy_table_path):
        parser = TaxonomyParser()
        count = parser.load_from_file(taxonomy_table_path)
        assert count == 2
        result = parser.parse("Sample1")
        assert result["domain"] == "Bacteria"

    def test_load_from_file_multi_column(self, fixtures_dir):
        path = fixtures_dir / "taxonomies" / "taxonomy_multicol.tsv"
        parser = TaxonomyParser()
        count = parser.load_from_file(path)
        assert count == 2

    def test_check_label_consistency(self, taxonomy_table_path):
        parser = TaxonomyParser()
        parser.load_from_file(taxonomy_table_path)
        matched, unmatched_tree, unmatched_table = parser.check_label_consistency(
            ["Sample1", "Sample2"]
        )
        assert len(matched) == 2
        assert not unmatched_tree
        assert not unmatched_table
