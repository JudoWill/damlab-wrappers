"""Integration tests for hiv/vadr-genbank.

Verifies that version suffixes are stripped correctly from the LOCUS /
ACCESSION / VERSION lines in GenBank files and from FASTA headers, while all
other content is preserved unchanged.
"""

import re
from pathlib import Path

ACCESSION  = "K03455"
TEST_ROOT  = Path(__file__).parent
OUT_DIR    = TEST_ROOT / "test_output"
OUT_GB     = OUT_DIR / f"{ACCESSION}.gb"
OUT_FASTA  = OUT_DIR / f"{ACCESSION}.fasta"


class TestGenBankNormalisation:
    def test_output_genbank_exists(self):
        assert OUT_GB.exists(), f"Normalised GenBank not found: {OUT_GB}"

    def test_locus_name_is_accession(self):
        """LOCUS name token must equal the bare accession (not the legacy name)."""
        for line in OUT_GB.read_text().splitlines():
            if line.startswith("LOCUS"):
                locus_name = line.split()[1]
                assert locus_name == ACCESSION, (
                    f"LOCUS name was not normalised to accession: {line!r}"
                )
                return
        raise AssertionError("No LOCUS line found in normalised GenBank file")

    def test_version_line_no_suffix(self):
        """VERSION line must contain the bare accession with no .N suffix."""
        for line in OUT_GB.read_text().splitlines():
            if line.startswith("VERSION"):
                assert re.search(r"\bK03455\.\d+", line) is None, (
                    f"VERSION line still contains a version suffix: {line!r}"
                )
                assert ACCESSION in line, (
                    f"VERSION line does not contain accession: {line!r}"
                )
                return
        raise AssertionError("No VERSION line found in normalised GenBank file")

    def test_features_and_sequence_preserved(self):
        """Content outside LOCUS and VERSION lines must not be altered."""
        original   = (
            Path(__file__).parent.parent.parent / "vadr-test" / f"{ACCESSION}.gb"
        ).read_text()
        normalised = OUT_GB.read_text()

        def drop_locus_and_version(text):
            return "\n".join(
                l for l in text.splitlines()
                if not l.startswith("LOCUS") and not l.startswith("VERSION")
            )

        assert drop_locus_and_version(original) == drop_locus_and_version(normalised), (
            "Content outside LOCUS/VERSION lines was unexpectedly modified"
        )


class TestFastaCopy:
    """FASTA must be copied unchanged — v-build.pl requires the .N version suffix."""

    def test_output_fasta_exists(self):
        assert OUT_FASTA.exists(), f"FASTA output not found: {OUT_FASTA}"

    def test_fasta_identical_to_input(self):
        """The output FASTA must be byte-for-byte identical to the input."""
        original = (
            Path(__file__).parent.parent.parent / "vadr-test" / f"{ACCESSION}.fasta"
        ).read_text()
        output = OUT_FASTA.read_text()
        assert original == output, (
            "FASTA was modified — it must be copied unchanged so that "
            "v-build.pl can find the required accession.version header"
        )

    def test_header_retains_version_suffix(self):
        """Sanity check: the header still has the .N suffix."""
        for line in OUT_FASTA.read_text().splitlines():
            if line.startswith(">"):
                token = line[1:].split()[0]
                assert re.match(r"K03455\.\d+$", token), (
                    f"FASTA header is missing required version suffix: {line!r}"
                )
                return
        raise AssertionError("No FASTA header found in output FASTA")
