"""Retain paired identities for counted singleton secondary STAR alignments."""
import sys
from collections import Counter
from pathlib import Path

import pysam
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from bulk_rna import assigned_pair


@pytest.mark.parametrize("selected_r1", [True, False])
@pytest.mark.parametrize("mate_reverse", [True, False])
@pytest.mark.parametrize("unmapped_primary", ["same", "opposite", "neither"])
def test_singleton_hit_recovers_mate_from_other_primary_alignment(selected_r1, mate_reverse, unmapped_primary):
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chr1", LN=1000)]))

    def record(r1, secondary, position, gene, reverse=False):
        r = pysam.AlignedSegment(header)
        r.query_name = "same-pair"
        r.flag = 1 | (64 if r1 else 128) | (256 if secondary else 0) | (16 if reverse else 0)
        r.reference_id = r.next_reference_id = 0
        r.reference_start = position
        r.next_reference_start = 300
        r.query_sequence = "AACCGT" if not reverse else "ACGGTT"
        r.query_qualities = [10, 11, 12, 13, 14, 15] if not reverse else [15, 14, 13, 12, 11, 10]
        r.cigarstring = "6M"
        r.set_tags([("GX", gene), ("CB", "physical-pool"), ("XP", "physical-pool"),
                    ("UB", "ACGT"), ("UM", "ATGT"), ("SB", "sample"), ("XI", "library-cell"),
                    ("ZI", 42), ("NH", 2), ("HI", 1 if secondary else 2)])
        return r

    primary1 = record(True, False, 300, "-", mate_reverse if not selected_r1 else False)
    primary2 = record(False, False, 320, "-", mate_reverse if selected_r1 else False)
    if unmapped_primary != "neither":
        primary_is_r1 = selected_r1 if unmapped_primary == "same" else not selected_r1
        primary = primary1 if primary_is_r1 else primary2
        primary.is_unmapped = True
        primary.reference_id = primary.reference_start = -1
        primary.cigarstring = None
    selected = record(selected_r1, True, 100, "gene")
    selected.mate_is_unmapped = True
    selected.next_reference_id = selected.next_reference_start = -1
    key, pair, complete = assigned_pair([primary1, primary2, selected], None, Counter())
    assert key == ("gene", "physical-pool", "ACGT") and not complete
    assert len(pair) == 2 and {r.is_read1 for r in pair} == {True, False}
    mapped, unmapped = pair
    assert not mapped.is_secondary and not unmapped.is_secondary
    assert mapped.reference_start == 100 and unmapped.is_unmapped
    assert mapped.mate_is_unmapped and not unmapped.mate_is_unmapped
    assert unmapped.query_sequence == "AACCGT"
    assert list(unmapped.query_qualities) == [10, 11, 12, 13, 14, 15]
    assert unmapped.next_reference_start == 100 and unmapped.reference_start == 100
    assert all(r.get_tag("ZI") == 42 and r.get_tag("XP") == "physical-pool" for r in pair)
    assert mapped.query_name == unmapped.query_name
    assert unmapped.get_tag("UM") == "ATGT" and unmapped.get_tag("UB") == "ACGT"

    # A truly missing input mate still fails explicitly.
    with pytest.raises(ValueError, match="missing mate"):
        assigned_pair([selected], None, Counter())
