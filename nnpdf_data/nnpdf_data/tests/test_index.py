"""
Tests for the dataset index module (nnpdf_data.cli.index).

Exercises the brace-expansion pattern matcher, the glob-based filter,
the sort-by-attribute function, and the full index builder using real datasets.
"""

from nnpdf_data.cli.index import _expand_pattern, build_index, filter_index, sort_index


def test_expand_pattern():
    """Test brace expansion (e.g. {ATLAS,CMS}* -> ['ATLAS*', 'CMS*'])."""
    # Simple glob (no braces) should pass through
    assert _expand_pattern("ATLAS*") == ["ATLAS*"]
    # Brace expansion should produce two patterns
    assert _expand_pattern("{ATLAS,CMS}*") == ["ATLAS*", "CMS*"]
    # Unmatched brace should fall back to the literal pattern
    assert _expand_pattern("unmatched{") == ["unmatched{"]


def test_filter_index(real_index):
    """Test that real datasets are correctly filtered by glob patterns."""
    # ATLAS filter should match multiple datasets
    atlas = filter_index(real_index, "ATLAS*")
    assert len(atlas) > 0
    # All matched datasets should have ATLAS experiment
    for name, ds in atlas.items():
        assert ds.experiment == "ATLAS" or name.startswith("ATLAS")

    # Brace expansion should match all ATLAS and CMS datasets
    atlas_cms = filter_index(real_index, "{ATLAS,CMS}*")
    assert len(atlas_cms) > len(atlas)  # should have more entries

    # Nonexistent pattern should match nothing
    assert len(filter_index(real_index, "NONEXISTENT_12345*")) == 0

    # Default wildcard should return all entries
    assert len(filter_index(real_index)) == len(real_index)


def test_sort_index(real_index):
    """Test that real entries are sorted by a given attribute key."""
    # Filter to a small subset for predictable sorting
    subset = {
        name: ds
        for name, ds in real_index.items()
        if ds.experiment in ("ATLAS", "CMS") and ds.process == "DY"
    }
    # Need at least 2 entries
    if len(subset) < 2:
        # Fallback to any ATLAS and CMS datasets
        atlas_entries = {
            name: ds for name, ds in list(real_index.items()) if ds.experiment == "ATLAS"
        }
        cms_entries = {name: ds for name, ds in list(real_index.items()) if ds.experiment == "CMS"}
        subset = {**dict(list(atlas_entries.items())[:1]), **dict(list(cms_entries.items())[:1])}

    result = sort_index(subset, "experiment")
    keys = list(result.keys())

    # ATLAS should come before CMS alphabetically
    experiments = [result[k].experiment for k in keys]
    assert experiments == sorted(experiments)


def test_build_index():
    """Test that the full index builder returns a sorted, deterministic dictionary."""
    idx = build_index()

    # Should return a non-empty dict
    assert isinstance(idx, dict) and len(idx) > 0
    # Keys should already be sorted
    assert list(idx.keys()) == sorted(idx.keys())
    # Multiple calls should be deterministic
    assert list(idx.keys()) == list(build_index().keys())
