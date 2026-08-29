"""
Tests for chart construction in app/core/plots.py.

Wiring pairwise alignment into the UI (see tests/test_alignment.py) finally
exposed views that had never been rendered end to end. Screenshotting the
populated analysis page against real sample data surfaced three real defects
in the charts themselves, all pinned here so they can't come back silently:

1. create_composition_plot / create_kmer_plot: outside-positioned bar labels
   (e.g. "25.3%") were clipped by the plot's top edge because the y-axis had
   no headroom reserved above the tallest bar.
2. create_orf_plot: frame labels like "+1" / "-1" parse as numbers, so
   Plotly's default axis-type inference plotted the six reading frames on a
   continuous number line instead of as six discrete categories.
3. create_length_distribution / create_orf_plot's length histogram: when
   every sequence (or ORF) is exactly the same length — the common case for
   fixed-cycle sequencer reads — the axis autoranged so tightly around that
   one value that it labeled ticks in fractional base pairs ("39.6", "40.2").
"""

from app.core.orf import ORF
from app.core.plots import create_composition_plot, create_kmer_plot, create_orf_plot, create_length_distribution


def _make_orf(length_nt: int, frame: int, strand: str) -> ORF:
    return ORF(
        start=0, end=length_nt, length_nt=length_nt, length_aa=length_nt // 3 - 1,
        frame=frame, strand=strand, start_codon="ATG", stop_codon="TAA",
        sequence="A" * length_nt, protein="M" * (length_nt // 3 - 1), gc_content=50.0,
    )


class TestBarLabelHeadroom:
    """An outside text label must have room to render above its bar."""

    def test_composition_plot_leaves_headroom_above_tallest_bar(self):
        composition = {"A": 253, "T": 247, "C": 247, "G": 252, "N": 1}
        fig = create_composition_plot(composition)

        percentages = fig.data[0].y
        y_range = fig.layout.yaxis.range
        assert y_range is not None
        assert y_range[1] > max(percentages), (
            "y-axis must extend past the tallest bar or its outside label clips"
        )

    def test_kmer_plot_leaves_headroom_above_tallest_bar(self):
        counts = {f"kmer{i}": v for i, v in enumerate([81, 70, 69, 68, 65])}
        fig = create_kmer_plot(counts, k=3)

        y_range = fig.layout.yaxis.range
        assert y_range is not None
        assert y_range[1] > 81


class TestOrfFramesAreCategorical:
    """Reading frames ('+1', '-2', ...) are six discrete bins, not a number line."""

    def test_frame_axis_is_category_typed(self):
        orfs = [_make_orf(150, 1, "+"), _make_orf(120, 2, "+"), _make_orf(90, 2, "-")]
        fig = create_orf_plot(orfs)

        # Second subplot ("ORFs by Frame") is xaxis2 in a 1x2 make_subplots figure.
        assert fig.layout.xaxis2.type == "category"

    def test_frame_labels_are_biologically_ordered(self):
        # Deliberately constructed out of frame order to check the plot
        # doesn't just echo whatever order Counter first saw them in.
        orfs = [_make_orf(90, 2, "-"), _make_orf(150, 1, "+"), _make_orf(120, 2, "+")]
        fig = create_orf_plot(orfs)

        frame_trace = fig.data[1]
        assert list(frame_trace.x) == ["+1", "+2", "-2"]


class TestDegenerateLengthAxis:
    """Uniform-length input must not autorange into fractional-bp ticks."""

    def test_uniform_lengths_get_whole_bp_ticks(self):
        fig = create_length_distribution([40, 40, 40, 40, 40])

        assert fig.layout.xaxis.dtick == 1
        x_range = fig.layout.xaxis.range
        assert x_range is not None
        assert x_range[1] - x_range[0] > 1, "a single-bp-wide window is what produced fractional ticks"

    def test_varied_lengths_are_left_to_autorange(self):
        # The non-degenerate case must keep its normal nice round-number
        # ticks (320, 340, ...) rather than always forcing dtick=1.
        fig = create_length_distribution([320, 400, 405, 406, 410, 480])

        assert fig.layout.xaxis.dtick is None

    def test_orf_length_subplot_also_gets_whole_bp_ticks_when_uniform(self):
        orfs = [_make_orf(120, f, "+") for f in (1, 2, 3)]
        fig = create_orf_plot(orfs)

        assert fig.layout.xaxis.dtick == 1
