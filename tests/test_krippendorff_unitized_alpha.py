from krippendorff_unitized_alpha.krippendorffs_unitized_alpha import (
    Segment,
    EMPTY_TAG,
    KrippendorffUnitizedAlpha,
    impute_empty_segments
)

import pytest
from math import isclose

TEST_DATA = {
    'coder1': [
        Segment(EMPTY_TAG, 0, 2),
        Segment("keyphrase", 3, 14),
        Segment(EMPTY_TAG, 15, 21),
        Segment("keyphrase", 22, 52),
        Segment(EMPTY_TAG, 53, 55)
    ],
    'coder2': [
        Segment(EMPTY_TAG, 0, 2),
        Segment("keyphrase", 3, 9),
        Segment(EMPTY_TAG, 10, 21),
        Segment("keyphrase", 22, 52),
        Segment(EMPTY_TAG, 53, 55)
    ],
    'coder3': [
        Segment(EMPTY_TAG, 0, 21),
        Segment("keyphrase", 22, 45),
        Segment(EMPTY_TAG, 46, 55)
    ],
}


TEST_DATA_WITHOUT_EMPTY_SEGMENTS = {
    'coder1': [
        Segment("keyphrase", 3, 14),
        Segment("keyphrase", 22, 52),
    ],
    'coder2': [
        Segment("keyphrase", 3, 9),
        Segment("keyphrase", 22, 52),
    ],
    'coder3': [
        Segment("keyphrase", 22, 45),
    ],
}


def test_works_as_expected():
    kua = KrippendorffUnitizedAlpha(TEST_DATA)
    expected = 0.56
    assert isclose(kua.result, expected, abs_tol=0.01)


@pytest.mark.parametrize("coder", TEST_DATA.keys())
def test_impute_empty_segments_works_as_expected(coder: str):
    test_text = "An example text about keyphrase identification coding..."
    result = impute_empty_segments(TEST_DATA_WITHOUT_EMPTY_SEGMENTS[coder], len(test_text))
    expected = TEST_DATA[coder]
    assert result == expected
