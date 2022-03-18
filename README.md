# Krippendorff Unitized Alpha

This calculates the Krippendorff Unitized Alpha Inter Annotator Agreement.
https://link.springer.com/article/10.1007/s11135-015-0266-1

## Usage

```python
from krippendorff_unitized_alpha import *
data = {
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

kua = KrippendorffUnitizedAlpha(data)
print(kua.result)
# >>> 0.5671992508994135
```

If you don't have the *empty segments* then you can use a utility function to impute/compute them.

```python

from krippendorff_unitized_alpha import *

test_text = "An example text about keyphrase identification coding..."
data_without_empty_segments = {
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

kua = KrippendorffUnitizedAlpha(
    {
        annotator: impute_empty_segments(annotations, len(test_text))
        for annotator, annotations in data_without_empty_segments.items()
    }
)
print(kua.result)
# >>> 0.5671992508994135
```

