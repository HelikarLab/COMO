import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

from como.peak_finder import (
    _encode_signs,
    _enforce_minimum_peak_distance,
    _validate_args,
    find_peaks,
)


class TestValidateArgs:
    def test_multidimensional_array_raises_error(self):
        x: npt.NDArray[float] = np.array([[1, 2], [3, 4]])
        with pytest.raises(ValueError, match="Expected a 1D array, got 2D array instead"):
            _validate_args(x, 1, 1, "0", 0.0, 1, 0.0)

    def test_nan_values_raise_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, np.nan, 4.0])
        with pytest.raises(ValueError, match="Input x contains NaNs"):
            _validate_args(x, 1, 1, "0", 0.0, 1, 0.0)

    def test_negative_nups_raises_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError, match="Argument 'nups' must be non-negative"):
            _validate_args(x, -1, 1, "0", 0.0, 1, 0.0)

    def test_negative_ndowns_raises_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError, match="Argument 'ndowns' must be non-negative"):
            _validate_args(x, 1, -1, "0", 0.0, 1, 0.0)

    def test_invalid_zero_raises_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError, match="Argument 'zero' must be '0', '\\+', or '-'"):
            _validate_args(x, 1, 1, "x", 0.0, 1, 0.0)  # type: ignore[invalid-type-argument]

    def test_negative_min_peak_height_raises_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError, match="Argument 'min_peak_height' must be non-negative"):
            _validate_args(x, 1, 1, "0", -1.0, 1, 0.0)

    def test_negative_min_peak_distance_raises_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError, match="Argument 'minpeakdistance' must be non-negative"):
            _validate_args(x, 1, 1, "0", 0.0, -1, 0.0)

    def test_negative_threshold_raises_error(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        with pytest.raises(ValueError, match="Argument 'threshold' must be non-negative"):
            _validate_args(x, 1, 1, "0", 0.0, 1, -1.0)

    def test_valid_args_does_not_raise(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        _validate_args(x, 1, 1, "0", 0.0, 1, 0.0)


class TestEncodeSigns:
    def test_increasing_sequence(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0])
        result: str = _encode_signs(x, "0")
        assert result == "+++"

    def test_decreasing_sequence(self):
        x: npt.NDArray[float] = np.array([4.0, 3.0, 2.0, 1.0])
        result: str = _encode_signs(x, "0")
        assert result == "---"

    def test_flat_sequence_with_zero(self):
        x: npt.NDArray[float] = np.array([1.0, 1.0, 1.0])
        result: str = _encode_signs(x, "0")
        assert result == "00"

    def test_flat_sequence_with_plus(self):
        x: npt.NDArray[float] = np.array([1.0, 1.0, 1.0])
        result: str = _encode_signs(x, "+")
        assert result == "++"

    def test_flat_sequence_with_minus(self):
        x: npt.NDArray[float] = np.array([1.0, 1.0, 1.0])
        result: str = _encode_signs(x, "-")
        assert result == "--"

    def test_flat_sequence_with_dollarsign(self):
        x: npt.NDArray[float] = np.array([1.0, 1.0, 1.0])
        result: str = _encode_signs(x, "$")
        assert result == "$$"

    def test_mixed_sequence(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 2.0, 3.0, 2.0])
        result: str = _encode_signs(x, "0")
        assert result == "+0+-"


class TestEnforceMinimumPeakDistance:
    def test_inplace_removes_close_peaks(self):
        df: pd.DataFrame = pd.DataFrame(
            {
                "height": [10.0, 8.0, 5.0],
                "peak_idx": [0, 2, 10],
                "start_idx": [0, 1, 9],
                "end_idx": [1, 3, 11],
            }
        )
        _enforce_minimum_peak_distance(df, min_peak_distance=5, inplace=True)
        assert len(df) == 2
        assert df["peak_idx"].tolist() == [0, 10]

    def test_not_inplace_returns_new_dataframe(self):
        df: pd.DataFrame = pd.DataFrame(
            {
                "height": [10.0, 8.0, 5.0],
                "peak_idx": [0, 2, 10],
                "start_idx": [0, 1, 9],
                "end_idx": [1, 3, 11],
            }
        )
        result = _enforce_minimum_peak_distance(df, min_peak_distance=5, inplace=False)
        assert len(result) == 2
        assert len(df) == 3  # original unchanged
        assert result["peak_idx"].tolist() == [0, 10]

    def test_all_peaks_sufficiently_spaced(self):
        df: pd.DataFrame = pd.DataFrame(
            {
                "height": [10.0, 8.0],
                "peak_idx": [0, 10],
                "start_idx": [0, 9],
                "end_idx": [1, 11],
            }
        )
        _enforce_minimum_peak_distance(df, min_peak_distance=5, inplace=True)
        assert len(df) == 2


class TestFindPeaks:
    def test_simple_peak(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, min_peak_height=0.0)
        assert len(result) == 1
        assert result.iloc[0]["peak_idx"] == 3
        assert result.iloc[0]["height"] == 3.0

    def test_multiple_peaks(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, min_peak_height=0.0)
        assert len(result) == 3

    def test_no_peaks(self):
        x: npt.NDArray[float] = np.asarray([1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)
        result: pd.DataFrame = find_peaks(x)
        assert len(result) == 0
        assert list(result.columns) == ["height", "peak_idx", "start_idx", "end_idx"]

    def test_min_peak_height_filter(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 0.0, 5.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, min_peak_height=2.0)
        assert len(result) == 1
        assert result.iloc[0]["height"] == 5.0

    def test_nups_and_ndowns(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, nups=2, ndowns=2, min_peak_height=0.0)
        assert len(result) == 1

    def test_npeaks_limits_output(self):
        x: npt.NDArray[float] = np.asarray([0.0, 5.0, 0.0, 3.0, 0.0, 1.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, npeaks=2, min_peak_height=0.0)
        assert len(result) == 2

    def test_sortstr_orders_by_height(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 0.0, 5.0, 0.0, 3.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, sortstr=True, min_peak_height=0.0)
        heights = result["height"].tolist()
        assert heights == sorted(heights, reverse=True)

    def test_min_peak_distance(self):
        x: npt.NDArray[float] = np.asarray([0.0, 10.0, 0.0, 8.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, min_peak_distance=3, min_peak_height=0.0)
        assert len(result) == 1
        assert result.iloc[0]["height"] == 10.0

    def test_threshold_filter(self):
        x: npt.NDArray[float] = np.asarray([1.0, 2.0, 1.5, 5.0, 1.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, threshold=2.0, min_peak_height=0.0)
        assert len(result) == 1
        assert result.iloc[0]["height"] == 5.0

    def test_accepts_list_input(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 2.0, 1.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, min_peak_height=0.0)
        assert len(result) == 1

    def test_accepts_numpy_array(self):
        x: npt.NDArray[float] = np.array([0.0, 1.0, 2.0, 1.0, 0.0])
        result: pd.DataFrame = find_peaks(x, min_peak_height=0.0)
        assert len(result) == 1

    def test_custom_peak_pattern(self):
        x: npt.NDArray[float] = np.asarray([0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0], dtype=float)
        result: pd.DataFrame = find_peaks(x, peak_pattern=r"[+]{3,}[-]{3,}", min_peak_height=0.0)
        assert len(result) == 1
