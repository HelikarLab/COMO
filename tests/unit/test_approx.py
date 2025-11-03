import numpy as np
import pytest

from como.approx import _coerce_to_float_array, _regularize_values, approx


class TestCoerceToFloatArray:
    """Tests for the _coerce_to_float_array helper function."""

    def test_converts_int_list_to_array(self):
        result = _coerce_to_float_array([1, 2, 3])
        assert isinstance(result, np.ndarray)
        assert result.dtype == float
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])

    def test_converts_float_list_to_array(self):
        result = _coerce_to_float_array([1.0, 2.0, 3.0])
        assert isinstance(result, np.ndarray)
        assert result.dtype == float
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])


class TestRegularizeValues:
    """Tests for the _regularize_values function."""

    def test_removes_na_pairs_when_na_rm_true(self):
        x: list[float] = np.array([1.0, 2.0, np.nan, 4.0])
        y: list[float] = np.array([10.0, np.nan, 30.0, 40.0])
        result = _regularize_values(x, y, na_rm=True, ties="mean")
        np.testing.assert_array_equal(result.x, [1.0, 4.0])
        np.testing.assert_array_equal(result.y, [10.0, 40.0])

    def test_raises_error_with_na_in_x_when_na_rm_false(self):
        x: list[float] = np.array([1.0, np.nan, 3.0])
        y: list[float] = np.array([10.0, 20.0, 30.0])
        with pytest.raises(ValueError, match="NA values in x are not allowed"):
            _regularize_values(x, y, na_rm=False, ties="mean")

    def test_sorts_by_x(self):
        x: list[float] = np.array([3.0, 1.0, 2.0])
        y: list[float] = np.array([30.0, 10.0, 20.0])
        result = _regularize_values(x, y, na_rm=True, ties="mean")
        np.testing.assert_array_equal(result.x, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(result.y, [10.0, 20.0, 30.0])

    def test_aggregates_duplicates_with_mean(self):
        x: list[float] = np.array([1.0, 1.0, 2.0])
        y: list[float] = np.array([10.0, 20.0, 30.0])
        result = _regularize_values(x, y, na_rm=True, ties="mean")
        np.testing.assert_array_equal(result.x, [1.0, 2.0])
        np.testing.assert_array_equal(result.y, [15.0, 30.0])

    def test_aggregates_duplicates_with_first(self):
        x: list[float] = np.array([1.0, 1.0, 2.0])
        y: list[float] = np.array([10.0, 20.0, 30.0])
        result = _regularize_values(x, y, na_rm=True, ties="first")
        np.testing.assert_array_equal(result.x, [1.0, 2.0])
        np.testing.assert_array_equal(result.y, [10.0, 30.0])

    def test_aggregates_duplicates_with_last(self):
        x: list[float] = np.array([1.0, 1.0, 2.0])
        y: list[float] = np.array([10.0, 20.0, 30.0])
        result = _regularize_values(x, y, na_rm=True, ties="last")
        np.testing.assert_array_equal(result.x, [1.0, 2.0])
        np.testing.assert_array_equal(result.y, [20.0, 30.0])

    def test_handles_empty_arrays(self):
        x: list[float] = np.array([])
        y: list[float] = np.array([])
        result = _regularize_values(x, y, na_rm=True, ties="mean")
        assert result.x.size == 0
        assert result.y.size == 0
        assert result.not_na.size == 0

    def test_callable_ties_function(self):
        x: list[float] = np.array([1.0, 1.0, 2.0])
        y: list[float] = np.array([10.0, 20.0, 30.0])
        result = _regularize_values(x, y, na_rm=True, ties=np.sum)
        np.testing.assert_array_equal(result.x, [1.0, 2.0])
        np.testing.assert_array_equal(result.y, [30.0, 30.0])


class TestApprox:
    """Tests for the main approx function."""

    def test_basic_linear_interpolation(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [1.5, 2.5]
        result = approx(x, y, xout=xout)
        np.testing.assert_array_almost_equal(result.y, [15.0, 25.0])

    def test_one_argument_form(self):
        y: list[float] = [10.0, 20.0, 30.0]
        result = approx(y, xout=[1.5, 2.5])
        np.testing.assert_array_almost_equal(result.y, [15.0, 25.0])

    def test_default_n_points(self):
        x: list[float] = [1.0, 5.0]
        y: list[float] = [10.0, 50.0]
        result = approx(x, y)
        assert len(result.x) == 50
        assert len(result.y) == 50

    def test_custom_n_points(self):
        x: list[float] = [1.0, 5.0]
        y: list[float] = [10.0, 50.0]
        result = approx(x, y, n=10)
        assert len(result.x) == 10

    def test_extrapolation_rule_1(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [0.5, 3.5]
        result = approx(x, y, xout=xout, rule=1)
        assert np.isnan(result.y[0])
        assert np.isnan(result.y[1])

    def test_extrapolation_rule_2(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [0.5, 3.5]
        result = approx(x, y, xout=xout, rule=2)
        assert result.y[0] == 10.0
        assert result.y[1] == 30.0

    def test_yleft_yright_parameters(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [0.5, 3.5]
        result = approx(x, y, xout=xout, yleft=5.0, yright=35.0)
        assert result.y[0] == 5.0
        assert result.y[1] == 35.0

    def test_constant_method_f_0(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [1.5]
        result = approx(x, y, xout=xout, method="constant", f=0.0)
        assert result.y[0] == 10.0

    def test_constant_method_f_1(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [1.5]
        result = approx(x, y, xout=xout, method="constant", f=1.0)
        assert result.y[0] == 20.0

    def test_constant_method_f_05(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [1.5]
        result = approx(x, y, xout=xout, method="constant", f=0.5)
        assert result.y[0] == 15.0

    def test_exact_match_returns_exact_value(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [2.0]
        result = approx(x, y, xout=xout)
        assert result.y[0] == 20.0

    def test_na_in_xout_returns_nan(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [np.nan, 2.0]
        result = approx(x, y, xout=xout)
        assert np.isnan(result.y[0])
        assert result.y[1] == 20.0

    def test_method_numeric_codes(self):
        x: list[float] = [1.0, 2.0]
        y: list[float] = [10.0, 20.0]
        xout: list[float] = [1.5]
        result1 = approx(x, y, xout=xout, method=1)
        result2 = approx(x, y, xout=xout, method=2, f=0.0)
        assert result1.y[0] == 15.0
        assert result2.y[0] == 10.0

    def test_raises_error_different_length_xy(self):
        x: list[float] = [1.0, 2.0]
        y: list[float] = [10.0]
        with pytest.raises(ValueError, match="x and y must have same length"):
            approx(x, y)

    def test_raises_error_invalid_method(self):
        x: list[float] = [1.0, 2.0]
        y: list[float] = [10.0, 20.0]
        with pytest.raises(ValueError, match="invalid interpolation method"):
            approx(x, y, method="invalid")

    def test_raises_error_invalid_method_code(self):
        x: list[float] = [1.0, 2.0]
        y: list[float] = [10.0, 20.0]
        with pytest.raises(ValueError, match="invalid interpolation method"):
            approx(x, y, method=3)

    def test_raises_error_invalid_f(self):
        x: list[float] = [1.0, 2.0]
        y: list[float] = [10.0, 20.0]
        with pytest.raises(ValueError, match="invalid f value"):
            approx(x, y, method="constant", f=2.0)

    def test_raises_error_need_two_points_for_linear(self):
        x: list[float] = [1.0]
        y: list[float] = [10.0]
        with pytest.raises(ValueError, match="need at least two non-NA values"):
            approx(x, y, method="linear")

    def test_raises_error_zero_points(self):
        x: list[float] = []
        y: list[float] = []
        with pytest.raises(ValueError, match="zero non-NA points"):
            approx(x, y, method="constant")

    def test_handles_ties_mean(self):
        x: list[float] = [1.0, 1.0, 2.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [1.5]
        result = approx(x, y, xout=xout, ties="mean")
        assert result.y[0] == 22.5

    def test_rule_as_list(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, 20.0, 30.0]
        xout: list[float] = [0.5, 3.5]
        result = approx(x, y, xout=xout, rule=[1, 2])
        assert np.isnan(result.y[0])
        assert result.y[1] == 30.0

    def test_na_rm_false_with_na_in_y(self):
        x: list[float] = [1.0, 2.0, 3.0]
        y: list[float] = [10.0, np.nan, 30.0]
        xout: list[float] = [2.5]
        result = approx(x, y, xout=xout, na_rm=False)
        # After filtering out NA, should interpolate between 1.0->10.0 and 3.0->30.0
        assert result.y[0] == 25.0
