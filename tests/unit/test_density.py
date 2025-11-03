from typing import Literal, cast

import numpy as np
import numpy.typing as npt
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from como.density import DensityResult, bin_distance, density, dnorm, nrd0

KERNEL_TYPE = Literal["gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"]


class TestBinDistance:
    def test_basic_binning(self):
        x: npt.NDArray[float] = np.array([0.5, 1.5, 2.5])
        weights: npt.NDArray[float] = np.array([1.0, 1.0, 1.0])
        result: npt.NDArray[float] = bin_distance(x, weights, lo=0, up=3, n=4)

        assert result.shape == (8,)
        assert np.all(result >= 0)

    def test_weighted_binning(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0])
        weights: npt.NDArray[float] = np.array([2.0, 1.0])
        result: npt.NDArray[float] = bin_distance(x, weights, lo=0, up=3, n=4)

        assert result.shape == (8,)
        assert result.sum() > 0

    def test_empty_array(self):
        x: npt.NDArray[float] = np.array([])
        weights: npt.NDArray[float] = np.array([])
        result: npt.NDArray[float] = bin_distance(x, weights, lo=0, up=1, n=2)

        assert result.shape == (4,)
        assert_array_equal(result, np.zeros(4))

    def test_out_of_bounds_handling(self):
        x: npt.NDArray[float] = np.array([10.0])
        weights: npt.NDArray[float] = np.array([1.0])
        result: npt.NDArray[float] = bin_distance(x, weights, lo=0, up=1, n=2)

        assert result.shape == (4,)


class TestDnorm:
    def test_standard_normal(self):
        # dnorm(0) for standard normal should be ~0.3989
        result: float = dnorm(0.0, mean=0.0, sd=1.0)
        expected: float = 1.0 / np.sqrt(2 * np.pi)
        assert_allclose(result, expected, rtol=1e-10)

    def test_with_mean_and_sd(self):
        result: float = dnorm(5.0, mean=5.0, sd=2.0)
        expected: float = 1.0 / (2.0 * np.sqrt(2 * np.pi))
        assert_allclose(result, expected, rtol=1e-10)

    def test_log_density(self):
        result: float = dnorm(0.0, mean=0.0, sd=1.0, log=True)
        expected: float = np.log(1.0 / np.sqrt(2 * np.pi))
        assert_allclose(result, expected, rtol=1e-10)

    def test_nan_inputs(self):
        assert np.isnan(dnorm(np.nan))
        assert np.isnan(dnorm(0.0, mean=np.nan))
        assert np.isnan(dnorm(0.0, sd=np.nan))

    def test_negative_sd(self):
        result: float = dnorm(0.0, sd=-1.0)
        assert np.isnan(result)

    def test_infinite_sd(self):
        result: float = dnorm(0.0, sd=np.inf)
        assert result == 0.0 or result == -np.inf

    def test_zero_sd_at_mean(self):
        result: float = dnorm(5.0, mean=5.0, sd=0.0)
        assert np.isinf(result)

    def test_zero_sd_away_from_mean(self):
        result: float = dnorm(5.0, mean=3.0, sd=0.0, log=False)
        assert result == 0.0

    def test_fast_dnorm_flag(self):
        result_fast: float = dnorm(1.0, fast_dnorm=True)
        result_slow: float = dnorm(1.0, fast_dnorm=False)
        assert_allclose(result_fast, result_slow, rtol=1e-5)

    def test_large_values(self):
        result: float = dnorm(100.0, mean=0.0, sd=1.0)
        assert result >= 0.0
        assert result < 1e-10


class TestNrd0:
    def test_basic_calculation(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result: float = nrd0(x)
        assert result > 0
        assert np.isfinite(result)

    def test_with_constant_values(self):
        x: npt.NDArray[float] = np.array([5.0, 5.0, 5.0, 5.0])
        result: float = nrd0(x)
        assert result > 0
        assert np.isfinite(result)

    def test_single_nonzero_value(self):
        x: npt.NDArray[float] = np.array([0.0, 7.0])
        result: float = nrd0(x)
        assert result > 0

    def test_all_zeros(self):
        # if the input array is changed from a shape of `(3,)`, the result will change
        # this is because `nrd0` takes the length of the input into account when calculating bandwidth
        x: npt.NDArray[float] = np.array([0.0, 0.0, 0.0], dtype=float)
        result: float = nrd0(x)
        assert result == 0.7224674055842076

    def test_insufficient_data(self):
        with pytest.raises(ValueError, match="need at least 2 data points"):
            nrd0(np.array([1.0]))

    def test_empty_array(self):
        with pytest.raises(ValueError, match="need at least 2 data points"):
            nrd0(np.array([]))


class TestDensity:
    def test_basic_density(self):
        x: npt.NDArray[float] = np.random.normal(0, 1, 100)
        result: DensityResult = density(x)

        assert isinstance(result, DensityResult)
        assert len(result.x) == 512
        assert len(result.y) == 512
        assert result.bw > 0
        assert np.all(result.y >= 0)

    def test_custom_bandwidth(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result: DensityResult = density(x, bw=0.5)

        assert_allclose(result.bw, 0.5, rtol=1e-10)

    def test_with_weights(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        weights: npt.NDArray[float] = np.array([0.1, 0.2, 0.4, 0.2, 0.1])
        result: DensityResult = density(x, weights=weights, n=5)

        assert len(result.x) == 5
        assert np.all(result.y >= 0)

    def test_custom_grid_range(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result: DensityResult = density(x, from_=0, to_=6, n=100)
        assert result.x[0] == 0
        assert result.x[-1] == 6
        assert len(result.x) == 100

    def test_different_kernels(self):
        x: npt.NDArray[float] = np.random.normal(0, 1, 50)

        for kernel in ["epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"]:
            with pytest.raises(NotImplementedError, match=f"Only 'gaussian' kernel is implemented; got '{kernel}'"):
                density(x, kernel=cast(KERNEL_TYPE, kernel))

    def test_kernel_only_mode(self):
        result: npt.NDArray[float] = density([1, 2, 3], kernel="gaussian", kernel_only=True)
        expected: float = 1 / (2 * np.sqrt(np.pi))
        assert isinstance(result, float)
        assert_allclose(result, expected, rtol=1e-10)

    def test_na_removal(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
        result: DensityResult = density(x, remove_na=True)

        assert len(result.x) == 512
        assert np.all(np.isfinite(result.y))

    def test_na_without_removal(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, np.nan, 4.0, 5.0])

        with pytest.raises(ValueError, match="NA values found"):
            density(x, remove_na=False)

    def test_insufficient_data_for_nrd0(self):
        with pytest.raises(ValueError, match="at least two points"):
            density([1.0], bw="nrd0")

    def test_adjust_parameter(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result1: DensityResult = density(x, bw=1.0, adjust=1.0)
        result2: DensityResult = density(x, bw=1.0, adjust=2.0)

        assert_allclose(result2.bw, 2.0 * result1.bw, rtol=1e-10)

    def test_invalid_bandwidth_string(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])

        with pytest.raises(TypeError, match="must be a number or 'nrd0'"):
            density(x, bw="invalid")

    def test_negative_weights(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0])
        weights: npt.NDArray[float] = np.array([1.0, -1.0, 1.0])

        with pytest.raises(ValueError, match="Negative values found"):
            density(x, weights=weights)

    def test_infinite_values_in_x(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, np.inf, 4.0, 5.0])
        result: DensityResult = density(x)

        assert np.all(np.isfinite(result.y))

    def test_result_named_tuple_fields(self):
        x: npt.NDArray[float] = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result: DensityResult = density(x, n=100)

        assert hasattr(result, "x")
        assert hasattr(result, "y")
        assert hasattr(result, "x_grid")
        assert hasattr(result, "y_grid")
        assert hasattr(result, "bw")
        assert hasattr(result, "n")
