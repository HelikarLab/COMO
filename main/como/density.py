from __future__ import annotations

import sys
from collections.abc import Sequence
from typing import Literal, NamedTuple, cast

import numpy as np
import numpy.typing as npt

from como.approx import Approx, approx


class DensityResult(NamedTuple):
    x: npt.NDArray[float]
    y: npt.NDArray[float]
    x_grid: npt.NDArray[float]
    y_grid: npt.NDArray[float]
    bw: float
    n: int


NUMBER = int | float | np.number


def bin_distance(x: npt.NDArray[float], weights: npt.NDArray[float], lo: NUMBER, up: NUMBER, n: int) -> npt.NDArray[float]:
    """Bin weighted distances."""
    ixmin: int = 0
    ixmax: int = n - 2
    delta: float = (up - lo) / (n - 1)

    y: npt.NDArray[float] = np.zeros((2 * n,), dtype=float)
    for i in range(x.size):
        i: int
        xpos: float = (x[i] - lo) / delta
        if xpos > sys.maxsize or xpos < -sys.maxsize:  # avoid integer overflows (taken from R's massdist.c)
            continue
        ix: int = np.floor(xpos).astype(int)
        fx: float = xpos - ix
        wi: float = weights[i]
        if ixmin <= ix <= ixmax:
            y[ix] += (1 - fx) * wi
            y[ix + 1] += fx * wi
        elif ix == -1:
            y[0] += fx * wi
        elif ix == ixmax + 1:
            y[ix] += (1 - fx) * wi
    return y


def dnorm(x: float, mean: NUMBER = 0.0, sd: NUMBER = 1.0, log: bool = False, fast_dnorm: bool = False) -> float:
    """Density function for the normal distribution.

    This is a reproduction of R's `density` function.
    Neither SciPy nor NumPy are capable of producing KDE values that align with R, and as a result,
        a manual translation of R's KDE implementation was necessary.

    References:
        1) Documentation: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/Normal
        2) Source code (2025-OCT-30): https://github.com/wch/r-source/blob/3f7e2528990351bc4b0d1f1b237714668ab4c0c5/src/nmath/dnorm.c

    Args:
        x: Value at which to evaluate the density.
        mean: Mean of the normal distribution.
        sd: Standard deviation of the normal distribution.
        log: If True, return the log density.
        fast_dnorm: If True, use a faster but less accurate calculation for small `x`.


    """
    # Constants
    m_ln2: float = 0.693147180559945309417232121458  # ln(2)
    m_1_sqrt_2pi: float = 0.398942280401432677939946059934  # 1/sqrt(2pi)
    m_ln_sqrt_2pi: float = 0.918938533204672741780329736406  # log(sqrt(2*pi))
    r_d__0 = float(-np.inf) if log else 0.0  # R's R_D__0: (log_p ? ML_NEGINF : 0.)

    # 11 bits exponent, where one bit is used to indicate special numbers (e.g. NaN and Infinity),
    #   so the maximum exponent is 10 bits wide (2^10 == 1024).
    # dbl_min_exp = -1022
    dbl_min_exp = -1074

    # The mantissa is 52 bits wide, but because numbers are normalized the initial binary 1 is represented
    #   implicitly (the so-called "hidden bit"), which leaves us with the ability to represent 53 bits wide mantissa.
    dbl_mant_dig = 53

    mean: float = float(mean)
    sd: float = float(sd)

    if np.isnan(x) or np.isnan(mean) or np.isnan(sd):
        return x + mean + sd
    if sd < 0.0:
        return float(np.nan)
    if not np.isfinite(sd):
        return r_d__0
    if not np.isfinite(x) and x == mean:
        return float(np.nan)
    if sd == 0.0:
        return float(np.inf) if x == mean else r_d__0

    # From this point on, dividing by `sd` is safe because we know it is not 0
    z = (x - mean) / sd
    if (not np.isfinite(z)) and (x == mean):
        return float(np.nan)

    if not np.isfinite(z):
        return r_d__0

    a = np.fabs(z)
    if a >= 2 * np.sqrt(np.finfo(float).max):
        return r_d__0
    if log:
        return -(m_ln_sqrt_2pi + 0.5 * a * a + np.log(sd))
    if a < 5 or fast_dnorm:  # for `x < 5`, this is more accurate but less fast
        return m_1_sqrt_2pi * np.exp(-0.5 * a * a) / sd

    # underflow boundary
    boundary = np.sqrt(-2.0 * m_ln2 * (dbl_min_exp + 1 - dbl_mant_dig))
    if a > boundary:
        return 0.0

    # Now, to get full accuracy, split x into two parts,
    #   x = x1+x2, such that |x2| <= 2^-16.
    #   Assuming that we are using IEEE doubles, that means that
    #   x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).
    #   If we do not have IEEE this is still an improvement over the naive formula.
    a1 = np.ldexp(np.rint(np.ldexp(a, 16)), -16)
    a2 = a - a1
    return (m_1_sqrt_2pi / sd) * (np.exp(-0.5 * a1 * a1) * np.exp((-a1 * a2) - (0.5 * a2 * a2)))


def nrd0(x: npt.NDArray[float]) -> float:
    """Calculate nrd0 from R source.

    This bandwidth calculation matches R's

    References:
        1) Documentation: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/bandwidth
        2) Source code (as of 2025-OCT-30, copied below): https://github.com/wch/r-source/blob/trunk/src/library/stats/R/bandwidths.R
            ```R
            bw.nrd0 <- function (x)
            {
                if(length(x) < 2L) stop("need at least 2 data points")
                hi <- sd(x)
                if(!(lo <- min(hi, IQR(x)/1.34)))# qnorm(.75) - qnorm(.25) = 1.34898
                    (lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1.)
                0.9 * lo * length(x)^(-0.2)
            }
            ```
    """
    if x.size < 2:
        raise ValueError("need at least 2 data points")

    hi = np.std(x, ddof=1)
    q75, q25 = cast(tuple[float, float], cast(object, np.percentile(x, [75, 25])))
    iqr_over_134 = (q75 - q25) / 1.34

    # We are using a cascading series of checks on `lo` to  make sure it is a non-zero value
    lo = min(hi, iqr_over_134)
    if lo == 0:
        if hi != 0:
            lo = hi
        elif abs(x[0]) != 0:
            lo = abs(x[0])
        else:
            lo = 1.0

    return float(0.9 * lo * x.size ** (-1 / 5))


def density(
    x: Sequence[NUMBER] | npt.NDArray[NUMBER],
    bw: NUMBER | Literal["nrd0"] = "nrd0",
    adjust: NUMBER = 1,
    kernel: Literal["gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"] = "gaussian",
    weights: Sequence[NUMBER] | None = None,
    n: int = 512,
    from_: NUMBER | None = None,
    to_: NUMBER | None = None,
    cut: int = 3,
    ext: int = 4,
    remove_na: bool = False,
    kernel_only: bool = False,
) -> DensityResult:
    """Compute kernel density estimates (KDE) using FFT method.

    This is a reproduction of R's `density` function.
    Neither SciPy nor NumPy are capable of producing KDE values that align with R, and as a result,
        a manual translation of R's KDE implementation was necessary.

    Args:
        x: Input data points.
        bw: Bandwidth for the kernel. If "nrd0", uses R's nrd0 method.
        adjust: Adjustment factor for the bandwidth.
        kernel: Kernel type to use.
        weights: Optional weights for each data point.
        n: Number of points in the output grid.
        from_: Start of the grid (calculated automatically if not provided).
        to_: End of the grid (calculated automatically if not provided).
        cut: Number of bandwidths to extend the grid on each side.
        ext: Number of bandwidths to extend the grid for FFT calculation.
        remove_na: Whether to remove NA values from `x`.
        kernel_only: If True, returns only the integral of the kernel function.

    Returns:
        DensityResult: A named tuple containing:
            x: The x-coordinates of the density estimate.
            y: The estimated density values at the x-coordinates.
            x_grid: The grid of x-coordinates used for FFT calculation.
            y_grid: The density values on the FFT grid.
            bw: The bandwidth used.
            n: The number of points in the output grid.
    """
    if not isinstance(x, Sequence | np.ndarray) or (len(x) < 2 and bw == "nrd0"):
        raise ValueError("Need at at least two points to select a bandwidth automatically using 'nrd0'")
    if kernel_only:
        if kernel == "gaussian":
            return 1 / (2 * np.sqrt(np.pi))
        elif kernel == "epanechnikov":
            return 3 / (5 * np.sqrt(5))
        elif kernel == "rectangular":
            return np.sqrt(3) / 6
        elif kernel == "triangular":
            return np.sqrt(6) / 9
        elif kernel == "biweight":
            return 5 * np.sqrt(7) / 49
        elif kernel == "cosine":
            return 3 / 4 * np.sqrt(1 / 3 - 2 / np.pi**2)
        elif kernel == "optcosine":
            return np.sqrt(1 - 8 / np.pi**2) * np.pi**2 / 16

    if kernel != "gaussian":
        raise NotImplementedError(f"Only 'gaussian' kernel is implemented; got '{kernel}'")

    x: npt.NDArray[float] = np.asarray(x, dtype=float)

    has_weights = weights is not None
    weights: npt.NDArray[float] | None = np.asarray(weights, float) if weights is not None else None
    if has_weights and (weights is not None and weights.size != x.size):
        raise ValueError(f"The length of provided weights does not match the length of x: {weights.size} != {x.size}")

    x_na: npt.NDArray[np.bool_] = np.isnan(x)
    if np.any(x_na):
        if remove_na:
            x: npt.NDArray[float] = x[~x_na]
            if has_weights and weights is not None:
                true_d = weights.sum().astype(float) == 1
                weights = weights[~x_na]
                if true_d:
                    weights = weights / weights.sum()
        else:
            raise ValueError("NA values found in 'x'. Set 'remove_na=True' to ignore them.")

    nx = n
    x_finite = np.isfinite(x)
    if np.any(~x_finite):
        x = x[x_finite]
        nx = x.size
    if not has_weights:
        weights: npt.NDArray[float] = np.full(shape=nx, fill_value=1 / nx, dtype=float)
        total_mass = nx / n
    else:
        weights = np.asarray(weights)
        if not np.all(np.isfinite(weights)):
            raise ValueError("Non-finite values found in 'weights'.")
        if np.any(weights < 0):
            raise ValueError("Negative values found in 'weights'.")
        wsum: float = weights.sum()
        if np.any(~x_finite):
            weights = weights[x_finite]
            total_mass = float(weights.sum() / wsum)
        else:
            total_mass = float(1)

    n_user: int = n
    n = max(n, 512)
    if n > 512:  # round n up to the next power of 2 (i.e., 2^8=512, 2^9=1024)
        n: int = int(2 ** np.ceil(np.log2(n)))

    if isinstance(bw, str) and bw != "nrd0":
        raise TypeError("Bandwidth 'bw' must be a number or 'nrd0'.")
    elif isinstance(bw, str) and bw == "nrd0":
        bw_calc = nrd0(x)
    else:
        bw_calc = float(bw)
    if not np.isfinite(bw_calc):
        raise ValueError("Calculated bandwidth 'bw' is not finite.")
    bw_calc *= adjust

    if bw_calc <= 0:
        raise ValueError("Bandwidth 'bw' must be positive.")

    # have to use `... if ... else` because `0` is falsey, resulting in the right-half being used instead of the user-provided value
    from_ = float(from_ if from_ is not None else min(x) - cut * bw_calc)
    to_ = float(to_ if to_ is not None else max(x) + cut * bw_calc)

    if not np.isfinite(from_):
        raise ValueError("'from_' is not finite.")
    if not np.isfinite(to_):
        raise ValueError("'to_' is not finite.")

    lo = float(from_ - ext * bw_calc)
    up = float(to_ + ext * bw_calc)

    y: npt.NDArray[float] = bin_distance(x, weights, lo, up, n) * total_mass
    kords: npt.NDArray[float] = np.linspace(start=0, stop=((2 * n - 1) / (n - 1) * (up - lo)), num=2 * n, dtype=float)
    kords[n + 1 : 2 * n] = -kords[n:1:-1]  # mirror/negate tail: R's kords[n:2] will index from the reverse if `n`>2

    # Initial diverge here (inside dnorm calculation)
    kords: npt.NDArray[float] = np.asarray([dnorm(i, sd=bw_calc) for i in kords], dtype=float)

    fft_y: npt.NDArray[np.complex128] = np.fft.fft(y)
    conj_fft_kords: npt.NDArray[np.complex128] = np.conjugate(np.fft.fft(kords))
    # Must multiply by `kords.size` because R does not produce a normalize inverse FFT, but NumPy normalizes by `1/size`
    kords: npt.NDArray[np.complex128] = np.fft.ifft(fft_y * conj_fft_kords) * kords.size
    kords: npt.NDArray[float] = (np.maximum(0, kords.real)[0:n]) / y.size  # for values of kords, get 0 or kords[i], whichever is larger
    xords: npt.NDArray[float] = np.linspace(lo, up, num=n, dtype=float)

    # xp=known x-coords, fp=known y-cords, x=unknown x-coords; returns interpolated (e.g., unknown) y-coords
    interp_x: npt.NDArray[float] = np.linspace(from_, to_, num=n_user)
    interp_y: Approx = approx(xords, kords, interp_x)

    return DensityResult(
        x=interp_x,
        y=interp_y.y,
        x_grid=xords,
        y_grid=kords,
        bw=bw_calc,
        n=n,
    )
