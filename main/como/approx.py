from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Literal, NamedTuple

import numpy as np
import numpy.typing as npt


class RegularizedArray(NamedTuple):
    x: npt.NDArray[np.number]
    y: npt.NDArray[np.number]
    not_na: npt.NDArray[bool]
    kept_na: bool


class Approx(NamedTuple):
    x: npt.NDArray[float]
    y: npt.NDArray[float]


def _coerce_to_float_array(a):
    """Helper to ensure input is a 1D float array."""
    arr = np.asarray(a, dtype=float)
    if arr.ndim != 1:
        arr = arr.ravel()
    return arr


def _regularize_values(x: np.ndarray, y: np.ndarray, ties, na_rm: bool) -> dict[str, Any]:
    """Minimal reimplementation of R's regularize.values() used by approx().

      - Removes NA pairs if na_rm is True (else requires x to have no NA).
      - Sorts by x (stable).
      - Collapses duplicate x via `ties` aggregator.
    Returns dict with:
      x: sorted unique x
      y: aggregated y aligned with x
      not_na: boolean mask of y that are not NaN
      kept_na: True iff any NaN remained in y after regularization
    """
    if na_rm:
        # Remove pairs where x or y is NA
        keep = ~np.isnan(x) & ~np.isnan(y)
        x = x[keep]
        y = y[keep]
        kept_na = False
    else:
        # R's C code errors if na_rm=False and any x is NA
        if np.isnan(x).any():
            raise ValueError("approx(x,y, .., na.rm=False): NA values in x are not allowed")
        kept_na = np.isnan(y).any()

    if x.size == 0:
        return RegularizedArray(x=x, y=y, not_na=np.array([], dtype=bool), kept_na=kept_na)

    # Use a stable sort (mergesort) to match R's order()
    order = np.argsort(x, kind="mergesort")
    x_sorted = x[order]
    y_sorted = y[order]

    # Handle the 'ties' argument
    if callable(ties):
        agg_fn = ties
    else:
        ties_str = "mean" if ties is None else str(ties).lower()
        if ties_str in ("mean", "avg", "average"):
            agg_fn = np.mean
        elif ties_str in ("first", "left"):

            def agg_fn(a):
                return a[0]
        elif ties_str in ("last", "right"):

            def agg_fn(a):
                return a[-1]
        elif ties_str == "min":
            agg_fn = np.min
        elif ties_str == "max":
            agg_fn = np.max
        elif ties_str == "median":
            agg_fn = np.median
        elif ties_str == "sum":
            agg_fn = np.sum
        else:
            raise ValueError("Unsupported `ties`; use a callable or one of 'mean', 'first', 'last', 'min', 'max', 'median', 'sum'.")

    # Find unique x values and their indices/counts
    unique_x, start_idx, counts = np.unique(x_sorted, return_index=True, return_counts=True)

    # Aggregate y values for tied x values
    y_agg = np.empty(unique_x.shape, dtype=float)
    for k, (start, count) in enumerate(zip(start_idx, counts, strict=True)):
        seg = y_sorted[start : start + count]
        # Note: aggregators like np.mean will return NaN if `seg` contains NaN,
        # matching R's default (mean(..., na.rm = FALSE)).
        y_agg[k] = agg_fn(seg)

    not_na = ~np.isnan(y_agg)
    # Check if any NaNs remain in the *aggregated* y values
    kept_na_final = np.any(~not_na) if np.any(np.isnan(y_agg)) else False
    return {"x": unique_x, "y": y_agg, "not_na": not_na, "kept_na": kept_na_final}
    return RegularizedArray(x=unique_x, y=y_agg, not_na=not_na, kept_na=kept_na_final)


def approx(
    x: Sequence[float],
    y: Sequence[float] | None = None,
    xout: Sequence[float] | None = None,
    method: str | int = "linear",
    n: int = 50,
    yleft: float | None = None,
    yright: float | None = None,
    rule: int | Sequence[int] = 1,
    f: float = 0.0,
    ties: str | Callable[[np.ndarray], float] = "mean",
    na_rm: bool = True,
) -> Approx:
    """Faithful Python port of R's `approx` function.

    This implementation aims to replicate the behavior of R's `approx`
    function, including its handling of `ties`, `rule`, `na_rm`, and
    `method`.

    Args:
        x: Coordinates, or y-values if `y` is None.
        y: y-coordinates. If None, `x` is treated as `y` and `x` becomes `1..n`.
        xout: Points at which to interpolate.
        method: Interpolation method. "linear" (1) or "constant" (2).
        n: If `xout` is not provided, interpolation happens at `n`
             equally spaced points spanning the range of `x`.
        yleft: Value to use for extrapolation to the left.
               Defaults to `NA` (np.nan) if `rule` is 1.
        yright: Value to use for extrapolation to the right.
                Defaults to `NA` (np.nan) if `rule` is 1.
        rule: Extrapolation rule.
              - 1: Return `np.nan` for points outside the `x` range.
              - 2: Use `yleft` and `yright` for points outside the range.
        f: For `method="constant"`, determines the split point.
           `f=0` is left-step, `f=1` is right-step, `f=0.5` is midpoint.
        ties: How to handle duplicate `x` values.
              Can be 'mean', 'first', 'last', 'min', 'max', 'median', 'sum',
              or a callable function.
        na_rm: If True, `NA` pairs are removed before interpolation.
               If False, `NA`s in `x` cause an error, `NA`s in `y`
               are propagated.

    Returns:
        `Approx` class with:
          - 'x': numpy array of xout used
          - 'y': numpy array of interpolated values
    """
    # One-argument form: approx(y) -> x := 1..n, y := y
    if y is None:
        y_arr = _coerce_to_float_array(x)
        x_arr = np.arange(1.0, y_arr.size + 1.0, dtype=float)
    else:
        x_arr = _coerce_to_float_array(x)
        y_arr = _coerce_to_float_array(y)
        if x_arr.size != y_arr.size:
            raise ValueError("x and y must have same length")

    # --- Method normalization ---
    # (matches R's pmatch() result: 1=linear, 2=constant)
    if isinstance(method, str):
        m = method.strip().lower()
        if m.startswith("lin"):
            method_code = 1
        elif m.startswith("const"):
            method_code = 2
        else:
            raise ValueError("invalid interpolation method")
    else:
        if method in (1, 2):
            method_code = int(method)
        else:
            raise ValueError("invalid interpolation method")

    # --- Rule normalization ---
    if isinstance(rule, list | tuple | np.ndarray):
        rlist = list(rule)
        if not (1 <= len(rlist) <= 2):
            raise ValueError("`rule` must have length 1 or 2")
        if len(rlist) == 1:
            rlist = [rlist[0], rlist[0]]
    else:
        rlist = [rule, rule]

    # --- Regularize values ---
    # Sort by x, collapse ties, and handle NAs like R's regularize.values()
    r = _regularize_values(x_arr, y_arr, ties, na_rm)
    x_reg = r["x"]
    y_reg = r["y"]
    not_na_mask = r["not_na"]
    r: RegularizedArray = _regularize_values(x_arr, y_arr, na_rm=na_rm, ties=ties)
    # no_na is True if we don't have to worry about NAs in y_reg
    no_na = na_rm or (not r["kept_na"])
    # nx is the number of *valid* (non-NA) points for interpolation
    nx = x_reg.size if no_na else int(np.count_nonzero(not_na_mask))

    if np.isnan(nx):
        raise ValueError("invalid length(x)")

    # R's C_ApproxTest logic
    if nx <= 1:
        if method_code == 1:
            raise ValueError("need at least two non-NA values to interpolate")
        if nx == 0:
            raise ValueError("zero non-NA points")

    # --- Set extrapolation values (yleft/yright) ---
    # This logic matches R's C code (R_approxtest)
    if no_na:
        y_first = y_reg[0]
        y_last = y_reg[-1]
    else:
        # Find first and last non-NA y values
        y_valid = y_reg[not_na_mask]
        y_first = y_valid[0]
        y_last = y_valid[-1]

    yleft_val = (np.nan if int(rlist[0]) == 1 else y_first) if yleft is None else float(yleft)
    yright_val = (np.nan if int(rlist[1]) == 1 else y_last) if yright is None else float(yright)

    # --- Fabricate xout if missing ---
    if xout is None:
        if n <= 0:
            raise ValueError("'approx' requires n >= 1")
        if no_na:
            x_first = x_reg[0]
            x_last = x_reg[nx - 1]
        else:
            x_valid = x_reg[not_na_mask]
            x_first = x_valid[0]
            x_last = x_valid[-1]
        xout_arr = np.linspace(x_first, x_last, num=int(n), dtype=float)
    else:
        xout_arr = _coerce_to_float_array(xout)

    # --- C_ApproxTest checks ---
    if method_code == 2 and (not np.isfinite(f) or f < 0.0 or f > 1.0):
        raise ValueError("approx(): invalid f value")
    if not no_na:
        # R re-filters x and y here if NAs remained
        x_reg = x_reg[not_na_mask]
        y_reg = y_reg[not_na_mask]

    # --- Interpolation ---
    # This section vectorized the logic from R's approx1 and R_approxfun
    yout = np.empty_like(xout_arr, dtype=float)
    mask_nan_xout = np.isnan(xout_arr)
    yout[mask_nan_xout] = np.nan
    mask_valid = ~mask_nan_xout

    if x_reg.size == 0:
        # No valid points to interpolate against
        yout[mask_valid] = np.nan
        return {"x": xout_arr, "y": yout}

    xv = xout_arr[mask_valid]
    left_mask = xv < x_reg[0]
    right_mask = xv > x_reg[-1]
    mid_mask = ~(left_mask | right_mask)

    res = np.empty_like(xv)
    res[left_mask] = yleft_val
    res[right_mask] = yright_val

    if np.any(mid_mask):
        xm = xv[mid_mask]

        # Find indices
        # j_right[k] = index of first x_reg > xm[k]
        j_right = np.searchsorted(x_reg, xm, side="right")
        # j_left[k] = index of first x_reg >= xm[k]
        j_left = np.searchsorted(x_reg, xm, side="left")

        # Points that exactly match an x_reg value
        eq_mask = j_left != j_right
        # Points that fall between x_reg values
        in_interval_mask = ~eq_mask

        res_mid = np.empty_like(xm)

        if np.any(eq_mask):
            # For exact matches, use the corresponding y_reg value
            # R's C code uses the 'j_left' index here
            res_mid[eq_mask] = y_reg[j_left[eq_mask]]

        if np.any(in_interval_mask):
            # R's C code sets i = j-1, where j is the 'right' index
            j = j_right[in_interval_mask]
            i = j - 1

            # Get the surrounding x and y values
            xi = x_reg[i]
            xj = x_reg[j]
            yi = y_reg[i]
            yj = y_reg[j]

            if method_code == 1:  # linear
                t = (xm[in_interval_mask] - xi) / (xj - xi)
                res_mid[in_interval_mask] = yi + (yj - yi) * t
            else:  # constant
                # This matches R_approxfun's constant logic
                if f == 0.0:
                    res_mid[in_interval_mask] = yi
                elif f == 1.0:
                    res_mid[in_interval_mask] = yj
                else:
                    # R computes (1-f)*yi + f*yj, but carefully
                    f1 = 1.0 - f
                    f2 = f
                    part = np.zeros_like(yi)
                    if f1 != 0.0:
                        part = part + yi * f1
                    if f2 != 0.0:
                        part = part + yj * f2
                    res_mid[in_interval_mask] = part

        res[mid_mask] = res_mid

    yout[mask_valid] = res
    return Approx(x=xout_arr, y=yout)
