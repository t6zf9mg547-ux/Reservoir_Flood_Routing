"""
condition_index.py
-------------------
Lives in Module/reliability/. Converts a Condition Index (CI) estimate
into a failure probability.

WHERE THE NUMBERS COME FROM -- this module does not invent, look up,
or infer CI values from anywhere. CI_mean and CI_sd are always supplied
by the caller (read from a case's own components.csv/
common_cause_events.csv), reflecting whatever evidence -- an inspection
report, a walk-down, or plain engineering judgement when no better
evidence exists -- the person who wrote that CSV had at the time. This
module's only job is the probability arithmetic once those two numbers
are given; see Data/Template/reliability/*.csv for the input format.

CI scale (0-100, per the reference methodology this engine implements):
    85-100  Excellent       10-24  Very poor
    70-84   Good             0-9   Failed
    55-69   Fair
    40-54   Marginal
    25-39   Poor

CI is treated as NORMALLY DISTRIBUTED around CI_mean with spread
CI_sd -- not as a single fixed truth -- because inspection judgement
and limited evidence are themselves uncertain, and collapsing that
uncertainty into one point estimate would overstate how well the
condition is actually known. A tight, well-evidenced estimate (a clear
photo, a direct measurement) should get a SMALL CI_sd; a rough
judgement call with no direct evidence should get a WIDE one -- that
distinction is the caller's to make, not this module's.

Given a failure threshold CI_failure (the CI value at or below which
the component is considered non-functional), the failure probability
is the fraction of that normal distribution at or below the
threshold:

    Pf = P(CI <= CI_failure) = Phi((CI_failure - CI_mean) / CI_sd)

using Phi, the standard normal CDF (implemented here via math.erf --
no scipy dependency needed for one function).
"""

from __future__ import annotations
import math


def pf_from_ci(ci_mean: float, ci_sd: float, ci_failure: float) -> float:
    """Failure probability implied by a CI estimate.

    ci_mean, ci_sd : the CI estimate, on the 0-100 scale (see module
        docstring) -- supplied by the caller, never computed here.
    ci_failure : the CI value at or below which this component is
        considered failed/non-functional. Typically the same threshold
        across a whole case's components.csv (a single documented
        assumption), but nothing here requires that -- each row can
        carry its own.

    Returns a probability in [0, 1]. If ci_sd <= 0 (a stated CI with no
    uncertainty at all -- unusual, but not disallowed), falls back to a
    step function: Pf = 1.0 if ci_mean <= ci_failure else 0.0.
    """
    if ci_sd <= 0:
        return 1.0 if ci_mean <= ci_failure else 0.0
    z = (ci_failure - ci_mean) / ci_sd
    # Phi(z) via the standard erf identity -- avoids a scipy dependency
    # for a single CDF evaluation.
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
