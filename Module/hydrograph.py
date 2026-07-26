"""
hydrograph.py
-------------
Defines the Hydrograph class: a time -> discharge series, linearly
interpolated between tabulated points -- exactly the behavior described
in the original manual ("le programme gere automatiquement les calculs
d'interpolation... a l'aide d'une fonction lineaire").

Used both for the incoming flood hydrograph (Qin) and for the withdrawal /
withdrawal schedule (Qwithdrawal). Time is in seconds internally; CSV input
may be given in hours (see time_unit).

Expected CSV format (Data/inflow_hydrograph.csv):
    time_h,Q_m3s
    0,30.00
    2,30.06
    ...
"""

from __future__ import annotations
import csv
import numpy as np


class Hydrograph:
    def __init__(self, csv_path: str, time_col: str = "time_h", q_col: str = "Q_m3s",
                 time_unit: str = "h"):
        t_list, q_list = [], []
        with open(csv_path, encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            for row in reader:
                t_list.append(float(row[time_col]))
                q_list.append(float(row[q_col]))

        t = np.asarray(t_list)
        if time_unit == "h":
            t = t * 3600.0
        elif time_unit == "s":
            pass
        else:
            raise ValueError(f"unknown time_unit '{time_unit}', expected 'h' or 's'")

        order = np.argsort(t)
        self.t = t[order]
        self.Q = np.asarray(q_list)[order]

    def discharge(self, t: float) -> float:
        """Q(t) [m^3/s] via linear interpolation. Holds the first/last
        tabulated value constant outside the given time range (i.e. no
        extrapolated ramp), matching the original tool's behavior of a
        defined hydrograph duration."""
        return float(np.interp(t, self.t, self.Q))

    @classmethod
    def constant(cls, value: float) -> "ConstantHydrograph":
        return ConstantHydrograph(value)


class ConstantHydrograph:
    """Convenience stand-in for a hydrograph that is simply zero or a
    fixed value (e.g. no withdrawal defined)."""
    def __init__(self, value: float = 0.0):
        self.value = value

    def discharge(self, t: float) -> float:
        return self.value
