""""
reservoir.py
------------
Defines the Reservoir class: holds the stage-surface (and optionally
stage-volume) relationship B(H), read from a CSV file.

The reservoir continuity equation dV/dt = Qin - Qout, combined with
dV = B(H)*dH, needs B(H) as a continuous function of level for the RK4
solver -- this class provides that via linear interpolation over a
tabulated curve. Linear interpolation is continuous and
monotonic-safe; if you have a smoother physical curve (e.g. from
bathymetry) you can swap the interpolation method without touching the
solver.

Expected CSV format (Data/<CaseName>/reservoir_curve.csv):
    H_masl,Surface_km2[,Volume_m3]
    90.0,0.5
    95.0,0.5
    ...

Surface area is entered in km^2 (more natural to read/edit for a real
reservoir than m^2) but is converted to m^2 immediately on load and
stored internally in m^2 -- so every other part of the tool (solver,
outlets, discharge equations, which all work in SI units of m, m^2,
m^3/s) is unaffected by this choice; only this file's CSV boundary
knows about km^2.

NOTE: replace the placeholder Surface_km2 values in Data/Template/
with your own dam's real stage-surface curve (from bathymetric survey
or design data) before using this tool for any real design or safety
decision.
"""

from __future__ import annotations
import csv
import numpy as np

KM2_TO_M2 = 1.0e6


class Reservoir:
    def __init__(self, csv_path: str):
        H_list, B_list, V_list = [], [], []
        has_volume = False

        with open(csv_path, encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            fieldnames = [fn.strip() for fn in reader.fieldnames]
            has_volume = "Volume_m3" in fieldnames
            for row in reader:
                H_list.append(float(row["H_masl"]))
                B_list.append(float(row["Surface_km2"]) * KM2_TO_M2)
                if has_volume:
                    V_list.append(float(row["Volume_m3"]))

        order = np.argsort(H_list)
        self.H = np.asarray(H_list)[order]
        self.B = np.asarray(B_list)[order]  # stored internally in m^2
        self.V = np.asarray(V_list)[order] if has_volume else None

        if np.any(np.diff(self.H) <= 0):
            raise ValueError("reservoir_curve.csv must have strictly increasing H_masl values")

    def surface(self, H: float) -> float:
        """Reservoir surface area B(H) [m^2] (converted from the CSV's
        km^2 at load time), linearly interpolated. Extrapolates (flat)
        beyond the table bounds rather than raising, so a transient RK4
        sub-step that slightly overshoots the table doesn't crash the
        simulation."""
        return float(np.interp(H, self.H, self.B))

    def volume(self, H: float) -> float | None:
        """Reservoir volume V(H) [m^3], if a Volume_m3 column was supplied."""
        if self.V is None:
            return None
        return float(np.interp(H, self.H, self.V))

    @property
    def H_min(self) -> float:
        return float(self.H[0])

    @property
    def H_max(self) -> float:
        return float(self.H[-1])