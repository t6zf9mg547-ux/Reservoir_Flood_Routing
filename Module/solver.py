"""
solver.py
---------
Implements the 4th-order Runge-Kutta integration of the reservoir
continuity equation:

    dH/dt = f(H, t) = (Qin(t) - Qwithdrawal(t) - sum(Qs_i(H, t))) / B(H)

This file should not normally need editing -- reservoir shape, inflow,
withdrawal, and outlet behavior are all supplied as objects with a
common interface (see reservoir.py, hydrograph.py, outlets.py).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List
import numpy as np

from reservoir import Reservoir
from hydrograph import Hydrograph
from outlets import Outlet


@dataclass
class SimulationResult:
    t_s: np.ndarray          # time [s]
    t_h: np.ndarray          # time [h]
    H: np.ndarray            # reservoir level [m a.s.l.]
    Qin: np.ndarray          # inflow [m3/s]
    Qwithdrawal: np.ndarray     # withdrawal [m3/s]
    Qout_by_outlet: dict     # {outlet_name: array of Q [m3/s]}
    Qout_total: np.ndarray   # sum of all outlet discharges [m3/s]
    V_rate_cmh: np.ndarray   # dH/dt expressed in cm/h (a common convention
                             # for reporting rate-of-rise in reservoir operations)


def f(H: float, t: float, reservoir: Reservoir, inflow: Hydrograph,
      withdrawal: Hydrograph, outlets: List[Outlet]) -> float:
    """dH/dt at a given (H, t), used for RK4 trial evaluations.
    Does NOT commit any outlet state (see outlets.py docstring)."""
    Qin = inflow.discharge(t)
    Qwithdrawal = withdrawal.discharge(t)
    Qout = sum(o.discharge(H, t) for o in outlets)
    B = reservoir.surface(H)
    return (Qin - Qwithdrawal - Qout) / B


def run_simulation(reservoir: Reservoir, inflow: Hydrograph, withdrawal: Hydrograph,
                    outlets: List[Outlet], H0: float, t_max: float, dt: float,
                    print_every: int = 1) -> SimulationResult:
    """
    Integrate the reservoir level from t=0 to t=t_max using classic RK4
    with fixed step dt (all in seconds), then record output every
    `print_every` steps -- lets you take small, accurate integration
    steps while writing a coarser, more readable output series.

    Parameters (all read from Data/<CaseName>/scalars.csv by run_case.py):
        H0          : initial reservoir level [m a.s.l.]
        t_max       : total simulation duration [s]
        dt          : integration timestep [s]
        print_every : record output every N accepted steps
    """
    for o in outlets:
        o.reset()

    n_steps = int(round(t_max / dt))
    t = 0.0
    H = H0

    rows = []  # (t, H) recorded at print frequency

    def record(t_, H_):
        Qin = inflow.discharge(t_)
        Qwithdrawal = withdrawal.discharge(t_)
        Qs = {o.name: o.discharge(H_, t_) for o in outlets}
        Qout_tot = sum(Qs.values())
        dHdt = (Qin - Qwithdrawal - Qout_tot) / reservoir.surface(H_)
        rows.append((t_, H_, Qin, Qwithdrawal, Qs, Qout_tot, dHdt * 100 * 3600))  # cm/h

    record(t, H)

    for step in range(1, n_steps + 1):
        k1 = f(H, t, reservoir, inflow, withdrawal, outlets)
        k2 = f(H + 0.5 * dt * k1, t + 0.5 * dt, reservoir, inflow, withdrawal, outlets)
        k3 = f(H + 0.5 * dt * k2, t + 0.5 * dt, reservoir, inflow, withdrawal, outlets)
        k4 = f(H + dt * k3, t + dt, reservoir, inflow, withdrawal, outlets)

        H = H + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + dt

        # commit outlet state exactly once, on the accepted level
        for o in outlets:
            o.commit(H)

        if step % print_every == 0 or step == n_steps:
            record(t, H)

    t_s = np.array([r[0] for r in rows])
    H_arr = np.array([r[1] for r in rows])
    Qin_arr = np.array([r[2] for r in rows])
    Qwithdrawal_arr = np.array([r[3] for r in rows])
    outlet_names = outlets[0].name if False else [o.name for o in outlets]
    Qout_by_outlet = {name: np.array([r[4][name] for r in rows]) for name in outlet_names}
    Qout_total = np.array([r[5] for r in rows])
    Vrate = np.array([r[6] for r in rows])

    return SimulationResult(
        t_s=t_s, t_h=t_s / 3600.0, H=H_arr, Qin=Qin_arr, Qwithdrawal=Qwithdrawal_arr,
        Qout_by_outlet=Qout_by_outlet, Qout_total=Qout_total, V_rate_cmh=Vrate,
    )
