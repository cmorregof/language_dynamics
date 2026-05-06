"""
Language Dynamical Equation (LDE) simulations for Old English -> Middle English.

v3 changes from v2:
- No-contact counterfactual now runs strictly n=4 throughout all three phases.
  When old_norse_fraction == 0.0, Phase 2 uses GRAMMARS_4 and the G5 column
  is padded with NaN for output consistency. This prevents G5 from appearing
  via the mutation term (a*u > 0) in a scenario where it should not exist.
- Added plot_no_contact_vs_moderate_comparison(): two-panel figure for §6.3.
- Rewrote __main__: uses Path(__file__) for relative output paths; generates
  all seven manuscript figures in a single run.

Unchanged from v2:
- RHS make_symmetric_lde_rhs: no clipping or renormalization inside the vector field.
- Mass conservation and positivity diagnostics.
- All scenario parameter values match Table 2 of the manuscript.
- rk4_convergence_check for appendix verification.
- All other plotting functions.

Author: Carlos Manuel Orrego Franco
Code cleaned/restructured with AI assistance, May 2026.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

Array = np.ndarray


GRAMMARS_4 = ["West Saxon (G1)", "Mercian (G2)", "Northumbrian (G3)", "Kentish (G4)"]
GRAMMARS_5 = GRAMMARS_4 + ["Old Norse (G5)"]


@dataclass(frozen=True)
class PhaseConfig:
    """Configuration for a single historical/modeling phase."""

    name: str
    q: float
    a: float
    t_start: float
    t_end: float
    n_eval: int
    labels: list[str]


@dataclass(frozen=True)
class SimulationDiagnostics:
    """Numerical diagnostics for one integration segment."""

    max_mass_error: float
    min_component: float
    success: bool
    message: str


def simplex_entropy(x: Array, normalized: bool = False) -> float:
    """Shannon entropy H(x) = -sum_j x_j log(x_j)."""
    x = np.asarray(x, dtype=float)
    positive = x[x > 0]
    h = -float(np.sum(positive * np.log(positive)))
    if normalized and len(x) > 1:
        h /= np.log(len(x))
    return h


def concentration_m2(x: Array) -> float:
    """Concentration/coherence proxy M2 = sum_j x_j^2."""
    x = np.asarray(x, dtype=float)
    return float(np.dot(x, x))


def make_symmetric_lde_rhs(q: float, a: float, n: int) -> Callable[[float, Array], Array]:
    """
    Build RHS for the fully symmetric language dynamical equation.

    dx_j/dt = (1-a)[(q-u)x_j^2 + u M2 - x_j M2] - a u n x_j + a u,
    where u = (1-q)/(n-1), M2 = sum_i x_i^2.

    No clipping or renormalization is performed here; the simplex is
    invariant under this flow by the positivity argument in Mitchener (2003).
    """
    if n < 2:
        raise ValueError("n must be at least 2.")
    if not (0 <= a <= 1):
        raise ValueError("a must lie in [0, 1].")
    if not (0 <= q <= 1):
        raise ValueError("q must lie in [0, 1].")

    u = (1.0 - q) / (n - 1)

    def rhs(t: float, x: Array) -> Array:
        x = np.asarray(x, dtype=float)
        m2 = np.dot(x, x)
        return (1.0 - a) * ((q - u) * x**2 + u * m2 - x * m2) - a * u * n * x + a * u

    return rhs


def validate_simplex(x0: Array, tol: float = 1e-10) -> Array:
    """Validate and normalize an initial condition on the simplex."""
    x0 = np.asarray(x0, dtype=float)
    if np.any(x0 < -tol):
        raise ValueError("Initial condition has negative components.")
    total = float(np.sum(x0))
    if abs(total - 1.0) > tol:
        x0 = x0 / total
    return x0


def simulate_segment(
    x0: Array,
    phase: PhaseConfig,
    method: str = "RK45",
    rtol: float = 1e-10,
    atol: float = 1e-12,
) -> tuple[Array, Array, SimulationDiagnostics]:
    """Simulate one phase using solve_ivp."""
    x0 = validate_simplex(x0)
    n = len(x0)
    if len(phase.labels) != n:
        raise ValueError(
            f"Phase '{phase.name}' has {len(phase.labels)} labels but x0 has length {n}."
        )

    rhs = make_symmetric_lde_rhs(q=phase.q, a=phase.a, n=n)
    t_eval = np.linspace(phase.t_start, phase.t_end, phase.n_eval)

    sol = solve_ivp(
        rhs,
        (phase.t_start, phase.t_end),
        x0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )

    y = sol.y.T
    diagnostics = SimulationDiagnostics(
        max_mass_error=float(np.max(np.abs(np.sum(y, axis=1) - 1.0))),
        min_component=float(np.min(y)),
        success=bool(sol.success),
        message=str(sol.message),
    )

    if not sol.success:
        raise RuntimeError(f"Integration failed in phase '{phase.name}': {sol.message}")

    return sol.t, y, diagnostics


def rk4_fixed_step(
    rhs: Callable[[float, Array], Array],
    x0: Array,
    t_start: float,
    t_end: float,
    h: float,
) -> tuple[Array, Array]:
    """
    Fixed-step classical RK4 integrator for appendix/convergence checks.

    Not used as the production integrator; see simulate_segment for RK45.
    """
    x0 = validate_simplex(x0)
    n_steps = int(np.ceil((t_end - t_start) / h))
    t = np.linspace(t_start, t_end, n_steps + 1)
    y = np.zeros((n_steps + 1, len(x0)))
    y[0] = x0

    for k in range(n_steps):
        dt = t[k + 1] - t[k]
        x = y[k]
        k1 = rhs(t[k], x)
        k2 = rhs(t[k] + dt / 2, x + dt * k1 / 2)
        k3 = rhs(t[k] + dt / 2, x + dt * k2 / 2)
        k4 = rhs(t[k] + dt, x + dt * k3)
        y[k + 1] = x + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return t, y


def build_three_phase_scenario(
    a: float = 0.9,
    x0_oe: Iterable[float] = (0.35, 0.35, 0.20, 0.10),
    q_phase1: float = 0.90,
    q_phase2: float = 0.60,
    q_phase3: float = 0.80,
    old_norse_fraction: float = 1.0 / 3.0,
) -> tuple[Array, list[PhaseConfig]]:
    """
    Define a three-phase Old English -> Middle English scenario.

    Returns the initial condition and three PhaseConfig objects.
    When old_norse_fraction == 0.0, Phase 2 is a no-contact baseline
    (n=4 labels); run_three_phase_scenario detects this and integrates
    with GRAMMARS_4 throughout.
    """
    if not (0 <= old_norse_fraction <= 1):
        raise ValueError("old_norse_fraction must lie in [0, 1].")

    x0 = validate_simplex(np.array(list(x0_oe), dtype=float))

    phase2_labels = GRAMMARS_4 if old_norse_fraction == 0.0 else GRAMMARS_5

    phases = [
        PhaseConfig("Phase 1: Old English dialects", q_phase1, a, 0.0, 20.0, 201, GRAMMARS_4),
        PhaseConfig("Phase 2: Scandinavian contact", q_phase2, a, 20.0, 80.0, 601, phase2_labels),
        PhaseConfig("Phase 3: Post-contact consolidation", q_phase3, a, 80.0, 120.0, 401, GRAMMARS_4),
    ]
    return x0, phases


def run_three_phase_scenario(
    a: float = 0.9,
    x0_oe: Iterable[float] = (0.35, 0.35, 0.20, 0.10),
    q_phase1: float = 0.90,
    q_phase2: float = 0.60,
    q_phase3: float = 0.80,
    old_norse_fraction: float = 1.0 / 3.0,
    method: str = "RK45",
) -> dict[str, object]:
    """
    Run the three-phase scenario and return concatenated results plus diagnostics.

    No-contact case (old_norse_fraction == 0.0):
        Phase 2 is integrated with n=4 throughout. G5 does not appear, so
        the mutation term cannot generate Old Norse speakers from zero.
        The output y array has 5 columns for consistency with the contact
        scenarios; column 4 (G5) is NaN throughout.

    Contact case (old_norse_fraction > 0):
        Phase 2 introduces G5 at the specified fraction. At the end of Phase
        2, G5 mass is redistributed equally among the four OE dialects
        (explicit heuristic described in §4.3 of the manuscript).
    """
    x0, phases = build_three_phase_scenario(
        a=a,
        x0_oe=x0_oe,
        q_phase1=q_phase1,
        q_phase2=q_phase2,
        q_phase3=q_phase3,
        old_norse_fraction=old_norse_fraction,
    )

    diagnostics: list[SimulationDiagnostics] = []

    # --- Phase 1: always n=4 ---
    t1, y1, d1 = simulate_segment(x0, phases[0], method=method)
    diagnostics.append(d1)

    if old_norse_fraction == 0.0:
        # No-contact baseline: remain strictly n=4 throughout Phase 2.
        # Using n=5 with G5 starting at zero would allow G5 to grow via
        # the mutation term (a*u > 0 when x5=0), which is not a valid
        # no-contact scenario.
        t2, y2, d2 = simulate_segment(validate_simplex(y1[-1]), phases[1], method=method)
        diagnostics.append(d2)

        x0_3 = validate_simplex(y2[-1])
        t3, y3, d3 = simulate_segment(x0_3, phases[2], method=method)
        diagnostics.append(d3)

        # Pad to 5 columns (G5 = NaN) for output shape consistency.
        nan_col = lambda y: np.column_stack([y, np.full(len(y), np.nan)])
        t_all = np.concatenate([t1, t2, t3])
        y_all = np.vstack([nan_col(y1), nan_col(y2), nan_col(y3)])
        phases_used = phases  # all three have GRAMMARS_4 in Phase 2 per build_three_phase_scenario

    else:
        # Contact phase: introduce G5 (Old Norse) at the specified fraction.
        x_end_1 = y1[-1]
        x0_2 = np.append((1.0 - old_norse_fraction) * x_end_1, old_norse_fraction)
        x0_2 = validate_simplex(x0_2)
        t2, y2, d2 = simulate_segment(x0_2, phases[1], method=method)
        diagnostics.append(d2)

        # Redistribute G5 mass equally among the four OE dialects.
        x_end_2 = y2[-1]
        x0_3 = x_end_2[:4] + x_end_2[4] / 4.0
        x0_3 = validate_simplex(x0_3)
        t3, y3, d3 = simulate_segment(x0_3, phases[2], method=method)
        diagnostics.append(d3)

        # Phases 1 and 3 are n=4; pad G5 column with NaN.
        nan_col = lambda y: np.column_stack([y, np.full(len(y), np.nan)])
        t_all = np.concatenate([t1, t2, t3])
        y_all = np.vstack([nan_col(y1), y2, nan_col(y3)])
        phases_used = phases

    entropy = np.array([simplex_entropy(row[~np.isnan(row)], normalized=True) for row in y_all])
    m2 = np.array([concentration_m2(row[~np.isnan(row)]) for row in y_all])

    return {
        "t": t_all,
        "y": y_all,
        "labels": GRAMMARS_5,
        "phases": phases_used,
        "diagnostics": diagnostics,
        "entropy_normalized": entropy,
        "m2": m2,
        "parameters": {
            "a": a,
            "x0_oe": tuple(x0_oe),
            "q_phase1": q_phase1,
            "q_phase2": q_phase2,
            "q_phase3": q_phase3,
            "old_norse_fraction": old_norse_fraction,
            "method": method,
        },
    }


def plot_frequencies(result: dict[str, object], output_path: str | Path | None = None) -> None:
    """Plot grammar frequencies over time."""
    t = result["t"]
    y = result["y"]
    labels = result["labels"]
    params = result["parameters"]

    plt.figure(figsize=(12, 6))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for i, (label, color) in enumerate(zip(labels, colors)):
        col = y[:, i]
        if not np.all(np.isnan(col)):
            plt.plot(t, col, label=label, color=color, linewidth=2)

    plt.axvline(20, linestyle="--", color="gray", alpha=0.7)
    plt.axvline(80, linestyle="--", color="gray", alpha=0.7)
    plt.xlabel("Time (model units)")
    plt.ylabel("Grammar frequency $x_j$")
    plt.title(
        f"a={params['a']}, q=({params['q_phase1']}, {params['q_phase2']}, {params['q_phase3']}), "
        f"Old Norse fraction={params['old_norse_fraction']:.2f}"
    )
    plt.ylim(0, 1)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(loc="best")
    plt.tight_layout()

    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")


def plot_diagnostics(result: dict[str, object], output_path: str | Path | None = None) -> None:
    """Plot entropy and M2 diagnostics on separate figures."""
    t = result["t"]

    plt.figure(figsize=(10, 4))
    plt.plot(t, result["entropy_normalized"], linewidth=2)
    plt.axvline(20, linestyle="--", color="gray", alpha=0.7)
    plt.axvline(80, linestyle="--", color="gray", alpha=0.7)
    plt.xlabel("Time (model units)")
    plt.ylabel("Normalized entropy")
    plt.title("Linguistic diversity proxy: normalized Shannon entropy")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()
    if output_path is not None:
        output_path = Path(output_path)
        plt.savefig(
            output_path.with_name(output_path.stem + "_entropy" + output_path.suffix),
            dpi=200,
            bbox_inches="tight",
        )

    plt.figure(figsize=(10, 4))
    plt.plot(t, result["m2"], linewidth=2)
    plt.axvline(20, linestyle="--", color="gray", alpha=0.7)
    plt.axvline(80, linestyle="--", color="gray", alpha=0.7)
    plt.xlabel("Time (model units)")
    plt.ylabel(r"$M_2 = \sum_j x_j^2$")
    plt.title("Coherence/concentration proxy")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()
    if output_path is not None:
        output_path = Path(output_path)
        plt.savefig(
            output_path.with_name(output_path.stem + "_m2" + output_path.suffix),
            dpi=200,
            bbox_inches="tight",
        )


def rk4_convergence_check(
    x0: Array,
    q: float,
    a: float,
    t_start: float,
    t_end: float,
    h_values: Iterable[float] = (0.1, 0.05, 0.025),
) -> list[dict[str, float]]:
    """Compare fixed-step RK4 final states across step sizes."""
    x0 = validate_simplex(x0)
    rhs = make_symmetric_lde_rhs(q=q, a=a, n=len(x0))
    outputs: list[tuple[float, Array]] = []
    for h in h_values:
        _, y = rk4_fixed_step(rhs, x0, t_start, t_end, h)
        outputs.append((h, y[-1]))

    h_ref, x_ref = outputs[-1]
    rows: list[dict[str, float]] = []
    for h, x_final in outputs:
        rows.append(
            {
                "h": float(h),
                "error_vs_smallest_h_linf": float(np.max(np.abs(x_final - x_ref))),
                "mass_error": float(abs(np.sum(x_final) - 1.0)),
                "min_component": float(np.min(x_final)),
            }
        )
    return rows


def run_named_scenarios(method: str = "RK45") -> dict[str, dict[str, object]]:
    """
    Run the five historically motivated scenarios from Table 2 of the manuscript.

    Parameter values are exploratory scenario regimes, not empirical measurements.
    """
    scenarios = {
        "canonical_mitchener_a05": dict(
            a=0.5,
            x0_oe=(0.35, 0.35, 0.20, 0.10),
            q_phase1=0.90,
            q_phase2=0.60,
            q_phase3=0.80,
            old_norse_fraction=1 / 3,
            method=method,
        ),
        "conservative_contact": dict(
            a=0.8,
            x0_oe=(0.35, 0.35, 0.20, 0.10),
            q_phase1=0.90,
            q_phase2=0.75,
            q_phase3=0.85,
            old_norse_fraction=0.15,
            method=method,
        ),
        "moderate_contact": dict(
            a=0.8,
            x0_oe=(0.35, 0.35, 0.20, 0.10),
            q_phase1=0.90,
            q_phase2=0.65,
            q_phase3=0.80,
            old_norse_fraction=0.25,
            method=method,
        ),
        "strong_contact": dict(
            a=0.75,
            x0_oe=(0.30, 0.35, 0.20, 0.15),
            q_phase1=0.90,
            q_phase2=0.55,
            q_phase3=0.75,
            old_norse_fraction=1 / 3,
            method=method,
        ),
        "no_norse_counterfactual": dict(
            a=0.8,
            x0_oe=(0.35, 0.35, 0.20, 0.10),
            q_phase1=0.90,
            q_phase2=0.90,
            q_phase3=0.90,
            old_norse_fraction=0.0,
            method=method,
        ),
    }
    return {name: run_three_phase_scenario(**params) for name, params in scenarios.items()}


_DISPLAY_NAMES: dict[str, str] = {
    "canonical_mitchener_a05": "Canonical comparison",
    "conservative_contact": "Conservative contact",
    "moderate_contact": "Moderate contact",
    "strong_contact": "Strong contact",
    "no_norse_counterfactual": "No contact",
}


def plot_scenario_comparison(
    results: dict[str, dict[str, object]],
    output_path: str | Path | None = None,
    diagnostic: str = "entropy_normalized",
) -> None:
    """
    Compare named scenarios using one scalar diagnostic over time.

    diagnostic: 'entropy_normalized' or 'm2'.
    Legend labels use publication display names from _DISPLAY_NAMES.
    """
    if diagnostic not in {"entropy_normalized", "m2"}:
        raise ValueError("diagnostic must be 'entropy_normalized' or 'm2'.")

    plt.figure(figsize=(11, 5))
    for name, result in results.items():
        display = _DISPLAY_NAMES.get(name, name)
        plt.plot(result["t"], result[diagnostic], linewidth=2, label=display)

    plt.axvline(20, linestyle="--", color="gray", alpha=0.5)
    plt.axvline(80, linestyle="--", color="gray", alpha=0.5)
    plt.xlabel("Time (model units)")
    ylabel = (
        "Normalized entropy"
        if diagnostic == "entropy_normalized"
        else r"$M_2 = \sum_j x_j^2$"
    )
    plt.ylabel(ylabel)
    plt.title(f"Scenario comparison: {ylabel}")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")


def plot_no_contact_vs_moderate_comparison(
    no_contact_result: dict[str, object],
    moderate_result: dict[str, object],
    output_path: str | Path | None = None,
) -> None:
    """
    Two-panel figure comparing no-contact baseline and moderate contact scenario.

    Left panel: four Old English grammar frequencies under no-contact
    (G5 absent throughout; column 4 is NaN and not plotted).
    Right panel: all five grammar frequencies under moderate contact,
    with G5 present only during the contact phase.

    Used in §6.3 of the manuscript to illustrate that the primary effect
    of contact is a transient displacement, not a final-state shift.
    """
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    panels = [
        (no_contact_result, "No-contact baseline"),
        (moderate_result, "Moderate contact"),
    ]

    for ax, (result, panel_title) in zip(axes, panels):
        t = result["t"]
        y = result["y"]
        labels = result["labels"]

        for i, (label, color) in enumerate(zip(labels, colors)):
            col = y[:, i]
            if not np.all(np.isnan(col)):
                ax.plot(t, col, label=label, color=color, linewidth=2)

        ax.axvline(20, linestyle="--", color="gray", alpha=0.7)
        ax.axvline(80, linestyle="--", color="gray", alpha=0.7)
        ax.set_xlabel("Time (model units)")
        ax.set_ylim(0, 1)
        ax.set_title(panel_title, fontsize=11)
        ax.grid(True, linestyle="--", alpha=0.35)
        ax.legend(loc="upper right", fontsize=8)

    axes[0].set_ylabel("Grammar frequency $x_j$")
    fig.suptitle(
        "Dashed lines mark contact-phase boundaries (t = 20 and t = 80)",
        fontsize=10,
        y=0.02,
    )
    plt.tight_layout()

    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")


def final_concentration_sweep(
    a_values: Iterable[float],
    q_contact_values: Iterable[float],
    old_norse_fraction: float = 0.25,
    x0_oe: Iterable[float] = (0.35, 0.35, 0.20, 0.10),
    q_phase1: float = 0.90,
    q_phase3: float = 0.80,
    method: str = "RK45",
) -> tuple[Array, Array, Array]:
    """
    Sweep over (a, q_contact) and return final M2 values.

    Used to locate the historically motivated scenarios in the full
    parameter space and to identify regimes of stronger disruption
    or bifurcation-like final-state divergence.
    """
    a_values = np.array(list(a_values), dtype=float)
    q_contact_values = np.array(list(q_contact_values), dtype=float)
    z = np.zeros((len(q_contact_values), len(a_values)))

    for i, q2 in enumerate(q_contact_values):
        for j, a in enumerate(a_values):
            result = run_three_phase_scenario(
                a=float(a),
                x0_oe=x0_oe,
                q_phase1=q_phase1,
                q_phase2=float(q2),
                q_phase3=q_phase3,
                old_norse_fraction=old_norse_fraction,
                method=method,
            )
            final_row = result["y"][-1]
            final_row = final_row[~np.isnan(final_row)]
            z[i, j] = concentration_m2(final_row)

    return a_values, q_contact_values, z


def plot_final_concentration_heatmap(
    a_values: Iterable[float],
    q_contact_values: Iterable[float],
    old_norse_fraction: float = 0.25,
    output_path: str | Path | None = None,
) -> tuple[Array, Array, Array]:
    """Plot a heatmap of final M2 over the (a, q_contact) parameter plane."""
    a_grid, q_grid, z = final_concentration_sweep(
        a_values=a_values,
        q_contact_values=q_contact_values,
        old_norse_fraction=old_norse_fraction,
    )

    plt.figure(figsize=(8, 5))
    image = plt.imshow(
        z,
        origin="lower",
        aspect="auto",
        extent=[a_grid.min(), a_grid.max(), q_grid.min(), q_grid.max()],
    )
    plt.colorbar(image, label=r"Final $M_2 = \sum_j x_j^2$")
    plt.xlabel(r"Cross-grammar intelligibility $a$")
    plt.ylabel(r"Contact-phase fidelity $q_2$")
    plt.title(
        f"Final concentration over parameter space (Old Norse fraction={old_norse_fraction:.2f})"
    )
    plt.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")
    return a_grid, q_grid, z


def compute_transient_displacement(
    contact_result: dict[str, object],
    no_contact_result: dict[str, object],
) -> dict[str, float]:
    """
    Transient displacement integrals D_H and D_M2 (Eqs. dH, dM2 in manuscript §5).

    D_H  = integral |H_norm_contact(t) - H_norm_no_contact(t)| dt
    D_M2 = integral |M2_contact(t)     - M2_no_contact(t)|     dt

    Both results must share the same time grid (identical phase structure).
    Uses the trapezoid rule over the full simulation window [0, 120].
    """
    t_c = contact_result["t"]
    t_nc = no_contact_result["t"]
    if t_c.shape != t_nc.shape or not np.allclose(t_c, t_nc):
        raise ValueError(
            "Time grids do not match. Both scenarios must use the same phase structure."
        )
    dH = float(np.trapezoid(
        np.abs(contact_result["entropy_normalized"] - no_contact_result["entropy_normalized"]),
        t_c,
    ))
    dM2 = float(np.trapezoid(
        np.abs(contact_result["m2"] - no_contact_result["m2"]),
        t_c,
    ))
    return {"D_H": dH, "D_M2": dM2}


if __name__ == "__main__":
    out_dir = Path(__file__).parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Running all named scenarios...")
    results = run_named_scenarios()
    moderate = results["moderate_contact"]
    no_contact = results["no_norse_counterfactual"]

    # --- Figure 1: moderate scenario frequency trajectories ---
    print("  Fig 1: moderate frequencies")
    plot_frequencies(moderate, out_dir / "lde_moderate_frequencies.png")
    plt.close("all")

    # --- Figures 2–3: moderate scenario diagnostics (entropy, M2) ---
    print("  Figs 2-3: moderate diagnostics")
    plot_diagnostics(moderate, out_dir / "lde_moderate_diagnostics.png")
    plt.close("all")

    # --- Figures 4–5: scenario comparison (entropy, M2) ---
    print("  Fig 4: scenario entropy comparison")
    plot_scenario_comparison(
        results, out_dir / "lde_scenario_entropy.png", diagnostic="entropy_normalized"
    )
    plt.close("all")

    print("  Fig 5: scenario M2 comparison")
    plot_scenario_comparison(
        results, out_dir / "lde_scenario_m2.png", diagnostic="m2"
    )
    plt.close("all")

    # --- Figure 6: two-panel no-contact vs moderate comparison ---
    print("  Fig 6: no-contact vs moderate comparison")
    plot_no_contact_vs_moderate_comparison(
        no_contact, moderate, out_dir / "lde_comparison_nocontact_vs_moderate.png"
    )
    plt.close("all")

    # --- Figure 7: parameter sweep heatmap ---
    print("  Fig 7: parameter sweep heatmap (may take a moment)...")
    plot_final_concentration_heatmap(
        a_values=np.linspace(0.3, 0.95, 20),
        q_contact_values=np.linspace(0.3, 0.9, 20),
        old_norse_fraction=0.25,
        output_path=out_dir / "lde_heatmap_m2.png",
    )
    plt.close("all")

    # --- Print final diagnostics ---
    print("\n=== Final state diagnostics ===")
    for name, result in results.items():
        e = result["entropy_normalized"][-1]
        m = result["m2"][-1]
        print(f"{name}: entropy={e:.6f}  M2={m:.6f}")
        for phase, diag in zip(result["phases"], result["diagnostics"]):
            print(
                f"  [{phase.name}]  mass_err={diag.max_mass_error:.2e}"
                f"  min_x={diag.min_component:.6f}"
            )

    print(f"\nAll figures written to: {out_dir.resolve()}")

    # --- Transient displacement integrals (manuscript Table 4) ---
    print("\nTransient displacement integrals D_H and D_M2")
    print("(copy values into Table tab:displacement-integrals in borrador_v6_master_internal_review.tex)")
    contact_scenarios = [
        ("canonical_mitchener_a05", "Canonical comparison"),
        ("conservative_contact",    "Conservative contact"),
        ("moderate_contact",        "Moderate contact"),
        ("strong_contact",          "Strong contact"),
    ]
    for key, display in contact_scenarios:
        disp = compute_transient_displacement(results[key], no_contact)
        print(f"  {display}: D_H = {disp['D_H']:.4f},  D_M2 = {disp['D_M2']:.4f}")

    # --- RK4 convergence check (Phase 1 representative) ---
    print("\nRK4 convergence check (Phase 1, moderate scenario):")
    rows = rk4_convergence_check(
        x0=np.array([0.35, 0.35, 0.20, 0.10]),
        q=0.90,
        a=0.8,
        t_start=0.0,
        t_end=20.0,
        h_values=(0.1, 0.05, 0.025),
    )
    for row in rows:
        print(row)
