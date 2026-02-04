"""Lugiato-Lefever Equation (LLE) Solver using Leapfrog Method - Minimal Version.

Minimal version requiring only numpy and matplotlib.

This module simulates the normalized Lugiato-Lefever equation with pump:
    dA/dt = -(1 - i*delta)*A - i*|A|^2*A + id^2/dx^2A + s

Where:
    - A: complex field envelope
    - s: pump strength
    - delta: detuning parameter

03/02/2026 Jesús Yelo-Sarrión

"""

import numpy as np
import matplotlib.pyplot as plt


def func(A: np.ndarray, s: float) -> np.ndarray:
    """Compute nonlinear terms of the LLE.

    Args:
        A: Complex field amplitude
        s: Pump strength

    Returns:
        Nonlinear contribution i*|A|^2*A + s
    """
    return 1j * A * np.conj(A) * A + s


def calc_disp_lf(qx2: np.ndarray, delta: float, h: float) -> tuple:
    """Calculate dispersion matrices for leapfrog integration.

    Args:
        qx2: Squared frequency array (-i*q^2/2)
        delta: Detuning parameter
        h: Time step

    Returns:
        Tuple of (eq, feq, eq2, feq2) matrices for leapfrog scheme
    """
    disp = qx2 - 1j * delta - 1  # Dispersion operator
    eq = np.exp(disp * h)
    feq = (eq - 1) / disp
    eq2 = np.exp(disp * 2 * h)
    feq2 = (eq2 - 1) / disp
    return eq, feq, eq2, feq2


def leapfrog(
    A: np.ndarray,
    h: float,
    s: float,
    eq: np.ndarray,
    feq: np.ndarray,
    eq2: np.ndarray,
    feq2: np.ndarray,
    func=func,
) -> np.ndarray:
    """Leapfrog integration step for one round trip.

    Args:
        A: Field at current time
        h: Time step
        s: Pump strength
        eq, feq, eq2, feq2: Dispersion matrices
        func: Nonlinear function

    Returns:
        Field after one round trip
    """
    t = 0.0
    A_f = np.fft.fft(A)
    while t <= 1:
        # Half step nonlinear
        A = func(A, s)
        A = np.fft.fft(A)

        # Full step linear in Fourier space
        A = eq * A_f + feq * A
        A = np.fft.ifft(A)

        # Half step nonlinear
        A = func(A, s)
        A = np.fft.fft(A)

        # Update field
        A_f = eq2 * A_f + feq2 * A
        A = np.fft.ifft(A_f)
        t += 2 * h
    return A


def loop(
    n: int,
    xx: float,
    h: float,
    rt: int,
    nt: int,
    delta: float,
    s: float,
    method: str = "lf",
) -> tuple:
    """Run LLE simulation for multiple round trips.

    Args:
        n: Number of spatial points
        xx: Half of spatial domain
        h: Time step size
        rt: Total number of round trips
        nt: Save solution every nt round trips
        delta: Detuning parameter
        s: Pump strength
        method: Integration method name

    Returns:
        Tuple of (x, A_final, qx, A_history) arrays
    """
    # Spatial grid setup
    x = np.linspace(-xx, xx, n)
    dx = abs(x[1] - x[0])
    qx = 2 * np.pi * np.fft.fftfreq(n, dx)
    qx2 = -1j * qx * qx  # Squared frequency for dispersion

    # Initial condition: sech soliton profile
    A1 = np.sqrt(2 * np.abs(delta)) / np.cosh(np.sqrt(np.abs(delta)) * x)

    # noise level (eg: 0.01 -> 1%)
    noise_level = 0.1 
    noise = noise_level * (np.random.randn(n) + 1j * np.random.randn(n))
    A1 = A1 + noise

    # Storage for trajectory
    A_ring = np.zeros((rt // nt, n), dtype=np.complex128)

    # Pre-calculate dispersion operators (constant throughout integration)
    eq, feq, eq2, feq2 = calc_disp_lf(qx2, delta, h)

    # Integration loop
    print(f"Running {rt} round trips with {method} method...")
    for j in range(1, rt + 1):
        # Simple progress indicator
        if j % max(1, rt // 10) == 0:
            print(f"  Round trip {j}/{rt}")

        A1 = leapfrog(A1, h, s, eq, feq, eq2, feq2)

        # Save solution every nt round trips
        if j % nt == 0:
            A_ring[j // nt - 1, :] = A1

    return x, A1, qx, A_ring


def _plot_results(
    x: np.ndarray,
    qx: np.ndarray,
    intensity: np.ndarray,
    final_field: np.ndarray,
    n_rounds: int,
    domain_half: float,
    delta: float,
    pump: float,
    time_step: float,
) -> None:
    """Plot simulation results in 2x2 grid.

    Args:
        x: Spatial coordinate array
        qx: Frequency array
        intensity: Time-space intensity matrix
        final_field: Final field profile
        n_rounds: Total round trips plotted
        domain_half: Half domain width
        delta: Detuning parameter
        pump: Pump strength
        time_step: Integration time step
    """
    # Figure setup with golden ratio
    fig_width_pt = 246.0
    inches_per_pt = 1.0 / 72.27
    golden_mean = (np.sqrt(5.0) - 1.0) / 2.0
    fig_width = 1.9 * fig_width_pt * inches_per_pt
    fig_height = 1.2 * fig_width * golden_mean

    fig = plt.figure(figsize=(fig_width, fig_height))
    plt.subplots_adjust(
        wspace=0.3, hspace=0.4, bottom=0.125, top=0.95, left=0.125, right=0.98
    )

    # Panel (a): 2D space-time trajectory
    ax1 = fig.add_subplot(2, 2, 1)
    Z = 10 * np.log10(intensity / np.max(intensity))
    yticks = np.linspace(0, n_rounds, Z.shape[0])
    cmesh = plt.pcolormesh(
        x, yticks, Z, cmap="Blues", shading="gouraud", rasterized=True
    )
    plt.ylabel(r"$\mathrm{Round-trip}$")
    plt.xlabel(r"$\tau$")
    plt.title(r"$\mathrm{(a)~Trajectory}$")
    plt.colorbar(cmesh, label="Intensity (dB)")

    # Panel (b): Final intensity profile
    ax2 = fig.add_subplot(2, 2, 2)
    plt.plot(x, np.abs(final_field) ** 2, linewidth=1.5, color="C0")
    plt.xlim(-domain_half, domain_half)
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$|A|^2$")
    plt.title(r"$\mathrm{(b)~Final~intensity}$")
    plt.grid(True, alpha=0.3)

    # Panel (c): Frequency spectrum
    ax3 = fig.add_subplot(2, 2, 3)
    sp = np.fft.fftshift(np.fft.fft(final_field))
    spectrum_normalized = np.abs(sp / sp.max()) ** 2
    spp = 10 * np.log10(np.clip(spectrum_normalized, 1e-20, 1.0))
    plt.plot(np.fft.fftshift(qx), spp, linewidth=1.5, color="C1")
    plt.ylabel(r"$\mathrm{Spectrum~(dB)}$")
    plt.xlabel(r"$\mathrm{Frequencies}$")
    plt.title(r"$\mathrm{(c)~Spectrum}$")
    plt.grid(True, alpha=0.3)

    # Panel (d): Re vs im profile
    ax4 = fig.add_subplot(2, 2, 4)
    plt.plot(x, np.real(final_field), linewidth=1.5, color="C2", label="Re")
    plt.plot(x, np.imag(final_field), linewidth=1.5, color="C3", label="Im")
    plt.xlim(-domain_half, domain_half)
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$\mathrm{Value}$")
    plt.title(r"$\mathrm{(d)~Real~and~Imaginary}$")
    plt.legend()
    plt.grid(True, alpha=0.3)

    fig.suptitle(rf"$\delta={delta}$, $s={pump}$, $h={time_step:.0E}$", fontsize=12)
    plt.tight_layout()
    plt.show()


def main() -> None:
    """Main simulation runner."""
    # Simulation parameters
    n_points = 2**10  # Number of spatial points
    time_step = 1e-2  # Integration time step
    domain_half = 40  # Half spatial domain width
    n_roundtrips = 1000  # Total round trips
    n_save = n_roundtrips//100  # Save every n_save round trips
    delta = 3.0  # Detuning parameter
    pump = 2.0  # Pump strength
    method = "lf"  # Integration method

    # Run simulation
    print("\n" + "=" * 60)
    print("Lugiato-Lefever Equation Solver - Leapfrog Method (Minimal)")
    print("=" * 60)
    print(f"Parameters: delta={delta}, s={pump}")
    print(f"Domain: [{-domain_half}, {domain_half}], n={n_points}")
    print(f"Round trips: {n_roundtrips}, Time step: {time_step}")
    print("=" * 60 + "\n")

    import time

    t0 = time.time()
    x, A_final, qx, A_ring = loop(
        n_points, domain_half, time_step, n_roundtrips, n_save, delta, pump, method
    )
    elapsed = time.time() - t0
    print(f"\nSimulation completed in {elapsed:.2f} seconds\n")

    # Compute results
    intensity = np.abs(A_ring) ** 2
    n_saved_rounds = intensity.shape[0] * n_save
    final_field = A_ring[-1, :]

    # Create visualization
    _plot_results(
        x,
        qx,
        intensity,
        final_field,
        n_saved_rounds,
        domain_half,
        delta,
        pump,
        time_step,
    )


if __name__ == "__main__":
    main()

