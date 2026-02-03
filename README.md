# Lugiato-Lefever Equation (LLE) Solver

A high-performance Python solver for the normalized Lugiato-Lefever equation using the leapfrog integration method. This repository provides both optimized and minimal implementations for flexible use.

## Description

### The Lugiato-Lefever Equation

The Lugiato-Lefever equation (LLE) is a fundamental model in nonlinear optics describing the dynamics of optical cavities with a single pumping wave. The normalized form is:

$$\frac{\partial A}{\partial t} = -(1 - i\delta)A - i|A|^2 A + s$$

Where:
- **A(τ,t)**: Complex field envelope (depends on space τ and time t)
- **δ**: Detuning parameter (normalized frequency offset)
- **s**: Pump strength (amplitude of the driving field)
- The equation describes the balance between:
  - **Cavity losses** and detuning: $-(1 - i\delta)A$
  - **Self-phase modulation** (Kerr nonlinearity): $-i|A|^2 A$
  - **Pump driving**: $+s$

### Physical Significance

The LLE is essential for understanding:
- **Optical frequency combs**: Generation of equally-spaced spectral lines
- **Temporal solitons**: Stable light pulses in nonlinear systems
- **Nonlinear wave dynamics**: Chaotic and turbulent regimes
- **Cavity QED**: Quantum effects in optical systems

## Files Overview

### `lle_3_2_leapfrog.py` (Optimized Version)
**Best for**: Production use, large simulations, best performance

**Features:**
- ✓ Numba JIT compilation for 10-50x speedup
- ✓ Intel MKL-FFT (when available) for faster Fourier transforms
- ✓ Progress bars with `tqdm`
- ✓ Publication-quality plots with `scienceplots`
- ✓ Pre-calculated dispersion operators for efficiency

**Requirements:**
```
numpy matplotlib scipy numba tqdm scienceplots mkl-fft (optional)
```

**Runtime:** ~1-5 seconds for 1000 round trips (typical CPU)

---

### `lle_3_2_leapfrog_minimal.py` (Minimal Version)
**Best for**: Quick tests, simple environments, educational use

**Features:**
- ✓ Only NumPy and Matplotlib
- ✓ Pure Python (no JIT compilation)
- ✓ Easy to understand and modify
- ✓ Same scientific accuracy as optimized version
- ✓ Works everywhere Python works

**Requirements:**
```
numpy matplotlib
```

**Runtime:** ~+5 seconds for 1000 round trips (typical CPU)

---

## Installation

### Option 1: Optimized Version (Recommended)

```bash
# Using pip
pip install numpy matplotlib scipy numba tqdm scienceplots

# Using conda
conda install -c conda-forge numpy matplotlib scipy numba tqdm scienceplots
```

### Option 2: Minimal Version Only

```bash
pip install numpy matplotlib
```

### Option 3: Both Versions

```bash
# All dependencies (includes both versions)
pip install numpy matplotlib scipy numba tqdm scienceplots
```

## Usage

### Quick Start - Minimal Version

```bash
python lle_3_2_leapfrog_minimal.py
```

No dependencies beyond numpy and matplotlib!

### Optimized Version

```bash
python lle_3_2_leapfrog.py
```

Much faster for large simulations.

### Expected Output

Both versions produce:
1. **Console output** showing progress and simulation time
2. **2×2 Figure** with four panels:
   - **(a) Trajectory**: Space-time evolution of field intensity (color map)
   - **(b) Final Intensity**: $|A|^2$ profile at end of simulation
   - **(c) Spectrum**: Frequency content in dB (log scale)
   - **(d) Real & Imaginary**: Real and imaginary parts of final field

## Parameters

Edit these in the `main()` function:

```python
# Spatial discretization
n_points = 2**10            # Number of spatial points (1024)
domain_half = 40            # Half domain width [-40, 40]

# Time integration
time_step = 1e-2            # Δt for leapfrog steps
n_roundtrips = 250          # Total cavity round trips

# Saving/Plotting
n_save = 10                 # Save every 10 round trips

# Physics parameters
delta = 3.0                 # Detuning parameter
pump = 2.0                  # Pump strength
```

### Physical Parameter Guide

| Parameter | Range | Effect |
|-----------|-------|--------|
| `delta` (δ) | -10 to 10 | Cavity detuning; controls oscillation frequency |
| `pump` (s) | 0 to 5 | Driving strength; higher → more energy → chaos |
| `time_step` | 1e-3 to 1e-2 | Smaller → more accurate but slower |
| `n_points` | 2^8 to 2^12 | Grid resolution; higher → more detail but slower |

### Common Scenarios

**Stable Soliton Formation:**
```python
delta = 3.0
pump = 2.0
n_roundtrips = 500
```

**Chaotic/Turbulent Dynamics:**
```python
delta = 1.0
pump = 4.0
n_roundtrips = 1000
```

**Quick Test:**
```python
n_points = 2**9
time_step = 1e-1
n_roundtrips = 50
```

## The Leapfrog Integration Method

### How It Works

The leapfrog algorithm splits the LLE into:

1. **Nonlinear part** (real space):
   $$\frac{\partial A}{\partial t} = -i|A|^2 A + s$$

2. **Linear part** (Fourier space):
   $$\frac{\partial A}{\partial t} = -(1 - i\delta)A$$

Each round trip:
```
A = func(A, s)           # Half-step nonlinear
A_k = FFT(A)             # Transform to Fourier
A_k = exp(L*h) * A_k     # Apply linear operator
A = iFFT(A_k)            # Transform back
A = func(A, s)           # Half-step nonlinear
```

### Why Leapfrog?

- ✓ **2nd-order accurate** (better than forward Euler)
- ✓ **Symplectic integrator** (energy conserving for Hamiltonian systems)
- ✓ **Simple to implement**
- ✓ **Stable for stiff problems**

## Performance Comparison

| Version | Optimizations | Runtime (250 RT) | Memory |
|---------|---------------|------------------|--------|
| Optimized | Numba + MKL-FFT | ~2-5 sec | ~64 MB |
| Minimal | None | ~30-60 sec | ~64 MB |

**Speedup**: 10-30x faster with optimizations

### Optimization Details

**Optimized Version:**
- Numba JIT compiles `func()` and `calc_disp_lf()` to machine code
- Dispersion operators computed once, reused 250 times
- Possible MKL-FFT integration for faster transforms

**Minimal Version:**
- Pure Python loops
- NumPy's built-in FFT (still fast)
- Same algorithm, less overhead

## Understanding the Output

### (a) Trajectory Panel
- **X-axis**: Spatial coordinate τ (normalized)
- **Y-axis**: Round trip number
- **Color**: Intensity $|A|^2$ in dB (normalized to peak)
- **Interpretation**: Watch how the field evolves, solitons form/stabilize, or chaos emerges

### (b) Final Intensity
- **X-axis**: Spatial coordinate τ
- **Y-axis**: Intensity $|A|^2$
- **Interpretation**: Shape of field at simulation end (usually sech-like soliton or complex structure)

### (c) Spectrum
- **X-axis**: Frequencies (Fourier components)
- **Y-axis**: Power spectral density (dB)
- **Interpretation**: Frequency comb structure; discrete peaks indicate soliton trains

### (d) Real & Imaginary
- **X-axis**: Spatial coordinate τ
- **Y-axis**: Value
- **Red (Re)**: Real part of A
- **Green (Im)**: Imaginary part of A
- **Interpretation**: Both components contribute to the physical field

## Code Structure

### Main Functions

```python
func(A, s)                    # Nonlinear term: i|A|^2*A + s
calc_disp_lf(qx2, delta, h)  # Dispersion operator matrices
leapfrog(A, h, s, eq, ...)   # One round trip integration
loop(n, xx, h, rt, nt, ...)  # Full simulation
_plot_results(...)            # Generate 2×2 figure
main()                        # Entry point
```

### Data Flow

```
main()
  ↓
loop() - Integration loop
  ├→ calc_disp_lf() - Pre-calculate operators
  ├→ for each round trip:
  │   └→ leapfrog() - One integration step
  │       └→ func() - Nonlinear part
  │       └→ FFT/iFFT - Linear part
  └→ returns x, A_final, qx, A_history
  ↓
_plot_results() - Visualization
  └→ Show 2×2 figure
```

## Examples & Tips

### Example 1: Find Soliton Region
```python
# Scan pump values to find solitons
for pump_val in [1.5, 2.0, 2.5, 3.0]:
    # Modify pump = pump_val in main()
    main()
```

### Example 2: Higher Resolution
```python
n_points = 2**12     # 4096 points instead of 1024
time_step = 5e-3     # Finer time stepping
n_roundtrips = 1000
```

### Example 3: Check Convergence
Run with `n_points = 2**9` and `2**10` to verify results don't change.

### Tips for Success
- Start with default parameters (they work!)
- If plot looks weird, increase `n_points` or decrease `time_step`
- For chaos/complex dynamics, increase `n_roundtrips` (need more iterations to settle)
- Monitor the simulation time - if too slow, reduce `n_points`

## Scientific Background

### Frequency Combs

When `pump` exceeds a threshold, parametric instabilities generate equally-spaced spectral lines - a **frequency comb**:
- Applications: Precision spectroscopy, optical clocks, LiDAR
- Observable in the spectrum plot (c)

### Soliton Solutions

Stable localized structures emerge at specific parameter regions:
- Profile: hyperbolic secant $\operatorname{sech}(x)$
- See panel (b) for the soliton envelope

### Chaos & Turbulence

High pump values lead to complex dynamics:
- Trajectory (a) shows intricate structure
- Spectrum (c) has continuous background noise

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Too slow | Use `lle_3_2_leapfrog.py` (optimized) |
| Out of memory | Reduce `n_points` (e.g., 2^9 instead of 2^10) |
| Weird plot | Increase `n_points` or decrease `time_step` |
| Flat intensity profile | Increase `pump` or try different `delta` |
| Segfault on import | Update numba: `pip install --upgrade numba` |

## Advanced Usage

### Modify Initial Conditions

Edit the `loop()` function:
```python
# Instead of sech soliton:
A1 = np.random.randn(n)  # Random noise
A1 = np.sin(np.pi * x / xx)  # Sinusoidal
```

### Custom Visualization

Modify `_plot_results()` to add:
- Energy evolution plots
- Spectral animation
- Phase portraits

### Batch Simulations

```python
for delta_val in np.linspace(0, 5, 10):
    # Run with different delta values
    # Save results to file
    pass
```

## Performance Tips

1. **For fast iteration**: Use minimal version during development
2. **For production**: Use optimized version with Numba
3. **Parallelization**: Run multiple parameter sweeps in parallel


### Key Papers
1. Parra-Rivas, P., Gomila, D., Gelens, L. & Knobloch, E. Bifurcation structure of localized states in the Lugiato-Lefever equation with anomalous dispersion. Phys. Rev. E 97, 042204 (2018).


### Related Topics
- Nonlinear optics & optical cavities
- Frequency comb generation
- Soliton physics
- Chaos in dynamical systems

## Author & Citation

**Repository**: https://github.com/ysarri/Leapfrog-LLE

If you use this code, please cite:
```bibtex
@software{lle_solver_2026,
  title={Leapfrog-LLE: Lugiato-Lefever Equation Solver},
  author={Yelo-Sarr\'{i}on, Jes\'{u}s},
  year={2026},
  url={https://github.com/ysarri/Leapfrog-LLE}
}
```

Or in plain text:
```
Yelo-Sarrión, J. (2026). Leapfrog-LLE: Lugiato-Lefever Equation Solver. 
Retrieved from https://github.com/ysarri/Leapfrog-LLE
```

## License

MIT License - See LICENSE file for details

## Contact

For questions or issues, please open a GitHub issue at: https://github.com/ysarri/Leapfrog-LLE/issues

---

**Repository**: [github.com/ysarri/Leapfrog-LLE](https://github.com/ysarri/Leapfrog-LLE)  
**Last Updated**: February 3, 2026  
**Python Version**: 3.12 
**Status**: Active Development
