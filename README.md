# Leapfrog-LLE: Lugiato-Lefever Equation Solver

A Python solver for the normalized Lugiato-Lefever equation using the leapfrog integration method.

## Quick Start

### Installation

**Minimal version** (numpy + matplotlib only):
```bash
pip install numpy matplotlib
python lle_3_2_leapfrog_minimal.py
```

**Optimized version** (faster):
```bash
pip install numpy matplotlib scipy numba tqdm scienceplots
python lle_3_2_leapfrog.py
```

## What It Does

Simulates the Lugiato-Lefever equation:
$$\frac{\partial A}{\partial t} = -(1 - i\delta)A + i\frac{\partial^2 A}{\partial x^2} - i|A|^2 A  + s$$

- **A**: Complex field envelope
- **δ (delta)**: Detuning parameter
- **s (pump)**: Pump strength

Outputs a 2×2 figure showing:
- **(a)** Space-time trajectory
- **(b)** Final intensity profile  
- **(c)** Frequency spectrum
- **(d)** Real and imaginary parts

## Customize

Edit parameters in `main()`:

```python
n_points = 2**10      # Grid resolution (1024)
delta = 3.0           # Detuning
pump = 2.0            # Pump strength
n_roundtrips = 250    # Simulation length
time_step = 1e-2      # Integration step
```

## Two Versions

| Feature | lle_3_2_leapfrog.py | lle_3_2_leapfrog_minimal.py |
|---------|---------------------|---------------------------|
| Speed | ~2-5 sec | ~+5 sec |
| Dependencies | numpy, matplotlib, scipy, numba, tqdm | numpy, matplotlib |
| Optimization | Numba JIT + MKL-FFT | None |

## Citation

```bibtex
@software{lle_solver_2026,
  title={Leapfrog-LLE: Lugiato-Lefever Equation Solver},
  author={Yelo-Sarrión, Jesús},
  year={2026},
  url={https://github.com/ysarri/Leapfrog-LLE}
}
```

## License

MIT License

## References

1. Parra-Rivas, P., Gomila, D., Gelens, L. & Knobloch, E. Bifurcation structure of localized states in the Lugiato-Lefever equation with anomalous dispersion. Phys. Rev. E 97, 042204 (2018).
2. Montagne, R., Hernández-García, E., Amengual, A. & San Miguel, M. Wound-up phase turbulence in the complex Ginzburg-Landau equation. Phys. Rev. E 56, 151–167 (1997).
3. Yelo Sarrión, J. Nonlinear Dynamics in driven-dissipative photonic dimers. https://difusion.ulb.ac.be/vufind/Record/ULB-DIPOT:oai:dipot.ulb.ac.be:2013/360823/Holdings (2023).

