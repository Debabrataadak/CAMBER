# CAMBER: <code><b>CAMB</b> with <b>E</b>volving <b>R</b>otation</code>

CAMBER is a modified version of [CAMB](https://github.com/cmbant/CAMB) 2.0.5 (Lewis & Challinor) that computes the
CMB power spectra in the presence of **cosmic birefringence**, a rotation of the CMB polarization plane by an
angle $\beta$. The angle can be constant or vary with time, $\beta(z)$, for example when it is produced by an
axion-like field that evolves during recombination.

The code was developed for [paper title, arXiv:XXXX.XXXXX].

## What is added

- **Constant rotation**: a single angle $\beta$ (`birefringence_angle`, in degrees).
- **Time-varying rotation**: a table $\beta(z)$ (`set_birefringence_table`). At each time step of the
  line-of-sight integration, the polarization source is rotated by $\beta$ at that redshift, so the rotation is
  correct when $\beta$ changes during recombination or reionization.
- **Axion-like field**: `camb.birefringence.axion_birefringence_table` solves
  $\ddot\phi + 3H\dot\phi + m^2\phi = 0$ for a quadratic potential in the background cosmology and returns
  $\beta(z) = \beta_{\rm ini}[f(z) - f(0)]$, with $f = \phi/\phi_{\rm ini}$ (convention of Namikawa, Murai &
  Naokawa, arXiv:2506.20824).
- **Spectra**: `CAMBdata.get_birefringence_cls` returns TT, EE, BB, TE, EB and TB, unlensed or lensed. Lensing
  uses CAMB's curved-sky correlation-function method (Challinor & Lewis 2005).

Everything else is unchanged CAMB.

## Installation

A Fortran compiler is needed (tested with gfortran 15.2 on macOS).

```bash
git clone https://github.com/<user>/CAMBER.git
cd CAMBER
pip install -e . --no-build-isolation
```

Check that the modified version is the one loaded:

```python
import camb
from camb.birefringence import axion_birefringence_table
print(camb.__file__)
```

## Usage

```python
import camb
from camb.birefringence import axion_birefringence_table

pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.965, lmax=3000)

# Constant rotation, beta = 0.3 deg
pars.birefringence_angle = 0.3
cl = camb.get_results(pars).get_birefringence_cls(CMB_unit='muK', lensed=True)
# cl[:, 0:6] = TT, EE, BB, TE, EB, TB  (D_ell = ell(ell+1)C_ell/2pi in muK^2)

# Axion-like field with m = 1e-28 eV and beta_ini = 0.3 deg
pars.birefringence_angle = 0.0
z, beta = axion_birefringence_table(pars, mass_ev=1e-28, beta_ini=0.3)
pars.set_birefringence_table(z, beta)
cl = camb.get_results(pars).get_birefringence_cls(CMB_unit='muK', lensed=True)
```

A table set with `set_birefringence_table` replaces `birefringence_angle`; call it with no arguments to remove
the table.

## Changes with respect to CAMB

Only 11 files differ from CAMB 2.0.5, and every changed or added line is marked with a comment
`! modified` (Fortran) or `# modified` (Python):

| File | Change |
|---|---|
| `fortran/cmbmain.f90` | rotation of the polarization source at each line-of-sight step |
| `fortran/lensing.f90` | lensing of the rotated BB, EB, TB spectra |
| `fortran/model.f90`, `fortran/results.f90`, `fortran/camb.f90`, `fortran/camb_python.f90` | parameters, storage and Python interface |
| `camb/model.py`, `camb/results.py`, `camb/_ini.py` | Python parameters and `get_birefringence_cls` |
| `camb/birefringence.py` | new: axion field solver and $\beta(z)$ table |
| `camb/tests/camb_test.py` | tests |

To list all changes:

```bash
grep -rn -E '(!|#) *modified' camb fortran
```

All other files, including the upstream development files (`.devcontainer/`, `AGENTS.md`, `CONTRIBUTING.md`,
`CHANGELOG.md`), are unchanged from CAMB.

## Citation

If you use this code, please cite:

- [Authors], *Constraints on the mass of an axion-like field from the shape of the parity-violating signal in the
  ACT DR6 polarization spectra*, arXiv:XXXX.XXXXX (this code)
- T. Namikawa, K. Murai & F. Naokawa, *Planck constraints on axion-like particles through isotropic cosmic
  birefringence*, arXiv:2506.20824 (time-dependent rotation and the convention for $\beta(z)$)
- A. Lewis, A. Challinor & A. Lasenby, ApJ 538, 473 (2000), arXiv:astro-ph/9911177 (CAMB)
- A. Challinor & A. Lewis, Phys. Rev. D 71, 103010 (2005), arXiv:astro-ph/0502425 (lensed spectra)

BibTeX:

```bibtex
@article{ThisPaper,
  author        = "[Authors]",
  title         = "{Constraints on the mass of an axion-like field from the shape of the parity-violating signal in the ACT DR6 polarization spectra}",
  eprint        = "XXXX.XXXXX",
  archivePrefix = "arXiv",
  primaryClass  = "astro-ph.CO",
  year          = "2026"
}

@article{NMN2025,
  author        = "Namikawa, Toshiya and Murai, Kai and Naokawa, Fumihiro",
  title         = "{Planck constraints on axion-like particles through isotropic cosmic birefringence}",
  eprint        = "2506.20824",
  archivePrefix = "arXiv",
  primaryClass  = "astro-ph.CO",
  year          = "2025"
}

@article{CAMB,
  author        = "Lewis, Antony and Challinor, Anthony and Lasenby, Anthony",
  title         = "{Efficient computation of CMB anisotropies in closed FRW models}",
  journal       = "Astrophys. J.",
  volume        = "538",
  pages         = "473",
  year          = "2000",
  doi           = "10.1086/309179",
  eprint        = "astro-ph/9911177",
  archivePrefix = "arXiv"
}

@article{Challinor2005,
  author        = "Challinor, Anthony and Lewis, Antony",
  title         = "{Lensed CMB power spectra from all-sky correlation functions}",
  journal       = "Phys. Rev. D",
  volume        = "71",
  pages         = "103010",
  year          = "2005",
  doi           = "10.1103/PhysRevD.71.103010",
  eprint        = "astro-ph/0502425",
  archivePrefix = "arXiv"
}
```

## Licence

CAMBER is distributed under the CAMB licence (`LICENCE.txt`, LGPL with additional conditions). The CAMB
copyright notice and licence apply to this modified version. Among the conditions: publications using results
of the code must be submitted to arXiv at the same time as, or before, submission to a journal.
