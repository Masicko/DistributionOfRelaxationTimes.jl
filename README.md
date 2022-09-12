# DistributionOfRelaxationTimes.jl

<info about the package>

## Installation
The package can be installed using 
```julialang
] add https://github.com/Masicko/DistributionOfRelaxationTimes.jl
```

## Usage
Given the frequencies `f_range` and impedance data `Z_range`, one can call a function 
```julialang
DRT_data = get_DRT(f_range, Z_range)
```

or specify the details by `DRT_control` struct

```julialang
DRT_data = get_DRT(f_range, Z_range, DRT_control(lambda=0.4, tau_range_fac=3.0))
```

The control parameters are

- `f_range` : desired range of reconstructed impedance data.
- `lambda` : smoothing (regularization parameters) - better conditioned problem, but further from the original data.
- `tau_max_fac` : adjusting maximal value of tau for DRT (value 1.0 is somehow a normal electrochemical value).
- `tau_min_fac`: adjusting minimal value of tau for DRT (value 1.0 is somehow a normal electrochemical value).
- `tau_range_fac` : how many tau nodes for DRT.
- `peak_merge_tol` : merges two peaks with characteristic time closer then this value

Output is stored in `DRT_struct` which has fields

- `EIS_df` : storing reconstructed impedance data from DRT approximation.
- `tau_range` : domain of DRT function
- `h` : DRT function
- `R_ohm` : fitted ohmic resistance
- `L` : fitted inductance
- `peaks_df` : DataFrame with results of peak analysis having each peak in each row with columns `tau_c` (characteristic frequency of appripriate RC element), `R` resistance and `C` capacitance.
- `control` : the control struct used for DRT evaluation.

## Plotting
Data can be plotted with functions which are shown in examples.

