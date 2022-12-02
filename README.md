# DistributionOfRelaxationTimes.jl

<info about the package>

## Installation
The package can be installed using 
```julialang
] add https://github.com/Masicko/DistributionOfRelaxationTimes.jl
```

## Usage
Given the frequencies `f_list` and impedance data `Z_list`, one can call a function 
```julialang
DRT_data = get_DRT(f_list, Z_list)
```

or using DataFrame e.g. `EIS_df` with columns `f` frequency and `Z` impedance
```julialang
DRT_data = get_DRT(EIS_df)
```

Details can be specified by `DRT_control` struct, for example

```julialang
DRT_data = get_DRT(f_list, Z_list, DRT_control(lambda=0.4, tau_sampling_fac=3.0))
```

The control parameters are

- `lambda` : smoothing (regularization parameters) - better conditioned problem, but further from the original data.
- `tau_max_fac` : adjusting maximal value of tau for DRT (1 / (2*pi*EIS_df.f[end]) / control.tau_min_fac).
- `tau_min_fac`: adjusting minimal value of tau for DRT (value 1.0 means that the lowest tau is a dicrect counterpart to the highest frequency in EIS data).
- `tau_max_abs`: absolute value for the highest tau in DRT analysis (with no respect to EIS data)
- `tau_min_abs`: absolute value for the lowest tau in DRT analysis (with no respect to EIS data)
- `tau_sampling_fac` : how many tau nodes for DRT.
- `peak_merge_tol` : merges two peaks with characteristic time closer then this value
- `f_list` : desired list of frequencies for reconstructed impedance data in frequency domain.

Output is stored in `DRT_struct` which has fields

- `EIS_df` : a DataFrame storing reconstructed impedance data from DRT approximation.
- `tau_list` : a domain of DRT function
- `h` : DRT function
- `R_ohm` : fitted ohmic resistance
- `L` : fitted inductance
- `peaks_df` : DataFrame with results of peak analysis having each peak in each row with columns `tau_c` (characteristic frequency of appripriate RC element), `R` resistance and `C` capacitance.
- `control` : the control struct used for DRT evaluation.

## Plotting
Data can be plotted with functions which are shown in Example_.

