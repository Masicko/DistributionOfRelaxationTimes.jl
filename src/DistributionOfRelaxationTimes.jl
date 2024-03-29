module DistributionOfRelaxationTimes

# [ ] radeji udelat obsirnejsi strukturu, kterou pak budou moct vyuzivat
#       ... DRT_tools
#         ... DistributionOfRelaxationTimes
#         ... Z_file_IO
#         ... Nyquist_plots


# TODO
# [ ] KK relations 
# [ ] separate DRT vs EIS_reconstruction vs peak_analysis in calling get_DRT

using PyPlot
using DataFrames

using NonNegLeastSquares
using LinearAlgebra

include("DRT_base.jl")
export DRT_control

include("DRT_peak_analysis.jl")
include("DRT_evaluation.jl")
export get_DRT

include("DRT_plotting.jl")
export plot_Nyquist
export plot_DRT_h
export plot_DRT_RC
export plot_DRT_Rtau
export plot_DRT_Rf

include("Z_file_IO.jl")
export read_Z_file
end # module
