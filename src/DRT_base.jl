Base.@kwdef mutable struct DRT_control
  lambda::Float32 = 0.0
  #
  tau_min_fac::Float32 = 10.0
  tau_max_fac::Float32 = 10.0
  tau_sampling_fac::Float32 = 100.0
  tau_min_abs::Union{Nothing, Float32} = nothing
  tau_max_abs::Union{Nothing, Float32} = nothing
  #
  peak_merge_tol::Float32 = 0.0
  f_range = nothing
  #
  HF_adjustment::Bool = true
end

Base.@kwdef mutable struct DRT_struct
  EIS_df::DataFrame
  tau_list::Array{Float64}
  h::Array{Float64}
  R_ohm::Float64
  L::Float64
  peaks_df::DataFrame
  control::DRT_control
end

function geomspace(A, B; N)
  res = zeros(0)
  if A == B
    append!(res, A)
    return res
  end
  if N <= 1 || A > B
    println("ERROR: N <= 1 || A > B")
    return throw(Exception)
  end
  for i in 1:N
    append!(res, A*((Float64(B)/A)^((i - 1)/(N - 1))))
  end
  res
end

function geomspace_by_fac(A, B, q_fac)
  # frequency nodes to compare
  # must be sorted upwards
  f_list = zeros(0)
  f = A
  if A == B
    append!(f_list, A)
    return f_list
  end
  if q_fac <= 1 || A > B
    println("ERROR: q_fac <= 1 || A > B")
    return throw(Exception)
  end
  while f < B
    append!(f_list, f)
    f *= q_fac
  end    
  return f_list
end
