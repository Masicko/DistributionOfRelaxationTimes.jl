
"""

Costruction of both matrix A and the right hand side b.
"""

function construct_A_matrix_and_b(f_nodes, tau_range, Z, lambda)
  N_f = size(f_nodes, 1)
  N_tau = size(tau_range, 1)
  
  # h, R_ohm, L
  n_cols = N_tau + 2  
  # real(Z), imag(Z), regularization
  if lambda == 0.0
    n_rows = 2*N_f
  else
    n_rows = 2*N_f + n_cols
  end
#   @show tau_range
#   @show f_nodes
    
  A = Matrix{Float64}(undef, n_rows, n_cols)
  b = Vector{Float64}(undef, n_rows)
  
  # assemble A and b
  for (i, f) in enumerate(f_nodes)
    #RC
    for (j, tau) in enumerate(tau_range)
      A[i, j]       = real(1/(1 + im*(2*pi*f)*tau))
      A[N_f + i, j] = imag(1/(1 + im*(2*pi*f)*tau))
    end
    # R_ohm
    A[i, N_tau + 1] = 1
    A[N_f + i, N_tau + 1] = 0
    # L
    A[i, N_tau + 2] = 0
    A[N_f + i, N_tau + 2] = 2*pi*f
    
    # b
    if Z != nothing
      b[i] = real(Z[i])
      b[N_f + i] = imag(Z[i])
    end
  end
  
  if lambda != 0.0
    # assemble "regularization" part of A and b
    A[2*N_f + 1 : end, :] .= Diagonal([lambda for i in 1:n_cols])
    b[2*N_f + 1 : end] .= 0
  end

  return (A, b, N_f, N_tau)
end

"""

When DRT has nonzero the lowest frequency (the highest tau) value, 
this value is used as an estiamte for R_ohm resistance 
and optimization is run once more.
"""

function high_frequency_adjustment!(solution, A, b, N_tau)
  if solution[1] > 0.000001
    println(" !  High frequency adjustment ---------------------- ")
    R_ohm_estimate = solution[N_tau + 1] + solution[1]
  
    next_row = zeros(N_tau + 2)
    next_row[N_tau + 1] = 1
    A_with_R_ohm = [A; transpose(next_row)]
    b_with_R_ohm = [b; R_ohm_estimate]
    
    solution .= nonneg_lsq(A_with_R_ohm, b_with_R_ohm; alg=:nnls)
  end
  return
end

"""

Optimization procedure sometimes returns high-frequency (and low-frequency) peaks which
do not make physical sense. These artifacts are removed from the DRT solution.
"""

function delete_corner_effects!(solution, control)
  if control.lambda == 0.0 && false
    # solution[end-1] is ohmic resistance R_ohm 
    # solution[end] is inductance L
    solution[end-7 : end-2] .= 0
    solution[1 : 3] .=0
  end
end

"""

The impedance data are generated using the extracted DRT information. 
It is useful in order to compare DRT approximation to the original impedance data.
"""

function reconstruct_EIS_from_DRT(A, N_f, solution, control, tau_range, EIS_df)
  if control.f_range != nothing  
    f_nodes = geomspace_by_fac(control.f_range...)
        
    # assemble A_new w.r.t. specified f_range
    (A_new, b, N_f, N_tau) = construct_A_matrix_and_b(f_nodes, tau_range, nothing, 0.0)
  
    b_new = A_new*solution
    
    Z_new = Array{Complex}(undef, N_f)
    for i in 1:N_f
      Z_new[i] = Complex(b_new[i], b_new[N_f + i])
    end
    
    EIS_new = DataFrame(f = f_nodes, Z = Z_new)
  else      
    b_new = A*solution
    
    EIS_new = deepcopy(EIS_df)
    for i in 1:N_f
      EIS_new.Z[i] = b_new[i] + im*b_new[N_f + i]
    end
  end
  return EIS_new
end

# function capture_corner_effects()
#   
# end


"""

Wrapper which do not use DataFrames for the input impedance data.
"""

function get_DRT(f_range, Z_range, control::DRT_control=DRT_control())
  EIS_df = DataFrame(f = f_range, Z = Z_range)
  return get_DRT(EIS_df::DataFrame, control)
end

"""

The function computes DRT spectrum from impedance data and returns the DRT struct with the DRT function h, reconstructed impedance data and dataframe with peaks.
"""

function get_DRT(EIS_df::DataFrame, control::DRT_control=DRT_control())  
  tau_min = 1.0/(2*pi*EIS_df.f[end]) / control.tau_min_fac
  tau_max = 1.0/(2*pi*EIS_df.f[1]) * control.tau_max_fac
  tau_range = geomspace(tau_min, tau_max, N=control.tau_range_fac*size(EIS_df.f, 1)-2)
  
  (A, b, N_f, N_tau) = construct_A_matrix_and_b(EIS_df.f, tau_range, EIS_df.Z, control.lambda)
  
  
  solution= nonneg_lsq(A, b; alg=:nnls)
  
  # processing  
  
  ## TODO ! ! ! !
  ## 
  ## enlarge tau domain to capture corner effects
  ## 
  #capture_corner_effects(solution)
  
  high_frequency_adjustment!(solution, A, b, N_tau)
  
  
  delete_corner_effects!(solution, control)
  
  # reconstruction of EIS from DRT
  EIS_new = reconstruct_EIS_from_DRT(A, N_f, solution, control, tau_range, EIS_df)
  
  DRT_out = DRT_struct(EIS_new, tau_range, solution[1:end - 2], solution[end-1], solution[end], DataFrame(), control)
  evaluate_RC_peaks_from_DRT!(DRT_out)
  
  return DRT_out
end
