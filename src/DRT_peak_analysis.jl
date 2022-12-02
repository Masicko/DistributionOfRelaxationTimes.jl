function merge_peaks(DRT::DRT_struct)
  if DRT.control.peak_merge_tol > 0.0
    work_to_be_done = true
    new_loop = false
    while work_to_be_done || new_loop
      previous_peak_R = 0
      previous_peak_tau = 10.0^(-20)
      new_loop = false
      for (i, tau) in enumerate(DRT.peaks_df.tau_c)
        if abs(log10(previous_peak_tau) - log10(tau)) < DRT.control.peak_merge_tol
          
          #
          DRT.peaks_df.tau_c[i-1] = 10^((log10(DRT.peaks_df.tau_c[i-1] * DRT.peaks_df.tau_c[i]))/2)
          DRT.peaks_df.R[i-1] += DRT.peaks_df.R[i]
          DRT.peaks_df.C[i-1] = DRT.peaks_df.tau_c[i-1]/DRT.peaks_df.R[i-1]
          deleterows!(DRT.peaks_df, i)
          #
          
          new_loop = true
        else
          previous_peak_tau = tau
          previous_peak_R = DRT.peaks_df.R[i]
        end
      end
      work_to_be_done = false
    end
  end
end


function evaluate_RC_peaks_from_DRT!(DRT::DRT_struct)
  
  h_max = maximum(DRT.h)
  threshold = h_max/50.0
  
  function peak_assesement(temp_peak_df)
    R = 0.0
    tau_c = 0.0
    for i in 1:size(temp_peak_df, 1)
      R += temp_peak_df.h[i]
      # only an intermediate step ... weighted average of peak points
      tau_c += temp_peak_df.h[i] * temp_peak_df.tau[i]
    end
    tau_c = tau_c/R
    C = tau_c/R
    return tau_c, R, C
  end
  
  #peak search
  peaks = DataFrame(tau_c = [], R = [], C = [])
  
  temp_peak_df = DataFrame(tau = [], h = [])
  in_peak_bool = false
  for (i, tau) in enumerate(DRT.tau_list)
    if DRT.h[i] > threshold
      push!(temp_peak_df, (tau, DRT.h[i]))
    else
      if size(temp_peak_df, 1) > 0
        push!(peaks, peak_assesement(temp_peak_df))
        
        temp_peak_df = DataFrame(tau = [], h = [])
      end
    end
  end
  
  DRT.peaks_df = peaks
  
  merge_peaks(DRT)
  
  return
end
