# default number of figures
const DRT_standard_figure = 33

function plot_Nyquist(Z; fignum=1, label="")
    s = subplot(111)
    title("Nyquist plot")
    xlabel("Re\$(Z) \\ [\\Omega]\$")
    ylabel("-Im\$(Z) \\ [\\Omega]\$")
    PyPlot.plot(real(Z), -imag(Z), label=label)
    if label!=""
      legend()
    end
    grid(true)
    s.set_aspect(1.0)  
end


function plot_DRT_h(DRT::DRT_struct, to_standard_figure=true, print_bool=false, plot_lambda=true; label="")
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  
  title("DRT")
  plot(log10.(DRT.tau_range), DRT.h, "-x", label=(label=="" ? "\$\\lambda\$=$(DRT.control.lambda)" : label))
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("\$h(\\tau)\$ [Ohm]")
  if (to_standard_figure || plot_lambda)
    legend()
  end
  if print_bool
    println("non-DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)")  
  end
end
  
function plot_DRT_RC(DRT::DRT_struct, to_standard_figure=true, print_bool=true)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
#   
  title("RC diagram")
  xlabel("C [F]")
  ylabel("R [Ohm]")
  plot(peaks_df.C, peaks_df.R, "x")
  grid(true)
  
  if print_bool
    println("non-DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)") 
    for i in 1:size(peaks_df,1)
      println(">> f$i = $(1/(2*pi*peaks_df.C[i]*peaks_df.R[i]))    R$i = $(peaks_df.R[i])   C$i = $(peaks_df.C[i])   ... (tau_c$i = $(peaks_df.tau_c[i]))")
    end
  end
end

function plot_DRT_Rtau(DRT::DRT_struct, to_standard_figure=true, print_bool=false)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
  
  title("Rtau diagram")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("R [Ohm]")
  xlim(log10.([DRT.tau_range[1], DRT.tau_range[end]])...)
  
  
  previos_maximum = (PyPlot.gca()).get_ylim()[2]
  
  ylim(0, max(
    (length(peaks_df.R) > 0 ? maximum(peaks_df.R)*1.1 : 0.0001), 
    previos_maximum)  )
  plot(log10.(peaks_df.C.*peaks_df.R), peaks_df.R, "o")
  grid(true)

end

function plot_DRT_Rf(DRT::DRT_struct, to_standard_figure=true, print_bool=false)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
  
  title("Rf diagram")
  xlabel("log10(\$f\$ [Hz])")
  ylabel("R [Ohm]")
  xlim(log10.(1 ./(2*pi*[DRT.tau_range[end], DRT.tau_range[1]]))...)
  plot(log10.(1 ./(2*pi*peaks_df.C.*peaks_df.R)), peaks_df.R, "o")
  grid(true)
  
  if print_bool
    for i in 1:size(peaks_df,1)
      println(">> f$i = $(1/(2*pi*peaks_df.C[i]*peaks_df.R[i]))    R$i = $(peaks_df.R[i])    ... (C$i = $(peaks_df.C[i]))")
    end
  end
end
