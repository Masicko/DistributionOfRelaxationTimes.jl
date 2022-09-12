module Example_02_lambda

using DistributionOfRelaxationTimes
using PyPlot

function main()
  f_range = [0.10, 0.21, 0.43, 0.89, 1.83, 3.79, 7.85, 16.24, 33.60, 69.52, 143.84, 297.64, 615.85, 1274.27, 2636.65, 5455.59, 11288.38, 23357.21, 48329.30, 100000.00]
  
  Z_range = [5919.90 - 15.79im, 5919.58 - 32.68im, 5918.18 - 67.58im, 5912.24 - 139.49im, 5887.12 - 285.74im, 5785.04 - 566.88im, 5428.94 - 997.19im, 4640.21 - 1257.83im, 3871.84 - 978.97im, 3537.68 - 564.96im, 3442.94 - 315.40im, 3418.14 - 219.69im, 3405.51 - 242.57im, 3373.90 - 396.07im, 3249.67 - 742.03im, 2808.42 - 1305.92im, 1779.41 - 1698.97im, 701.96 - 1361.47im, 208.29 - 777.65im, 65.93 - 392.51im]
  
  f_nyq = figure(1)
  f_drt = figure(2)
  
  for lambda in collect(0.0 : 0.5 : 1.0)
    my_DRT = get_DRT(f_range, Z_range, DRT_control(lambda=lambda))
    
    figure(1)
    nyquistPlot(my_DRT.EIS_df.Z, label="l = $(lambda)")  
    
    figure(2)
    plot_DRT_h(my_DRT, false, label="l = $(lambda)")
  end
  
  
  figure(1)
  nyquistPlot(Z_range, label="data")
  
  
  return
end

end # module
