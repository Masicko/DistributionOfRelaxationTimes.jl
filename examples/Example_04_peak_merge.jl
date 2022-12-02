module Example_04_peak_merge

using DistributionOfRelaxationTimes
using PyPlot

function main()
  f_list = Any[0.01, 0.01412537544622754, 0.0199526231496888, 0.028183829312644536, 0.039810717055349734, 0.056234132519034905, 0.07943282347242814, 0.11220184543019636, 0.15848931924611132, 0.22387211385683395, 0.31622776601683794, 0.44668359215096315, 0.6309573444801932, 0.8912509381337456, 1.2589254117941673, 1.7782794100389228, 2.51188643150958, 3.548133892335755, 5.011872336272722, 7.079457843841379, 10.0, 14.12537544622754, 19.952623149688797, 28.183829312644534, 39.810717055349734, 56.23413251903491, 79.43282347242814, 112.2018454301963, 158.48931924611142, 223.872113856834, 316.2277660168379, 446.683592150963, 630.957344480193, 891.2509381337459, 1258.9254117941675, 1778.2794100389228, 2511.88643150958, 3548.133892335753, 5011.872336272725, 7079.457843841381, 10000.0]
  Z_list = Any[26.933719711876503 - 0.8178533448733899im, 26.868620000060016 - 1.1477379281106466im, 26.741245439858034 - 1.6004652095869243im, 26.496677185508922 - 2.204418844753519im, 26.043640530917536 - 2.9665162184033242im, 25.257632720147996 - 3.8293083246261506im, 24.037219574286116 - 4.617295289512633im, 22.434472779794664 - 5.053442597590856im, 20.73607665054383 - 4.940132735468191im, 19.300125929699494 - 4.354109118611884im, 18.29951460119249 - 3.570821956331194im, 17.6908283111082 - 2.8373108907033844im, 17.347116739540677 - 2.2764581794058576im, 17.15389438632848 - 1.921571053672256im, 17.03141958001423 - 1.7708053266126194im, 16.923702776769506 - 1.818925097610685im, 16.780324262531998 - 2.06875245315941im, 16.537748522301193 - 2.5286245978652326im, 16.10243766866905 - 3.19407859542437im, 15.34853553277152 - 4.006332312426719im, 14.165743263555786 - 4.796646468979229im, 12.583928555372625 - 5.290404752950778im, 10.86779297909743 - 5.2733147534589895im, 9.376119853612918 - 4.799213741286676im, 8.295206826688577 - 4.1418425024904835im, 7.579294787371817 - 3.5652005655471193im, 7.073326135705628 - 3.205488683448504im, 6.613502416914236 - 3.085109175349607im, 6.064799772007791 - 3.1456134797130244im, 5.348847896778272 - 3.26609458257111im, 4.477833739155986 - 3.3024442198978954im, 3.55028944167983 - 3.153736916352337im, 2.6994399050171864 - 2.803543934191091im, 2.0321743634853906 - 2.318665269046926im, 1.581528388868406 - 1.8061771750811506im, 1.3112823697819886 - 1.3497856802130221im, 1.1615924222871021 - 0.9837064423865067im, 1.0824771216409765 - 0.7069841391996191im, 1.0417223018228337 - 0.5043711478047129im, 1.0210091395537235 - 0.35845956114164496im, 1.0105544408727343 - 0.2542676670704162im]
  
  f_nyq = figure(1)
  f_drt = figure(2)
  
  for peak_merge_tol in [0.0, 0.5]
    my_DRT = get_DRT(f_list, Z_list, DRT_control(peak_merge_tol=peak_merge_tol))
    
    figure(1)
    plot_Nyquist(my_DRT.EIS_df.Z, label="tol = $(peak_merge_tol)")  
    
    figure(2)
    plot_DRT_h(my_DRT, false, label="tol = $(peak_merge_tol)")
    
    figure(3)
    plot_DRT_Rtau(my_DRT, false, label="tol = $(peak_merge_tol)")
  end
  
  
  figure(1)
  plot_Nyquist(Z_list, label="data")
  
    
  return true
end

function test()
  main()
end

end # module
