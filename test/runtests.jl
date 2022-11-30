using DistributionOfRelaxationTimes, Test

include("data_storage.jl")

# [ ] add real world examples

function R_RC_element_test()

  f_list, Z_list = standard_R_RC_element_data()
  
  DRT_output = get_DRT(f_list, Z_list)
  
  # R-L-RC-RC
  correct_answer = (1.0, 3.0e-6, 2.0, 1.0e-4, 5.0, 3.0e-3)

  return all(
           isapprox.(
             (DRT_output.R_ohm,
             DRT_output.L,
             DRT_output.peaks_df.R[1],
             DRT_output.peaks_df.C[1],
             DRT_output.peaks_df.R[2],
             DRT_output.peaks_df.C[2]
             ),
             correct_answer,
             rtol=1.0e-2
           )
  )
end

function test_extrapol_DRT_vs_answer(input, correct_answer, DRT_control=DRT_control())
  DRT_output = get_DRT(input..., DRT_control)
  return @test all(
        isapprox.(
          (
            DRT_output.R_ohm,
            #DRT_output.L,
            DRT_output.peaks_df.R[1],
            DRT_output.peaks_df.C[1],
            DRT_output.peaks_df.R[2],
            DRT_output.peaks_df.C[2]
          ),
          correct_answer,
          rtol=1.0e-2
        )
  )
end

function tau_list_tests()
  # R -  RC - RC
  correct_answer = (20, 20, 0.0005, 20, 0.025)
  

  
  test_extrapol_DRT_vs_answer(extrapol_all_data(), correct_answer)
  test_extrapol_DRT_vs_answer(extrapol_almost_peak(), correct_answer, 
      DRT_control(tau_min_fac=10, tau_max_fac=1, tau_range_fac=10)
  )
  test_extrapol_DRT_vs_answer(extrapol_hard_case(), correct_answer, 
      DRT_control(tau_min_fac=10, tau_max_fac=2, tau_range_fac=100)
  )
end

function frequency_extrapolation_tests()
  correct_answer = (20, 20, 0.0005, 20, 0.025)
  
  test_extrapol_DRT_vs_answer(extrapol_almost_peak(), correct_answer)
  test_extrapol_DRT_vs_answer(extrapol_hard_case(), correct_answer)
end




function run_all_tests()
  @testset "General_tests" begin
      # Check correct answer to R-L-RC-RC element
      @test R_RC_element_test()    
  end

  @testset "Tau_list_tests" begin
    tau_list_tests()
  end
  
  @testset "Frequency_extrapolation" begin
    frequency_extrapolation_tests()
  end
  
end

run_all_tests()
