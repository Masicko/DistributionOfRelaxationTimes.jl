using DistributionOfRelaxationTimes, Test

include("data/data_storage.jl")

# [ ] add real world examples
# [ ] add HF example

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
      DRT_control(tau_min_fac=10, tau_max_fac=1, tau_sampling_fac=10)
  )
  test_extrapol_DRT_vs_answer(extrapol_hard_case(), correct_answer, 
      DRT_control(tau_min_fac=10, tau_max_fac=2, tau_sampling_fac=100)
  )
  test_extrapol_DRT_vs_answer(extrapol_hard_case(), correct_answer, 
      DRT_control(tau_min_abs=1.0e-4, tau_max_abs=1.0e4, tau_sampling_fac=100)
  )
  
  
  # TODO peak merge test, lambda
end

function plotting_tests()
  @test begin
    DRT_output = get_DRT(extrapol_all_data()...)
  
    plot_Nyquist(DRT_output.EIS_df.Z)
    plot_DRT_h(DRT_output)
    plot_DRT_RC(DRT_output, true, false)
    plot_DRT_Rtau(DRT_output)
    plot_DRT_Rf(DRT_output)
    true
  end
end

modname(fname)=splitext(basename(fname))[1]

function run_tests_from_directory(testdir,prefix)
    println("Directory $(testdir):")
    examples=modname.(readdir(testdir))
    for example in examples
        if length(example)>=length(prefix) &&example[1:length(prefix)]==prefix
            println("  $(example):")
            path=joinpath(testdir,"$(example).jl")
            @eval begin
                include($path)
                # Compile + run test
                @test eval(Meta.parse("$($example).test()"))
                # Second run: pure execution time.
                @time eval(Meta.parse("$($example).test()"))
            end
        end
    end
end

function run_all_tests()
  @testset "General_tests" begin
      # Check correct answer to R-L-RC-RC element
      @test R_RC_element_test()    
  end

  @testset "Tau_list_tests" begin
    tau_list_tests()
  end
  
  @testset "Plotting" begin
    plotting_tests()
  end
  
  @testset "Examples" begin
    run_tests_from_directory("../examples", "Example")
  end
end

run_all_tests()
