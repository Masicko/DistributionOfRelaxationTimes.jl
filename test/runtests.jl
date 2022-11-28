using DistributionOfRelaxationTimes, Test

# [ ] standardní case
# [ ] 
# [ ] 


function R_RC_element_test()
  f_list = [0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0, 16384.0, 32768.0,        
  65536.0, 131072.0, 262144.0]
  Z_list = [7.9999566250062575 - 0.014804038716500654im, 7.999826504539605 - 0.029607311013954685im, 7.999306090382245 - 0.05920849147321719im, 7.997225516509134 - 0.11836796402074179im, 7.988920507331971 - 0.2363445912343561im, 7.955974654576409 - 0.46958432924037374im, 7.828430016246774 - 0.9151288224796308im, 7.377764326740351 - 1.6603747779009874im, 6.187578798483693 - 2.423483564393849im, 4.526395789326643 - 2.3428639092377055im, 3.4920252494595356 - 1.5733429341387992im, 3.120896373931658 - 0.9653675926047738im, 2.983681817667524 - 0.7228105806992193im, 2.8210053520607663 - 0.785109509619991im, 2.4166123325165447 - 1.00398585317664im, 1.7535928232882099 - 1.0015305103631473im, 1.2624850542384358 - 0.6624837626405642im, 1.0727779748610742 - 0.3101732217968406im, 1.0187045081983228 - 0.04452524954644659im, 1.0047091282613076 + 0.208681260966955im, 1.0011793628768881 + 0.5675016337277559im, 1.0002949710575257 + 1.2102334554501843im, 1.0000737509150293 + 2.4581021305355573im, 1.000018438238243 + 4.935024388053915im]
  
  
  DRT_output = get_DRT(f_list, Z_list)
  
  # R-L-RC-RC
  prms = (1.0, 3.0e-6, 2.0, 1.0e-4, 5.0, 3.0e-3)

  return all(
           isapprox.(
             (DRT_output.R_ohm,
             DRT_output.L,
             DRT_output.peaks_df.R[1],
             DRT_output.peaks_df.C[1],
             DRT_output.peaks_df.R[2],
             DRT_output.peaks_df.C[2]
             ),
             prms,
             rtol=1.0e-2
           )
  )
end

@testset "General_tests" begin
    # Check correct answer to R-L-RC-RC element
    @test R_RC_element_test()    
    
    
    
#     # Check the circuitevolution function.
#     @test library_fit == circuitevolution(measurements,frequencies,initial_population = library).Circuit
#     # Evaluate properties of the generated circuit encoding.
#     @test isoperation(encoding[1]) == true
#     # Only terminals in the encoding's tail.
#     @test all(isterminal.(collect(encoding[head+1:end]))) == true
#     # conversion of circuit object to user-readable circuit.
#     @test readablecircuit(example_circuit) == "[C1,P2]-R3-[R4,R5]"
#     # Check the output of circuit simulation.
#     @test example_function(exam_params,100) == 1789.4845259396134 - 10.86703970004648im
#     # Replace redundant CPEs and simplify
#     replace_redundant_cpes!(example_circuit_redundant_CPE) # converts the redundant CPE to an equivalent capacitor.
#     simplifycircuit!(example_circuit_redundant_CPE) # simplifies all parralelly or serially connected similar components.
#     # Circuit simplification should lead to a single resistor element R3.
#     @test readablecircuit(example_circuit_redundant_CPE) == "C1-R2"
#     # Parameteroptimisation : basic checks of solution lengths and bounds.
#     optparams = parameteroptimisation("[C1,P2]-R3",measurements,frequencies)
#     @test length(optparams) == 4
#     C1,P2_1,P2_2,R3 = optparams
#     @test 0<C1<10
#     @test 0<P2_1<1.0e9
#     @test 0<P2_2<1
#     @test 0<R3<1.0e9
#     # Subtrees length checking. There ought to be as many circuits as coding terminal elements.
#     @test length(subcircuits(example_circuit)) == 5 

end
