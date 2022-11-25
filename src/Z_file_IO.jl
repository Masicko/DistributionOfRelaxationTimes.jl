
function read_Z_file(f_name)    
    df = DataFrame(f = Float32[], Z = Complex[])
    line_is_valid=false
    ###
    imaginary_compensation = 1.0
    ###
    open(f_name) do file
        for ln in eachline(file)
            if line_is_valid
                splitted_row = [parse(Float32,el) for el in split(ln)]
                #push!(df, (splitted_row[1], splitted_row[5], splitted_row[6] ))
                push!(df, (splitted_row[1], splitted_row[5] +  splitted_row[6]*im*imaginary_compensation ))
            else
                if ln=="End Comments"
                    line_is_valid=true
                end
            end
        end
    end
    if length(df.f) > 1 && (df.f[1] > df.f[end])
      df.f = reverse(df.f)
      df.Z = reverse(df.Z)    
    end
    return df 
    #return df.f, df.Z
end
