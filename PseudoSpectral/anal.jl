using PyPlot, Random, DelimitedFiles, CSV, DataFrames

function MSD(xi,yi)
    aux = 0
    MSD = zeros(0)
    push!(MSD,0)
    for i in 2:length(xi)
        aux+= ( (yi[i]-yi[1])*(yi[i]-yi[1]) + (xi[i]-xi[1])*(xi[i]-xi[1]) )
        push!(MSD,aux/i)
    end
    return MSD
end
function main()

end
main()