for a in 4 10 20
do
cp func_lanczos2.jl func_makehamil2.jl fig.jl opcon.jl my_${a}U
cd my_${a}U
julia opcon.jl ${a}
julia fig.jl ${a}
cd ..
done