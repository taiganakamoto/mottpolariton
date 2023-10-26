for a in {0..8}
do
cd opcon_10e${a}p20U10w_20231019
julia fig.jl 10 ${a}
cd ..
done