for a in 10 20 30 40 50 60
do
cp fig.jl ge_10e2p${a}U$(($a/2))w_20231025
cd ge_10e2p${a}U$(($a/2))w_20231025
julia fig.jl ${a}
cd ..
done