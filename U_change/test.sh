for a in 2 4 6 8 10 12 14 16 18 
do
cp fig.jl Uc_8e2p${a}U$(($a/2))w_20231023
cd Uc_8e2p${a}U$(($a/2))w_20231023
julia fig.jl ${a}
cd ..
done