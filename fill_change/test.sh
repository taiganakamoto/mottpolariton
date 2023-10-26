for a in {4..12}
do
cp fig.jl fc_8s${a}e2p20U10w_20231023
cd fc_8s${a}e2p20U10w_20231023
julia fig.jl ${a}
cd ..
done