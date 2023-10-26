for a in {4..9}
do
cp fig.jl pre_${a}e2p20U20w20231022
cd pre_${a}e2p20U20w20231022
julia fig.jl ${a}
cd ..
done

#cp job.sh opcon_${a}e2p20U10w20231022
#cd opcon_${a}e2p20U10w20231022
#sed -e "s/#PBS -l/#PBS -l nodes=qcrt${a}.ad12.riken.jp:ppn=28/" ./job.sh > temp.sh
#mv temp.sh job.sh
#sed -e "s/opcon.jl/opcon.jl ${a} 2/" ./job.sh > temp.sh
#mv temp.sh job.sh
#qsub job.sh
#cd ..