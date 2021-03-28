simul_study_id=1

nrows=$(expr $(cat Simulation_parameters/Simulation_parameters_list_$simul_study_id.txt | wc -l) - 1)

PFER_thr=100
do_exp="FALSE"

echo ID of simulation study: $simul_study_id

for j in $(seq 1 1 $nrows)
do
echo $j
sed "s/{simul_study_id_input}/${simul_study_id}/g" template_comparisons.sh > run2.sh
sed "s/{do_exp_input}/${do_exp}/g" run2.sh > run3.sh
sed "s/{params_id_input}/${j}/g" run3.sh > run4.sh
sed "s/{PFER_thr_input}/${PFER_thr}/g" run4.sh > run.sh
qsub run.sh
done

rm run2.sh
rm run3.sh
rm run4.sh
