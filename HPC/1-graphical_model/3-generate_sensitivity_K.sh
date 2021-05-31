simul_study_id=1

nrows=$(expr $(cat Simulation_parameters/Simulation_parameters_list_$simul_study_id.txt | wc -l) - 1)

topology=random
do_exp="FALSE"

echo ID of simulation study: $simul_study_id
echo topology: $topology

for j in $(seq 1 1 $nrows)
do
echo $j
sed "s/{simul_study_id_input}/${simul_study_id}/g" template_sensitivity_K.sh > run1.sh
sed "s/{topology_input}/${topology}/g" run1.sh > run2.sh
sed "s/{do_exp_input}/${do_exp}/g" run2.sh > run3.sh
sed "s/{params_id_input}/${j}/g" run3.sh > run.sh
qsub run.sh
done

rm run1.sh
rm run2.sh
rm run3.sh
