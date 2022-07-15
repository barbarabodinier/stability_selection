PFER_thr=5

for j in $(seq 1 1 2)
do
echo ID of simulation study: $j
sed "s/{simul_study_id_input}/${j}/g" template_comparisons.sh > run1.sh
sed "s/{PFER_thr_input}/${PFER_thr}/g" run1.sh > run.sh
qsub run.sh
done

rm run1.sh
