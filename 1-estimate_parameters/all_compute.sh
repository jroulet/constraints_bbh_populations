for i in LVT151012 GW150914 GW170104 GW170608 GW170814 GW151226
do
	python compute_likelihood.py $i >> compute_likelihood.log || exit 0
done
