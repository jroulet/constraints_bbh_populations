# THIS WILL DELETE ALL THE RUNS!!

read -p "Are you sure? This will delete all the runs" -n 1 -r
echo "\n"
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    [[ "$0" = "$BASH_SOURCE" ]] && exit 1 || return 1
fi

rm -rf GW*/*/ LVT*/*/
rm */parameter_grid
rm */parameter_sets
rm */grid_metadata
