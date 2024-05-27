data_dir="../data_for_scripts/recombination_detection/data_for_filtration"
rdp="RDP_pre_filtered_fixed.tsv"

echo "Remove weak dichotomy"
for i in $(seq 1 3); do
  python3 filter_tree.py -t $data_dir/ML_Dom${i}.raxml.support -s 70 -o $data_dir/ml_dom${i}.70.raxml.support
done
echo "Removing done"

echo "Filtering trees"
for i in $(seq 1 3); do
  python3 recomb2tree.py -r $data_dir/$rdp -t $data_dir/ML_Dom${i}.raxml.support -s 70 -d ${i} -o $data_dir/dom${i}.filtered.recomb.raxml
  sed 's/,:0//g' $data_dir/dom${i}.filtered.recomb.raxml > $data_dir/dom${i}.filtered.recomb.zero_fixed.raxml
done
echo "Filtering trees done"

echo "Unpivoting table"
python3 unpivot_table.py -r $data_dir/$rdp -o $data_dir/rdp.unpivot.csv
echo "Unpivoting done"

echo "Filtering parents"
trees_args=""
for i in $(seq 1 3); do
  trees_args="$trees_args -t $data_dir/dom${i}.filtered.recomb.zero_fixed.raxml"
done

python3 filter_parents.py -r $data_dir/rdp.unpivot.csv $trees_args -s 70 -l 3 -o $data_dir/recomb.unpivot.filtered.s70.l3.csv
echo "Filtering done"

echo "Filtering events"
python3 filter_events.py -i $data_dir/recomb.unpivot.filtered.s70.l3.csv \
  -o $data_dir/recomb.unpivot.filtered_events.s70.l3.csv > $data_dir/filtering.log
echo "Filtering done"

echo "Merging tables"
python3 merge_filtered_with_original.py -i $data_dir/$rdp -f $data_dir/recomb.unpivot.filtered_events.s70.l3.csv \
  -o $data_dir/merged.filtered_events.s70.l3.csv
echo "Merging done"
