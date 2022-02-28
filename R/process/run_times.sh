methods=(aucell wmean wsum ulm mlm viper gsva ora gsea udt mdt);
path_raw=data/raw
time=0
max_time=300


echo "lang" "mat" "method" "time" "memr";
 For Python
for method in ${methods[*]}; do
  for i_mat in $(seq45); do
    if [ "$time" -lt "$max_time" ]; then
      test="$(/usr/bin/time --verbose python R/process/run_single_method.py \
      $path_raw/${i_mat}_mat.csv \
      $path_raw/net.csv \
      $method 2>&1 | \
      awk '{if(NR==5){print $8};if(NR==10){print $6}}')";
      time="$(echo $test | awk '{print $1}')";
      memr="$(echo $test | awk '{print $2}')";
      has_h="$(echo $time | grep '\(.*:\)\{2\}')"
      if [[ $has_h ]]; then
        hours="$(echo $time | awk -F: '{print $1 * 3600}')"
        minutes="$(echo $time | awk -F: '{print $2 * 60}')"
        seconds="$(echo $time | awk -F: '{print $3}')"
        time="$(echo $hours $minutes $seconds | awk '{print $1+$2+$3}')"
      else
        minutes="$(echo $time | awk -F: '{print $1 * 60}')"
        seconds="$(echo $time | awk -F: '{print int( $2 )}')"
        time="$(echo $minutes $seconds | awk '{print $1+$2}')"
      fi
      echo "Python" $i_mat $method $time $memr;
    else
      echo "Python" $i_mat $method "NA" "NA";
    fi
  done
  time=0;
done

# For R
for method in ${methods[*]}; do
  for i_mat in $(seq 4); do
    if [ "$time" -lt "$max_time" ]; then
      test="$(/usr/bin/time --verbose Rscript R/process/run_single_method.R \
      $path_raw/${i_mat}_mat.csv \
      $path_raw/net.csv \
      $method 2>&1 | \
      awk '{if(NR==5){print $8};if(NR==10){print $6}}')";
      time="$(echo $test | awk '{print $1}')";
      memr="$(echo $test | awk '{print $2}')";
      has_h="$(echo $time | grep '\(.*:\)\{2\}')"
      if [[ $has_h ]]; then
        hours="$(echo $time | awk -F: '{print $1 * 3600}')"
        minutes="$(echo $time | awk -F: '{print $2 * 60}')"
        seconds="$(echo $time | awk -F: '{print $3}')"
        time="$(echo $hours $minutes $seconds | awk '{print $1+$2+$3}')"
      else
        minutes="$(echo $time | awk -F: '{print $1 * 60}')"
        seconds="$(echo $time | awk -F: '{print int( $2 )}')"
        time="$(echo $minutes $seconds | awk '{print $1+$2}')"
      fi
      echo "R" $i_mat $method $time $memr;
    else
      echo "R" $i_mat $method "NA" "NA";
    fi
  done
  time=0;
done
