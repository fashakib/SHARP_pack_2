#!/bin/bash
#============================================================
#     JOB SCRIPT TO RUN SHARP PACK v2 ON LOCAL MACHINE
#     WHILE RUNNING PARALLEL, BASED ON n-CPU, IT CREATES
#     n-DIRECTORIES AND RUN JOBS PARALLELLY, AND FINALLY 
#     COMPUTES AVERAGE OF THE DISIRED POPULATATION (NEED TO
#     SPECIFY FOR INPUT, SEE EXAMPLE IN, WRITE_AVERAGE_INPUT)
#
#     authors    - D.K. Limbu & F.A. Shakib     
#
#     Method Development and Materials Simulation Laboratory
#     New Jersey Institute of Technology
#
#     USAGE:: bash job_script_hpc.sh
#
#============================================================

# Path to your executable
#============================================================
root2bin=~/Softwares/SHARP_pack2/bin
exe=${root2bin}/sharp.x
average_exe=${root2bin}/average.x

# ================================
# CPU configuration
# ================================
maxcore=1
ncore=$(awk '/ncore/{print $2}' param.in)
maxcore=$((ncore > maxcore ? ncore : maxcore))
echo "cpu: $maxcore"

root_dir="run"
param_file="param.in"
ncore=$maxcore

# Determine padding width
if [ "$ncore" -lt 10 ]; then
   pad=1
elif [ "$ncore" -lt 100 ]; then
   pad=2
else
   pad=3
fi

# ================================
# FUNCTIONS
# ================================
# ==== FUNCTION TO CREATE AVERAGE INPUT ====
write_average_input() {
  #fname :: File name of specific population data for averaging
  #df    :: No. of data column(s)(except first column)?

  fname="pop_diabat3.out"
  df=2

  echo "Calculating average results from $ncore parallel jobs."

  rm -f input

  echo ${root_dir} >> input
  echo ${fname} >> input
  echo ${ncore} >> input
  echo ${df} >> input
  echo ${pad} >> input

}

# ==== FUNCTION TO CREATE PARALLEL DIRECTORIES ====
paralleldir() {
  # Check if the first directory already exists
  first_dir=$(printf "${root_dir}%0${pad}g" 1)
  if [ -d "$first_dir" ]; then
    echo '**********************************'
    echo "$first_dir exists. Remove it first!"
    echo '**********************************'
    exit 1
  fi

  # Create directories and copy files
  for i in $(seq -f "%0${pad}g" 1 $ncore); do
    dir="${root_dir}${i}"
    mkdir "$dir"
    cp param.in "$dir/param.in"
    echo "$i" > "$dir/icpu.in"
  done

  echo "Created $ncore directories: ${root_dir}$(printf "%0${pad}d" 1) to ${root_dir}${ncore}"
}

# ==== FUNCTION TO RUN SERIAL JOB  ====
run_serial() {
  echo "Running serial job on 1 core..."
  $exe &
}

# ==== FUNCTION TO RUN n-PARALLEL JOBS  ====
run_parallel() {
  echo "Running $ncore jobs in parallel..."

  paralleldir

  # Launch jobs in background
  for i in $(seq -f "%0${pad}g" 1 $ncore); do
    dir="${root_dir}${i}"
    (
      cd "$dir" || exit
      $exe
    ) &
  done

  # Wait for all jobs to complete
  wait
  echo "All parallel jobs completed."

  # Prepare averaging input and run average job
  write_average_input
  echo "Running averaging step..."
  $average_exe
}

# ================================
# MAIN
# ================================
if [ $ncore -eq 1 ]; then
  run_serial
else
  run_parallel
fi

echo "All jobs finished successfully."

