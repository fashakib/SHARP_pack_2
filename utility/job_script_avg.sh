#!/bin/bash
#============================================================
#     JOB SCRIPT TO COMPUTE AVERAGE OF THE DISIRED POPULATATION 
#     (NEED TO SPECIFY FOR INPUT, SEE EXAMPLE IN, 
#     WRITE_AVERAGE_INPUT), AFTER PARALLEL RUN COMPLETE
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

# ==== FUNCTION TO CREATE AVERAGE INPUT ====
write_average_input() {
  #fname :: File name of specific population data for averaging
  #df    :: No. of data column(s)(except first column)?

  fname="pop_diabat2.out"
  df=5

  echo "Calculating average results from $ncore parallel jobs."

  rm -f input

  echo ${root_dir} >> input
  echo ${fname} >> input
  echo ${ncore} >> input
  echo ${df} >> input
  echo ${pad} >> input

}

# ================================
# CPU configuration
# ================================
maxcore=1
ncore=$(awk '/ncore/{print $2}' param.in)
maxcore=$((ncore > maxcore ? ncore : maxcore))
echo "ncores(s): $maxcore"

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
# MAIN
# ================================
if [ $ncore -eq 1 ]; then
   echo "ncore = 1,  serial job. No average necessary!"
   exit 1
else
  # Prepare averaging input and run average job
  write_average_input
  echo "Running averaging step..."
  $average_exe

fi

