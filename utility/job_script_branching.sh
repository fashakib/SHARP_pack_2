#!/bin/bash
#============================================================
#     JOB SCRIPT TO RUN SHARP PACK v2 ON HPC FOR MULTIPLE 
#     JOBS WITH DIFFERENT K-VALUES TO COMPUTE BRANCHING
#     PROBABLITY FOR TULLY MODEL
#
#     authors    - D.K. Limbu & F.A. Shakib     
#
#     Method Development and Materials Simulation Laboratory
#     New Jersey Institute of Technology
#
#     USAGE:: bash job_script_branching.sh
#
#============================================================

# Path to your executable
#============================================================
root2bin=~/Softwares/SHARP_pack2/bin
exe=${root2bin}/sharp.x

# Set desired range and step
k_min=5  # starting value
k_max=35  # ending value
dk=1.0    # increment

pad=2

base_folder="run_p"
kval=$(printf "%.${pad}f" "$k_min")
folder_name="${base_folder}${kval}"

echo $folder_name
# Create base folder if it doesn't exist
if [ ! -d "${folder_name}" ]; then
    echo "Creating base folder: $base_folder"
#    mkdir "$base_folder"
else
    echo "${folder_name} already exists..."
    exit 1
fi

ndir=0
for k in $(seq $k_min $dk $k_max); do
  kval=$(printf "%.${pad}f" "$k")
  folder_name="${base_folder}${kval}"
  mkdir "${folder_name}"

  cp base_param.in ${folder_name}/param.in

  sed -i  s/pval/${kval}/g ${folder_name}/param.in

  ndir=$((ndir + 1))
  echo "Created folder: $folder_name"
done

echo "No. of folder:" $ndir

mkdir -p logs

# ==== FUNCTION TO WRITE SLURM ARRAY SCRIPT ====
write_array_job() {

  cat <<EOF > submit_array.sh
#!/bin/bash
#SBATCH --job-name=array_run
#SBATCH --array=0-$(($ndir - 1))
#SBATCH --output=logs/myoutput_%A_%a.out
#SBATCH --error=logs/jobError_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM-shared
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=end,fail
#SBATCH --time=1:00:00

k_min=$k_min
dk=$dk

kk=\$(echo "\$k_min + \$SLURM_ARRAY_TASK_ID * \$dk" | bc)

dir_name=\$(printf "${base_folder}%.${pad}f" "\$kk")

echo \$dir_name

cd \$dir_name

 $exe

EOF

  chmod +x submit_array.sh
}
 
# ==== FUNCTION TO WRITE ANALYSIS JOB SCRIPT ====
write_combine_data() {
  cat <<EOF > combine_data.sh
#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --output=logs/myoutput%j.out
#SBATCH --error=logs/jobError.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p RM-shared
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=end,fail


output_file="pop_branch_all.out"
> "\$output_file"  # Clear previous output
dfile="pop_branch.out"

for kval in \$(seq $k_min $dk $k_max); do
    dir=\$(printf "${base_folder}%.${pad}f" "\$kval")

    if [ -f "\$dir/\$dfile" ]; then
        tail -n +3 "\$dir/\$dfile" >> "\$output_file"
    else
        echo "Warning: \$dir/\$dfile not found"
    fi
done

echo "Data combined into \$output_file"

EOF

  chmod +x combine_data.sh
}

write_array_job
write_combine_data

# Submit array job
array_jobid=$(sbatch --parsable submit_array.sh)

echo "Submitting post-processing job dependent on array job ${array_jobid}..."

# Submit average job to run only after array job completes
sbatch --dependency=afterok:$array_jobid combine_data.sh

echo "All jobs submitted successfully!"

