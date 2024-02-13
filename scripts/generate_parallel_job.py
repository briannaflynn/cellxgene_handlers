import argparse

def generate_slurm_script(job_name, queue, nodes, runtime, control_file, allocation):
    script_template = f"""#!/bin/bash

#SBATCH -J {job_name}                  # Job name
#SBATCH -o {job_name}.%j.o              # Name of stdout output file (%j expands to jobId)
#SBATCH -e {job_name}.%j.e              # Name of stderr output file
#SBATCH -p {queue}                      # Queue (partition) name
#SBATCH -N {nodes}                      # Total number of nodes requested
#SBATCH -n 1                            # Tasks per node
#SBATCH -t {runtime}                    # Run time (hh:mm:ss)
#SBATCH -A {allocation}                 # Allocation name

# Load launcher
module load launcher

# Configure launcher
EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
PRUN=$TACC_LAUNCHER_DIR/paramrun
CONTROL_FILE={control_file}
export LAUNCHER_JOB_FILE={control_file}
export LAUNCHER_WORKDIR=`pwd`
export LAUNCHER_SCHED=interleaved

# Start launcher
$PRUN $EXECUTABLE $CONTROL_FILE
"""
    return script_template

def main():
    parser = argparse.ArgumentParser(description="Generate a custom SLURM job script.")
    parser.add_argument("--job_name", help="Job name")
    parser.add_argument("--queue", help="Queue name on the system")
    parser.add_argument("--nodes", help="Total number of nodes requested")
    parser.add_argument("--runtime", help="Run time in hh:mm:ss format")
    parser.add_argument("--control_file", help="Path to the control file")
    parser.add_argument("--allocation", help="Allocation name (optional)")

    args = parser.parse_args()

    # Generate the SLURM script
    slurm_script = generate_slurm_script(args.job_name, args.queue, args.nodes, args.runtime, args.control_file, args.allocation)

    # Save the script to a file
    script_filename = f"{args.job_name}_slurm_script.sh"
    with open(script_filename, 'w') as script_file:
        script_file.write(slurm_script)
    
    print(f"SLURM job script '{script_filename}' has been generated.")

if __name__ == "__main__":
    main()
