# How to run Jupyter on MARS with a GPU

## 1. Submit a Jupyter job

Create a SLURM script, for example:

`jupyter_gpu_plus.slurm`

```bash
#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=project0062
#SBATCH --job-name=JL.GPU+
#SBATCH --partition=gpuplus
#SBATCH --time=08:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-user=robin.carvajal@glasgow.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

############# LOAD ENVIRONMENT #############

module load apps/miniforge
module load apps/java

eval "$(conda shell.bash hook)"
conda activate jupyter

############# WORKING DIRECTORY #############

cd /mnt/data/project0062

############# START JUPYTER #############

echo "Running on node: $(hostname)"

jupyter lab \
    --no-browser \
    --ip=127.0.0.1 \
    --port=55555 \
    --ServerApp.port_retries=0
```

Before submitting the job, make sure the `logs` directory exists:

```bash
mkdir -p logs
```

Then submit the job:

```bash
sbatch jupyter_gpu_plus.slurm
```

This will return a job ID, for example:

```text
Submitted batch job 123456
```

---

## 2. Find the node running Jupyter

Once the job starts, check which compute node Slurm assigned:

```bash
squeue -j 123456
```

or:

```bash
squeue -j 123456 -o "%.18i %.9T %.20R"
```

You may see something like:

```text
JOBID      STATE     NODELIST(REASON)
123456     RUNNING   gpu101
```

In this example, the allocated node is:

```text
gpu101
```

Do not assume it will always be `gpu101`.

You can also check the job output:

```bash
cat logs/JL.GPU+_123456.out
```

The script prints the hostname when Jupyter starts.

---

## 3. Create an SSH tunnel from your laptop

On your **local laptop**, open a terminal and run:

```bash
ssh -N -L 55556:localhost:55555 gpu101
```

Replace `gpu101` with the node actually assigned to your job.

The ports mean:

```text
Laptop                    HPC compute node
127.0.0.1:55556  ──────►  localhost:55555
```

Keep this terminal open while using Jupyter.

---

## 4. Open Jupyter

On your laptop, open:

```text
http://127.0.0.1:55556/lab
```

If Jupyter asks for authentication, use either your configured Jupyter password or the token printed when the Jupyter server started.

You can find the token in the SLURM output file:

```bash
cat logs/JL.GPU+_<JOBID>.out
```

Look for a URL similar to:

```text
http://127.0.0.1:55555/lab?token=...
```

Do **not** use that URL directly because port `55555` exists on the compute node.

Instead, copy the token and use:

```text
http://127.0.0.1:55556/lab?token=...
```

---

## 5. Stop the Jupyter session

When finished, cancel the Slurm job:

```bash
scancel 123456
```

You can check whether it is still running with:

```bash
squeue -u $USER
```

Closing the browser or SSH tunnel alone does not necessarily stop the Jupyter job.


Use this to connect 

```bash 

ssh -N -J mars -L 55556:localhost:55555 gpu102

```