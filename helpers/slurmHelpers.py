import os

def addSlurmTask(jobName, 
                 logName, 
                 chromosome, 
                 coordinate, 
                 geneName, 
                 fastaPath, 
                 tsvPath, 
                 partition = "AI", 
                 folder = '~/'
                 ):
    print(f'#!/bin/bash\
    #SBATCH --nodes=1\
    #SBATCH --ntasks=1\
    #SBATCH --job-name={jobName}\
    #SBATCH --output=logs/{logName}.out\
    #SBATCH --error=logs/{logName}.err\
    #SBATCH --partition={partition}\
    source activate primers\
    srun cd {folder} && srun python3 singleTask.py -chr {chromosome} -coord {coordinate} -gn {geneName} -fasta {fastaPath} -tsv {tsvPath}')