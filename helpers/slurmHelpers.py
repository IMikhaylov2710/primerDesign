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
    print(f'sbatch\
    --nodes=1\
    --ntasks=1\
    --job-name={jobName}\
    --output=logs/{logName}.out\
    --error=logs/{logName}.err\
    --partition={partition}\
    source activate primers\
    srun cd {folder} && srun python3 singleTask.py -chr {chromosome} -coord {coordinate} -gn {geneName} -fasta {fastaPath} -tsv {tsvPath}')