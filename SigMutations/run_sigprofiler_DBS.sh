#BSUB -o Myjob.%J.log
#BSUB -e Myjob.%J.err
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -W 48:00

sample=SAMPLE

module load singularity/3.0.1-to-ve-removed

singularity run --bind  ${PWD}:/data ~/lingqi_workspace/sigprofiler-py-cpu.sif python3 /app/SigProfilerHelper/run_sigprofiler.py -g GRCh37 -c 4 -i 100 -s 1 -e 60 -t /data/${sample}/output/DBS/${sample}.DBS78.all -d /data/${sample}/output
