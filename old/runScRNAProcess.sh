#$ -N scProcess
#$ -l week
#$ -pe smp 8
#$ -cwd
#$ -e /projects/PPC/analysis/ppc_pilot/logs/log_process_scRNA.err
#$ -o /projects/PPC/analysis/ppc_pilot/logs/log_process_scRNA.out

/home/matteo/software/R-3.5.1/bin/Rscript /frazer01/home/mdonovan/pancreas_scRNA_map/process_scRNAseq.R
