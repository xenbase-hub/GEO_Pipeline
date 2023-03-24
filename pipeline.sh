# TO RUN: ./pipeline.sh | bsub;
# TO RESUME: Enter any value after "RESUME="
# TO MONITOR SPACE: du -sch /data/xenbase
# TO MONITOR PROGRESS: Use bpeek on the "xenbase" job

RUNFILE=/data/xenbase/runfile_2023-03-21.txt
OUTDIR=/data/xenbase/Processed/
RESUME=

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 100:00
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -M 8000
#BSUB -e /data/xenbase/logs/%J.err
#BSUB -o /data/xenbase/logs/%J.out
#BSUB -J xenbase

module load nextflow
nextflow -C /data/xenbase/NF_Pipeline/nextflow.config \
run /data/xenbase/NF_Pipeline/main.nf \
-with-report report.html \
--runFile $RUNFILE \
--outDir $OUTDIR ${RESUME:+"-resume"}
EOF
