XENBASE GEO PIPELINE

## OVERVIEW
This pipeline takes a "run file" containing SRRs and metadata and uses it to produce the outputs needed for the Xenbase server. It can simultaneously handle RNA-seq, ChIP-seq, and ATAC-seq data for v10 tropicalis and laevis across multiple GSEs. The analysis (RSEM, Bowtie2, RUVSeq, edgeR, MACS2, etc.) remains unchanged compared to the original CSBB pipeline for consistency and compatibility.

## HOW TO USE
- Download and unzip the reference files into a "data" directory. Files are available here: https://drive.google.com/file/d/1gtfdp3bDVWutqaeUxj-GbMhUp7hjiX2K- /view?usp=sharing
- Add the paths for the run file and the output directory to pipeline.sh, then submit the job.
- The only rule for the run file are to include all lines with the same combination of GSE, species, and assay. Or for simplicity, just include all lines for a given GSE. This ensures all the data needed is present and will be grouped together in the final output. 

## IMPORTANT NOTES
- All results needed will be copied to the specified output directory. Once that has been checked the work directory can be safely deleted.
- If the pipeline fails check the logs/report as well as the work directory of any processes that failed. The most likely issues are either memory/space related or a problem with the run file.
- After a failure you can resume the pipeline part way through using the "RESUME" option in pipeline.sh. 
- The "nextflow.config" file controls how jobs are submitted for each process. Default values have been chosen that should optimize the speed and amount of jobs that can be run at once, but the file can be modified as needed to provide additional memory or CPUs. 
- This pipeline performs automatic cleanup of files when they are no longer needed. They may still appear in the work directories but are replaced by sparse files that take up no space. This workaround "tricks" Nextflow into thinking the cached results are still there so the resume function can still be used.
