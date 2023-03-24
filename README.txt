// XENBASE GEO PIPELINE
// By Konrad Thorner, January 2023

// OVERVIEW
// This pipeline takes a "run file" containing SRRs and metadata and uses it to produce the outputs needed for the Xenbase server.
// It can simultaneously handle RNA-seq, ChIP-seq, and ATAC-seq data for v10 tropicalis and laevis across multiple GSEs. 
// The analysis (RSEM, Bowtie2, RUVSeq, edgeR, MACS2, etc.) remains unchanged compared to the original CSBB pipeline for consistency and compatibility.

// HOW TO USE
// Add the paths for the run file and the output directory to pipeline.sh, then submit the job.
// The run file can be generated automatically using the runfile.sh script. Requires a list of GSEs to subset from the main file.
// The only rules for the run file are:
// 1. Include all lines with the same combination of GSE, species, and assay. Or for simplicity, just include all lines for a given GSE. This ensures all the data needed is present and will be grouped together in the final output. 
// 2. Avoid processing too many samples at once (keep under 100). The xenbase directory only has 5 TB of storage, so depending on how much is already in use, the limit can be reached and cause the pipeline to fail.

// IMPORTANT NOTES
// All results needed for Xenbase will be copied to the specified output directory. Once that has been checked the work directory can be safely deleted.
// If the pipeline fails check the logs/report as well as the work directory of any processes that failed. The most likely issues are either memory/space related or a problem with the run file.
// After a failure you can resume the pipeline part way through using the "RESUME" option in pipeline.sh. 
// The "nextflow.config" file controls how jobs are submitted for each process. Default values have been chosen that should optimize the speed and amount of jobs that can be run at once, but the file can be modified as needed to provide additional memory or CPUs. 
// This pipeline performs automatic cleanup of files when they are no longer needed. They may still appear in the work directories but are replaced by sparse files that take up no space. This workaround "tricks" Nextflow into thinking the cached results are still there so the resume function can still be used.