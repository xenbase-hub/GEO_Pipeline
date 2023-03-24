#!/usr/bin/env nextflow

process cleanWorkFiles {

  input:
  val(file)

  output:
  val(1), emit: IS_CLEAN

  script:
  """
    ${workflow.projectDir}/bin/clean_work_files.sh "${file}"
  """
}