import pandas as pd
import sys
import os
import glob
args = sys.argv

# Script is primarily an exercise in writing Python
# Main benefit relative to R is it can be easier for text manipulation

# Read in the runfile and other arguments
runfile = pd.read_csv(str(args[1]), sep = '\t')
gse = str(args[2])
species = str(args[3])
assay = str(args[4])
species_name = ""

if species == "xenopus-laevis-v10":
	species_name = "XENLA_10.1"
elif species == "xenopus-trop-v10":
	species_name = "XENTR_10.0"

# Slight differences depending on assay

if assay == "RNA-SEQ":
	runfile_subset = runfile.loc[runfile['GSE'] == gse]
	runfile_subset = runfile_subset.loc[runfile_subset['SPECIES'] == species]
	runfile_subset = runfile_subset.loc[runfile_subset['ASSAYTYPE'] == assay]
	readme = runfile_subset[["BIN_ID","TRACKNAME","BIGWIG","CTRL_BIN_ID","TREATMENT","SRR_TREATMENT"]].copy()
	readme["CTRL_BIN_ID"] = readme["CTRL_BIN_ID"].astype(pd.Int64Dtype())
	treatment_srr = readme["SRR_TREATMENT"]
	de_files = ["_vs_" + s.split(",")[0] + ".txt" for s in treatment_srr]
	controls = readme.loc[readme['BIN_ID'].isin(readme["CTRL_BIN_ID"])]
	control_srr = readme["CTRL_BIN_ID"].map(dict(zip(controls["BIN_ID"],controls["SRR_TREATMENT"]))).fillna("NOCONTROL")
	de_files = [i.split(",")[0] + j for i, j in zip(control_srr, de_files)]
	readme = readme.drop(columns=["SRR_TREATMENT"])
	readme.set_axis(["Track #","Track Name","BigWig","Control Track #","GSMs"], axis=1, inplace=True)
	readme.insert(4,"DE-Analysis",de_files)
	readme.loc[readme["DE-Analysis"].str.contains("NOCONTROL"),"DE-Analysis"] = ""
	new_order = readme.index[readme["Track #"].isin(controls["BIN_ID"])].tolist() + readme.index[-readme["Track #"].isin(controls["BIN_ID"])].tolist()
	readme = readme.reindex(new_order)
	readme.to_csv("gsm_to_track.txt", sep="\t", index=False)

	with open("gsm_to_track.txt") as gsm_to_track, open("README.txt", 'w') as output:
		row = runfile_subset.iloc[[0]]
		output.write("GSE: " + row["GSE"].to_string(index=False) + "\n")
		output.write("GSE Title: " + row["GSE_TITLE"].to_string(index=False) + "\n")
		output.write("Xenbase Article Id: " + row["XB_ART_ID"].to_string(index=False) + "\n")
		output.write("Pubmed Id: " + row["PUBMED_ID"].to_string(index=False) + "\n")
		output.write("Download Location: " + "https://bigFrog.xenbase.org/xenbase/genomics/GEO/" + row["GSE"].to_string(index=False) + "/" + species_name + "/RNA-Seq/" + "\n")
		output.write("NGS Data: RNA-Seq" + "\n")
		output.write("Folder Descriptions:" + "\n")
		output.write("\t" + "BigWigs: .bw files in the folder represent binary format of mapped read intensity in the genome. These files can be used for visualization purposes." + "\n")
		output.write("\t" + "ExpressionFiles: Genes_TPM_Matrix.txt and Genes_Counts_Matrix.txt provide TPM and raw counts of genes across samples respectively." + "\n")
		output.write("\t" + "DE_Analysis: Differential Expression results between conditions. LogFC and FDR values can used to obtain differentially expressed genes." + "\n\n")
		for line in gsm_to_track:
		    output.write(line)

elif assay == "ChIP-TF" or assay == "ChIP-Epigenetic":
	runfile_subset = runfile.loc[runfile['GSE'] == gse]
	runfile_subset = runfile_subset.loc[runfile_subset['SPECIES'] == species]
	runfile_subset = runfile_subset.loc[runfile_subset['ASSAYTYPE'] == assay]
	readme = runfile_subset[["BIN_ID","TRACKNAME","BIGWIG","CTRL_BIN_ID","TREATMENT"]].copy()
	readme["CTRL_BIN_ID"] = readme["CTRL_BIN_ID"].astype(pd.Int64Dtype())
	bigwigs = readme["BIGWIG"]
	assays = runfile_subset["ASSAYTYPE"]
	assays = assays.replace(["ChIP-TF","ChIP-Epigenetic"],["_narrowPeaks.bed","_broadPeaks.bed"])
	peak_files = [s[:-3] for s in bigwigs]
	peak_files = peak_files + assays
	controls = readme.loc[readme['BIN_ID'].isin(readme["CTRL_BIN_ID"])]
	readme.set_axis(["Track #","Track Name","BigWig","Control Track #","GSMs"], axis=1, inplace=True)
	readme.insert(4,"Called Peaks",peak_files)
	new_order = readme.index[readme["Track #"].isin(controls["BIN_ID"])].tolist() + readme.index[-readme["Track #"].isin(controls["BIN_ID"])].tolist()
	readme = readme.reindex(new_order)
	readme.to_csv("gsm_to_track.txt", sep="\t", index=False)

	with open("gsm_to_track.txt") as gsm_to_track, open("README.txt", 'w') as output:
		row = runfile_subset.iloc[[0]]
		output.write("GSE: " + row["GSE"].to_string(index=False) + "\n")
		output.write("GSE Title: " + row["GSE_TITLE"].to_string(index=False) + "\n")
		output.write("Xenbase Article Id: " + row["XB_ART_ID"].to_string(index=False) + "\n")
		output.write("Pubmed Id: " + row["PUBMED_ID"].to_string(index=False) + "\n")
		output.write("Download Location: " + "https://bigFrog.xenbase.org/xenbase/genomics/GEO/" + row["GSE"].to_string(index=False) + "/" + species_name + "/ChIP-Seq/" + "\n")
		output.write("NGS Data: ChIP-Seq" + "\n")
		output.write("Folder Descriptions:" + "\n")
		output.write("\t" + "BigWigs: .bw files in the folder represent binary format of mapped read intensity in the genome. These files can be used for visualization purposes." + "\n")
		output.write("\t" + "Called_Peaks: Peak calls per sample are obtained using MACS2" + "\n\n")
		for line in gsm_to_track:
		    output.write(line)

elif assay == "ATAC":
	runfile_subset = runfile.loc[runfile['GSE'] == gse]
	runfile_subset = runfile_subset.loc[runfile_subset['SPECIES'] == species]
	runfile_subset = runfile_subset.loc[runfile_subset['ASSAYTYPE'] == assay]
	readme = runfile_subset[["BIN_ID","TRACKNAME","BIGWIG","CTRL_BIN_ID","TREATMENT"]].copy()
	readme["CTRL_BIN_ID"] = readme["CTRL_BIN_ID"].astype(pd.Int64Dtype())
	bigwigs = readme["BIGWIG"]
	peak_files = [s[:-3] + "_narrowPeaks.bed" for s in bigwigs]
	controls = readme.loc[readme['BIN_ID'].isin(readme["CTRL_BIN_ID"])]
	readme.set_axis(["Track #","Track Name","BigWig","Control Track #","GSMs"], axis=1, inplace=True)
	readme.insert(4,"Called Peaks",peak_files)
	new_order = readme.index[readme["Track #"].isin(controls["BIN_ID"])].tolist() + readme.index[-readme["Track #"].isin(controls["BIN_ID"])].tolist()
	readme = readme.reindex(new_order)
	readme.to_csv("gsm_to_track.txt", sep="\t", index=False)

	with open("gsm_to_track.txt") as gsm_to_track, open("README.txt", 'w') as output:
		row = runfile_subset.iloc[[0]]
		output.write("GSE: " + row["GSE"].to_string(index=False) + "\n")
		output.write("GSE Title: " + row["GSE_TITLE"].to_string(index=False) + "\n")
		output.write("Xenbase Article Id: " + row["XB_ART_ID"].to_string(index=False) + "\n")
		output.write("Pubmed Id: " + row["PUBMED_ID"].to_string(index=False) + "\n")
		output.write("Download Location: " + "https://bigFrog.xenbase.org/xenbase/genomics/GEO/" + row["GSE"].to_string(index=False) + "/" + species_name + "/ATAC-Seq/" + "\n")
		output.write("NGS Data: ATAC-Seq" + "\n")
		output.write("Folder Descriptions:" + "\n")
		output.write("\t" + "BigWigs: .bw files in the folder represent binary format of mapped read intensity in the genome. These files can be used for visualization purposes." + "\n")
		output.write("\t" + "Called_Peaks: Peak calls per sample are obtained using MACS2" + "\n\n")
		for line in gsm_to_track:
		    output.write(line)