task ageatonset {

	File phen
	String chrom
	File gds
	File genefile
	File resultfile
	String unrelcol
	String noHFcol
	Int disk
	Float memory
	Int cpus

	String out_base = "Chr" + chrom + "_ageatonset_prevalence_results.RData"
	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

##### Perform Firth logistic regression
R CMD BATCH "--args ${phen} ${chrom} ${gds} ${genefile} ${resultfile} ${unrelcol} ${noHFcol} ${out_base} " ./TOPMed_AFib_pipeline/GENESIS/ageatonset/ageatonset_unrelated_noHF_before_AF.R ${out_base}.out

	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1.4.1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
		bootDiskSizeGb: 50
}

	output {
		File out_file0 = "${out_base}.out"
		File out_file1 = "${out_base}"
	}
}


task combine_result {

	Array[File] resultfiles
	Int disk
	Float memory
	Int cpus

	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

##### Perform combine results
R CMD BATCH "--args ${sep="," resultfiles}" ./TOPMed_AFib_pipeline/GENESIS/ageatonset/ageatonset_unrelated_noHF_before_AF_summary.R summary.out

	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1.4.1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
		bootDiskSizeGb: 50
}

	output {
	File summary_file = "summary.RData"
	File out_file = "summary.out"
	}
}


workflow ageatonset_workflow {
	File Input_arrayfiles
	Array[Array[File]] Inputfiles = read_tsv(Input_arrayfiles)
	File this_genefile
	File this_phen
	String this_unrelcol
	String this_noHFcol
	Int this_disk
	Int this_cpus
	Float this_memory

	scatter(chrfile in Inputfiles) {
		call ageatonset {
			input: chrom = chrfile[0],
			gds = chrfile[1],
			resultfile = chrfile[2],
			genefile = this_genefile,
			phen = this_phen,
			unrelcol = this_unrelcol,
			noHFcol = this_noHFcol,
			disk = this_disk,
			memory = this_memory,
			cpus = this_cpus

		}
	}

	call combine_result {
	input: resultfiles = ageatonset.out_file1,
	disk = this_disk,
	memory = this_memory,
	cpus = this_cpus

	}

	output {
		File out_files = combine_result.out_file
		File sum_file = combine_result.summary_file

	}
}
