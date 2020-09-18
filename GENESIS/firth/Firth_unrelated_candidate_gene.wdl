task firth_test {

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

	String out_base = "Chr" + chrom + "_Firth_logistic_results.RData"
	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

##### Perform Firth logistic regression
R CMD BATCH "--args ${phen} ${chrom} ${gds} ${genefile} ${resultfile} ${unrelcol} ${noHFcol} ${out_base} " ./TOPMed_AFib_pipeline/GENESIS/firth/Firth_logistic_regression_unrelated_noHF_before_AF.R ${out_base}.out

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


workflow firth_test_scatter {
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
		call firth_test {
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

	output {
		Array[File] out_files = firth_test.out_file0
		Array[File] sum_file = firth_test.out_file1

	}
}
