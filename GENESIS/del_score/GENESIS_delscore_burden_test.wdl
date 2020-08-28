task delscore_test {
	File dbnsfp
	String score
	String scorevar
	Float mintool
	String chrom
	File gds
	File varlist
	File group
	File phen
	File nulmod
	Float scutoff
	Float acutoff
	String stat
	Int disk
	Float memory
	Int cpus

	String out_base = "Chr" + chrom + "_delscore" + scutoff + "_" + stat + "_results.RData"

	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

##### make scores
R CMD BATCH "--args ${dbnsfp} ${score} ${scorevar} ${mintool} " ./TOPMed_AFib_pipeline/GENESIS/del_score/dbNSFP_delscore.R ${score}_delscore.out

#### perform burden test
R CMD BATCH "--args ${chrom} ${gds} ${varlist} ${group} ${phen} ${nulmod} ${score} ${scorevar} ${scutoff} ${acutoff} ${stat} ${out_base} " ./TOPMed_AFib_pipeline/GENESIS/del_score/dbNSFP_delscore_hclof_missense_burden.R ${out_base}.out

	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1.4.1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
		bootDiskSizeGb: 50
}

	output {
		File out_file0 = "${score}"
		File out_file1 = "${out_base}"
		File out_file2 = "${out_base}.out"
	}
}


task test_summary {
	Array[File] resultfiles
	Float cmaccutoff
	Int disk
	Float memory
	Int cpus

	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

#### perform collapsed test
R CMD BATCH "--args ${cmaccutoff} ${sep="," resultfiles}" ./TOPMed_AFib_pipeline/GENESIS/collapse/TOPMed_freeze8_af_hclof_collapsed_summary.R summary.out

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
		File manhattan_plot = "manhattan_plot.png"
		File qq_plot = "qqplot.png"

	}
}

workflow rare_variant_delscore_test {
	File Input_arrayfiles
	Array[Array[File]] Inputfiles = read_tsv(Input_arrayfiles)
	String this_scorevar
	Float this_mintool
	File this_phen
	File this_null
	Float this_scorecutoff
	Float this_allelecutoff
	String this_stat
	Float this_mincmac
	Int this_disk
	Int this_cpus
	Float this_memory

	scatter(chrfile in Inputfiles) {
		call delscore_test {
			input: chrom = chrfile[0],
			gds = chrfile[1],
			varlist = chrfile[2],
			group = chrfile[3],
			dbnsfp = chrfile[4],
			score = chrfile[5],
			scorevar = this_scorevar,
			scutoff = this_scorecutoff,
			acutoff = this_allelecutoff,
			mintool = this_mintool,
			phen = this_phen,
			nulmod = this_null,
			stat = this_stat,
			disk = this_disk,
			memory = this_memory,
			cpus = this_cpus


		}
	}

		call test_summary {
		input: resultfiles = delscore_test.out_file1,
		cmaccutoff = this_mincmac,
		disk = this_disk,
		memory = this_memory,
		cpus = this_cpus

		}

	output {
		Array[File] score_files = delscore_test.out_file0
		Array[File] result_files = delscore_test.out_file1
		Array[File] out_files = delscore_test.out_file2
		File sum_file = test_summary.summary_file
		File sumout_file = test_summary.out_file
		File man_file = test_summary.manhattan_plot
		File qq_file = test_summary.qq_plot

	}
}
