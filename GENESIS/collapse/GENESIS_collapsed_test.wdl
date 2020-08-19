task collapsed_test {
	String chrom
	File gds
	File varlist
	File group
	File phen
	File nulmod
	String stat
	Float cutoff
	Int disk
	Float memory
	Int cpus

	String out_base = chrom + "_collapsed_results.RData"


	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

#### perform collapsed test
R CMD BATCH "--args ${chrom} ${gds} ${varlist} ${group} ${phen} ${nulmod} ${stat} ${cutoff} ${out_base} " ./TOPMed_AFib_pipeline/GENESIS/collapse/TOPMed_freeze8_af_hclof_collapsed.R ${out_base}.out

	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1.4.1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
		bootDiskSizeGb: 50
}

	output {
		File out_file1 = "${out_base}"
		File out_file2 = "${out_base}.out"
	}
}

workflow rare_variant_test {
	File Input_arrayfiles
	Array[Array[File]] Inputfiles = read_tsv(Input_arrayfiles)
	File this_phen
	File this_null
	Float this_cutoff
	String this_stat
	Int this_disk
	Int this_cpus
	Float this_memory

	scatter(chrfile in Inputfiles) {
		call collapsed_test {
			input: chrom = chrfile[0],
			gds = chrfile[1],
			varlist = chrfile[2],
			group = chrfile[3],
			phen = this_phen,
			nulmod = this_null,
			cutoff = this_cutoff,
			stat = this_stat,
			disk = this_disk,
			memory = this_memory,
			cpus = this_cpus
		}
	}


	output {
		Array[File] result_files = collapsed_test.out_file1
		Array[File] out_files = collapsed_test.out_file2
	}
}
