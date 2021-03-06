task runGds {
	String chr
	File bed
	File fam
	File bim
	File script
	Int disk
	Float memory
	Int ncpus

	String out_base = basename(bed, ".bed")


	command {

#### onvert vcf to gds
R CMD BATCH "--args ${bed} ${fam} ${bim} ${ncpus}" ${script} ${out_base}.out

	}

	runtime {
		docker: "quay.io/biocontainers/bioconductor-seqarray:1.30.0--r40h5f743cb_0"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${ncpus}"
		bootDiskSizeGb: 200
}

	output {
		File out_file1 = "${out_base}.gds"
		File out_file2 = "${out_base}.out"
	}
}

workflow makegds {
	File Input_arrayfiles
	Array[Array[File]] Inputfiles = read_tsv(Input_arrayfiles)
	File R_script
	Int this_disk
	Float this_memory
	Int this_cpus
	scatter(chrfile in Inputfiles) {
		call runGds {
			input: chr = chrfile[0],
			bed = chrfile[1],
			fam = chrfile[2],
			bim = chrfile[3],
			script = R_script,
			disk = this_disk,
			memory = this_memory,
			ncpus = this_cpus
		}
	}


	output {
		Array[File] gds_files = runGds.out_file1
		Array[File] out_files = runGds.out_file2
	}
}
