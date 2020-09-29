task runGds {
	File bcf
	Int disk
	Float memory

	String out_base = basename(bcf, ".bcf")


	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

#### onvert vcf to gds
R CMD BATCH "--args ${bcf} ${out_base}" ./TOPMed_AFib_pipeline/vcf2gds/bcf2gds.R ${out_base}.out

	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1.4.1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		bootDiskSizeGb: 50
}

	output {
		File out_file1 = "${out_base}.gds"
		File out_file2 = "${out_base}.out"
	}
}

workflow makegds {
	Array[File] bcf_files
	Int this_disk
	Float this_memory

	scatter(this_file in vcf_files) {
		call runGds {
			input: vcf = this_file, disk = this_disk, memory = this_memory
		}
	}


	output {
		Array[File] gds_files = runGds.out_file1
		Array[File] out_files = runGds.out_file2
	}
}
