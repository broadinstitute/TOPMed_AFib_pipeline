task runGds {
	File vcf
	Int disk
	Float memory
	Int cpus

	String out_base = basename(vcf, ".vcf.gz")


	command {
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

R CMD BATCH "--args ${vcf} ${out_base} ${cpus}" ./TOPMed_AFib_pipeline/vcf2gds/vcf2gds.R > ${out_base}.out

	}

	runtime {
		docker: "gcr.io/broad-afib-pipeline/gmmat:v1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
	}

	output {
		File out_file1 = "${out_base}.gds"
		File out_file2 = "${out_base}.out"
	}
}

workflow makegds {
	Array[File] vcf_files
	Int this_disk
	Int this_cpus
	Float this_memory

	scatter(this_file in vcf_files) {
		call runGds {
			input: vcf = this_file, disk = this_disk, memory = this_memory, cpus = this_cpus
		}
	}


	output {
		Array[File] gds_files = runGds.out_file1
		Array[File] out_files = runGds.out_file2
	}
}
