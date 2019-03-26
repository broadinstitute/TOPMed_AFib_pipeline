workflow vcf_filter {
	Array[File] vcf_files
	File this_sample
  Float this_ancount
  Int this_account
  Int this_disk
	Float this_memory

	scatter(this_file in vcf_files) {
		call ac_an_filter {
			input: vcf = this_file, sample = this_sample, ancount = this_ancount, acount = this_account, disk = this_disk, memory = this_memory
		}
	}


	output {
		Array[File] gds_files = this_ancount.out_file
	}
}

task ac_an_filter {
	File vcf
	File sample
  Float ancount
  Int account
	String out_base = basename(vcf, ".bcf")
  Int disk
	Float memory

	command {

bcftools view -Ou -S ${sample} ${vcf} | bcftools view -Oz -o ${out_base}.vcf.gz --include "INFO/AC>${account} & INFO/AN>${ancount}"

	}

	runtime {
		docker: "biocontainers/bcftools:v1.5_cv2"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
	}

	output {
		File out_file = "${out_base}.vcf.gz"
	}
}
