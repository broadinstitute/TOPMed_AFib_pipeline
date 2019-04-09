task runGds {
	File vcf
	Int disk
	Float memory
	Int cpus

	String out_base = basename(vcf, ".vcf.gz")


	command {

#### clone the pipeline
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

#### call header
bcftools view -h ${vcf} > header.txt

#### modify the header
sed s/"BETA_IF,Number=%d"/"BETA_IF,Number=."/g header.txt > header2.txt

#### modify the vcf file
bcftools reheader -h header2.txt -o tmp.vcf.gz ${vcf}

#### move tmp file to vcf
mv tmp.vcf.gz ${vcf}

#### onvert vcf to gds
R CMD BATCH "--args ${vcf} ${out_base} ${cpus}" ./TOPMed_AFib_pipeline/vcf2gds/vcf2gds.R ${out_base}.out

	}

	runtime {
		docker: "gcr.io/broad-afib-pipeline/gmmat:v1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
		bootDiskSizeGb: 50
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
