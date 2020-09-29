task annote_loftee {

	File vcffile
  String outfile
  Int disk
	Float memory
	Int cpus

	command {

#### copy a large file
curl https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw -o /opt/vep/.vep/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw

  perl ./vep \
  -i ${vcffile} \
  -o ${outfile} \
  --format vcf \
  --compress_output gzip \
  --assembly GRCh38 --species homo_sapiens \
  --offline --cache \
  --no_stats \
  --canonical --tab \
  --dir_plugins /opt/vep/.vep/Plugins/loftee/ \
  --force_overwrite \
  --plugin LoF,loftee_path:/opt/vep/.vep/Plugins/loftee/,gerp_bigwig:/opt/vep/.vep/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/opt/vep/.vep/loftee_data/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee_data/loftee.sql


	}

	runtime {
		docker: "schoi/vep_loftee:101.0"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
		cpu : "${cpus}"
		bootDiskSizeGb: 70
}

	output {
		File annotated_file = "${outfile}"

	}
}

workflow vep_loftee {
	Array[File] vcf_files
  Int this_disk
	Int this_cpus
	Float this_memory

	scatter(chrfile in vcf_files) {
 		call annote_loftee {
			input: vcffile = chrfile,
			disk = this_disk,
			memory = this_memory,
			cpus = this_cpus

		}
	}
