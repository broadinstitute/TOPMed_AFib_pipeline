task sampleqc {

File gdsfile
File varfile
File phenfile
String chr
Int chunknum
Int nvariants
Int disk
Float memory
String out_base  = basename(varfile,".variant.tsv")

command {

git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

R CMD BATCH "--args ${gdsfile} ${varfile} ${phenfile} ${chr} ${chunknum} ${nvariants} " ./sampleqc/sampleQC_chunk.R > ${out_base}_${chr}_${chunknum}.out

}

runtime {
  docker: "analysiscommon/genesis_wdl:v1.4.1"
  disks: "local-disk " + select_first([disk,"500"]) + " HDD"
  memory: select_first([memory,"30"]) + " GB"
  cpu : "${cpus}"
}

output {
File logfile = ${out_base}_sampleqc_chr${chr}_${chunknum}.out
File outfile = ${out_base}_sampleqc_chr${chr}_${chunknum}.tsv
}

#### WDL

workflow wgs_sampleqc {

# sampleqc procedures
File this_gdsfile
File this_varfile
File this_phenfile
String this_chr
Int this_chunksize
Int this_nvariants
Int? this_disk
Float? this_memory

Int this_totalchunknum = ceil(this_nvariants/this_chunksize)
Array[Int] this_chunknums = range(this_totalchunknum)+1


# Workflow metadata
	meta {
		description: "Sample QC for WGS per chromosome"
		tags: "Statistics"
	  author: "Seng Hoan Choi (WDL)"
	  email: "schoi@broadinstitute.org"
	}

parameter_meta {

this_gdsfile: "gds file created from SeqArray"
this_varfile: "tab delimited file, includes all variant id in gds file"
this_phenfile: "phenotype file for subtracting samples"
this_chr:  "chromosome number"
this_chunksize: "the number of variants in a chunk"
this_nvariants: "the number of variants in a varfile"
this_disk: "the size of the disk"
this_memory: "the size of the memory"

}

scatter(this_chunknum in this_chunknums) {
  call sampleqc {
  input:
  gdsfile = this_gdsfile
  varfile = this_varfile
  phenfile = this_phenfile
  chr = this_chr
  chunknum = this_chunknum
  nvariants = this_nvariants
  disk = this_disk
  memory = this_memory
}
}


output {
		Array[File] sampleqc_log_files = sampleqc.logfile
		Array[File] sampleqc_out_files = sampleqc.outfile
	}
}


}

File gdsfile
Int totalvariant
Int chunksize
Int totalchunknum = ceil(totalvariant/chunksize)
Array[Int] chunknums = range(chunknum)+1
