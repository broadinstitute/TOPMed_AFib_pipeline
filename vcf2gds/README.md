# vcf2gds
### who am I ????
check the wdl scripts\n

#### check wether my wdl script is correct
java -jar /Users/schoi/cromwell/womtool-38.jar validate vcfToGds.wdl\n

#### what are my input files
java -jar /Users/schoi/cromwell/womtool-38.jar inputs vcfToGds.wdl > myWorkflow_inputs.json

#### excute!!!!
java -jar /Users/schoi/cromwell/cromwell-38.jar run vcfToGds.wdl --inputs test_inputs.json\n
