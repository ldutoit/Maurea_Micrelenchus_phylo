#Kirsten#





to **independly explore** phlogenies of two groups of NZ molluscs *Maurea* and *Micrelenchus* 


This build on knowledge gained from : 

```
#private repos need access granted by ldutoit
https://github.com/ldutoit/Micrelenchus
https://github.com/ldutoit/Maurea_mito
```

## Workflow in four steps:


1. adapter removal and renaming

2. Map onto the longest of each of the species reference genome and call the variants


3. Make consensus fasta


4. Filter for positions that do not have any coverahe on an individual basis.
Bonus:

use MITOBIM and try to get longer fragment, or the whole mitogenomes, repeat 2-6 then.


# 1 adapter removal and samples renaming

```python
# 1 1. adapter removal ans ashorten 
import os

input_folder = "source_files/28samples/OG4710/OG4710_fastq/"
#os.mkdir("trimmed_reads")
samples= list(set([ x.split("L001")[0] for x in os.listdir(input_folder) if x.endswith("gz")]))

i=0
for sample in samples:
	i+=1
	print (i)
	os.system("cutadapt -m 50  --length 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_reads/"+sample+"1.fq -p trimmed_reads/"+sample+"2.fq  "+input_folder+sample+"L001_R1_001.fastq.gz "+input_folder+sample+"L001_R2_001.fastq.gz") 


#renaming samples
matching_dict= {"4710-01":"Micrelenchus_tesselatus4710-01","4710-02":"Micrelenchus_tesselatus4710-02","4710-03":"Micrelenchus_tesselatus4710-03","4710-04":"Micrelenchus_tesselatus4710-04","4710-05":"Micrelenchus_tesselatus4710-05","4710-06":"Micrelenchus_huttonii4710-06","4710-07":"Micrelenchus_huttonii4710-07","4710-08":"Micrelenchus_huttonii4710-08","4710-09":"Micrelenchus_huttonii4710-09","4710-10":"Micrelenchus_huttonii4710-10","4710-11":"Micrelenchus_tenebrosus4710-11","4710-12":"Maurea_osbornei4710-12","4710-13":"Maurea_tigris4710-13","4710-14":"Maurea_tigris4710-14","4710-15":"Maurea_turnerarum4710-15","4710-16":"Alertalax_blacki4710-16","4710-17":"Maurea_sanguineus4710-17","4710-18":"Maurea_sanguineus4710-18","4710-20":"Maurea_foveauxanum4710-20","4710-21":"Maurea_blacki4710-21","4710-22":"Maurea_blacki4710-22","4710-23":"Calliostoma_conulus4710-23","4710-24":"Alertalax_blacki4710-24","4710-26":"Micrelenchus_tenebrosus4710-26","4710-28":"Maurea_osbornei4710-28","4710-29":"Maurea_tigris4710-29","4710-30":"Maurea_pellucida_spirata4710-30","4710-31":"Maurea_gibbsorum4710-31"}

os.chdir("trimmed_reads")
for key in matching_dict.keys():
	torename=[x for x in os.listdir(".") if key in x]
	for file in torename:
		os.rename(file,matching_dict[key]+file[-5:])
		
```
	
## 2 mapping

### 2a mapping Micrelenchus

```sh
module load SAMtools BWA picard/2.1.0

```
```python
import os
os.system("mkdir -p simplemapping/tesselatus")
os.chdir("simplemapping/tesselatus")
	

#sample [forwardreads, backwardreads, eference]
sample_dict = {"Micrelenchus_huttonii4710-06": ["../../trimmed_reads/Micrelenchus_huttonii4710-06_1.fq",
"../../trimmed_reads/Micrelenchus_huttonii4710-06_2.fq","ref.fa"],
"Micrelenchus_huttonii4710-07":[
"../../trimmed_reads/Micrelenchus_huttonii4710-07_1.fq",
"../../trimmed_reads/Micrelenchus_huttonii4710-07_2.fq","ref.fa"],
"Micrelenchus_huttonii4710-08": ["../../trimmed_reads/Micrelenchus_huttonii4710-08_1.fq",
"../../trimmed_reads/Micrelenchus_huttonii4710-08_2.fq","ref.fa"],
"Micrelenchus_huttonii4710-09":["../../trimmed_reads/Micrelenchus_huttonii4710-09_1.fq",
"../../trimmed_reads/Micrelenchus_huttonii4710-09_2.fq","ref.fa"],
"Micrelenchus_huttonii4710-10":["../../trimmed_reads/Micrelenchus_huttonii4710-10_1.fq",
"../../trimmed_reads/Micrelenchus_huttonii4710-10_2.fq","ref.fa"],
"Micrelenchus_tenebrosus4710-11":["../../trimmed_reads/Micrelenchus_tenebrosus4710-11_1.fq",
"../../trimmed_reads/Micrelenchus_tenebrosus4710-11_2.fq","ref.fa"],
"Micrelenchus_tenebrosus4710-26":["../../trimmed_reads/Micrelenchus_tenebrosus4710-26_1.fq",
"../../trimmed_reads/Micrelenchus_tenebrosus4710-26_2.fq","ref.fa"],
"Micrelenchus_tesselatus4710-01":["../../trimmed_reads/Micrelenchus_tesselatus4710-01_1.fq",
"../../trimmed_reads/Micrelenchus_tesselatus4710-01_2.fq","ref.fa"],
"Micrelenchus_tesselatus4710-02":["../../trimmed_reads/Micrelenchus_tesselatus4710-02_1.fq",
"../../trimmed_reads/Micrelenchus_tesselatus4710-02_2.fq","ref.fa"],
"Micrelenchus_tesselatus4710-03":["../../trimmed_reads/Micrelenchus_tesselatus4710-03_1.fq",
"../../trimmed_reads/Micrelenchus_tesselatus4710-03_2.fq","ref.fa"],
"Micrelenchus_tesselatus4710-04":["../../trimmed_reads/Micrelenchus_tesselatus4710-04_1.fq",
"../../trimmed_reads/Micrelenchus_tesselatus4710-04_2.fq","ref.fa"],
"Micrelenchus_tesselatus4710-05":["../../trimmed_reads/Micrelenchus_tesselatus4710-05_1.fq",
"../../trimmed_reads/Micrelenchus_tesselatus4710-05_2.fq","ref.fa"],"Micrelenchus_sanguineus4710_17":["../../trimmed_reads/Micrelenchus_sanguineus4710-17_1.fq",
"../../trimmed_reads/Micrelenchus_sanguineus4710-17_2.fq","ref.fa"]}

sample_dict ={"Micrelenchus_sanguineus4710-18":["../../trimmed_reads/Micrelenchus_sanguineus4710-18_1.fq",
"../../trimmed_reads/Micrelenchus_sanguineus4710-18_2.fq","ref.fa"]}
#function ...
for sample in (sample_dict.keys()):
	read1,read2,ref = sample_dict[sample][0],sample_dict[sample][1],sample_dict[sample][2]
	os.system("samtools faidx %s" % (ref)) #indexing reference
	os.system("java -jar $PICARD CreateSequenceDictionary R=%s O=%s.dict" %(ref,ref))
	command ="#!/bin/sh\nmodule load BWA/0.7.10-goolf-1.5.14 SAMtools/0.1.19-foss-2015a\n"
	command+="bwa index -a is %s  \n" % (ref)
	command+="bwa mem  -T 40 %s %s %s > %s.sam\n" % (ref,read1,read2,sample) 
	command+="samtools view -Sb %s.sam > %s.bam" % (sample,sample)
	output = open("mapping"+sample+".sh","w")
	output.write(command)
	output.close()
	#os.system("sbatch -A uoo00116 -t 1-00:00  -J"+ sample+ " "+ "mapping"+sample+".sh") 	
	print (command)
```


```python
import os
os.system("mkdir -p simplemapping/maurea")
os.chdir("simplemapping/maurea")
	

#sample [forwardreads, backwardreads, eference]
sample_dict ={"Maurea_blacki4710-21":["../../trimmed_reads/Maurea_blacki4710-21_1.fq","../../trimmed_reads/Maurea_blacki4710-21_2.fq","ref.fa"],
"Maurea_blacki4710-22":["../../trimmed_reads/Maurea_blacki4710-22_1.fq",
"../../trimmed_reads/Maurea_blacki4710-22_2.fq","ref.fa"],
"Maurea_foveauxanum4710-20":["../../trimmed_reads/Maurea_foveauxanum4710-20_1.fq",
"../../trimmed_reads/Maurea_foveauxanum4710-20_2.fq","ref.fa"],
"Maurea_gibbsorum4710-31":["../../trimmed_reads/Maurea_gibbsorum4710-31_1.fq",
"../../trimmed_reads/Maurea_gibbsorum4710-31_2.fq","ref.fa"],
"Maurea_osbornei4710-12":["../../trimmed_reads/Maurea_osbornei4710-12_1.fq",
"../../trimmed_reads/Maurea_osbornei4710-12_2.fq","ref.fa"],
"Maurea_osbornei4710-28":["../../trimmed_reads/Maurea_osbornei4710-28_1.fq",
"../../trimmed_reads/Maurea_osbornei4710-28_2.fq","ref.fa"],
"Maurea_pellucida_spirata4710-30":["../../trimmed_reads/Maurea_pellucida_spirata4710-30_1.fq",
"../../trimmed_reads/Maurea_pellucida_spirata4710-30_2.fq","ref.fa"],
"Maurea_tigris4710-13":["../../trimmed_reads/Maurea_tigris4710-13_1.fq",
"../../trimmed_reads/Maurea_tigris4710-13_2.fq","ref.fa"],
"Maurea_tigris4710-14":["../../trimmed_reads/Maurea_tigris4710-14_1.fq",
"../../trimmed_reads/Maurea_tigris4710-14_2.fq","ref.fa"],
"Maurea_tigris4710-29":["../../trimmed_reads/Maurea_tigris4710-29_1.fq",
"../../trimmed_reads/Maurea_tigris4710-29_2.fq","ref.fa"],
"Maurea_turnerarum4710-15":["../../trimmed_reads/Maurea_turnerarum4710-15_1.fq",
"../../trimmed_reads/Maurea_turnerarum4710-15_2.fq","ref.fa"],
"Alertalax_blacki4710-16":["../../trimmed_reads/Alertalax_blacki4710-16_1.fq",
"../../trimmed_reads/Alertalax_blacki4710-16_2.fq","ref.fa"],
"Alertalax_blacki4710-24":["../../trimmed_reads/Alertalax_blacki4710-24_1.fq",
"../../trimmed_reads/Alertalax_blacki4710-24_2.fq","ref.fa"],
"Calliostoma_conulus4710-23":["../../trimmed_reads/Calliostoma_conulus4710-23_1.fq",
"../../trimmed_reads/Calliostoma_conulus4710-23_2.fq","ref.fa"]}

#function ...
for sample in (sample_dict.keys()):
	read1,read2,ref = sample_dict[sample][0],sample_dict[sample][1],sample_dict[sample][2]
	os.system("samtools faidx %s" % (ref)) #indexing reference
	os.system("java -jar $PICARD CreateSequenceDictionary R=%s O=%s.dict" %(ref,ref))
	command ="#!/bin/sh\nmodule load BWA/0.7.10-goolf-1.5.14 SAMtools/0.1.19-foss-2015a\n"
	command+="bwa index -a is %s  \n" % (ref)
	command+="bwa mem  -T 40 %s %s %s > %s.sam\n" % (ref,read1,read2,sample) 
	command+="samtools view -Sb %s.sam > %s.bam" % (sample,sample)
	output = open("mapping"+sample+".sh","w")
	output.write(command)
	output.close()
	#os.system("sbatch -A uoo00116 -t 1-00:00  -J"+ sample+ " "+ "mapping"+sample+".sh") 	
	print (command)
```


# Make consensus

We used this script:

```sh
#!/bin/sh
#Written by Ludovic dutoit, April 2018
#Make a simple consensus fastasequence using bcftools from a bam file and the reference fasta sequence on the nesi cluster
#usage: 
#MakeConsensus.sh bamfile reference output
#NOTE: this is a haploid consensus for mitochondria
module purge
module load SAMtools/1.8-gimkl-2017a BCFtools/1.8-gimkl-2017a picard/2.1.0 # loading modules on nesi

#file="remapping/ill01.bam"
#reference="remapping/S27ref.fasta" 
#outputfile="consensus/ill01.fa"

file=$1
reference=$2
outputfile=$3

echo "processing $1. with reference $2, to make consensus $3"


samtools sort $file > sorted.bam
java -jar -Xms128m -Xmx128m  $PICARD MarkDuplicates INPUT=sorted.bam OUTPUT=sorted.rmdup.bam METRICS_FILE=sorted.rmdup2.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1 TMP_DIR="temp"
ln -s /share/easybuild/RHEL6.3/westmere/software/BCFtools/1.8-gimkl-2017a/bin/vcfutils.pl .
samtools mpileup -uf $reference sorted.rmdup.bam | bcftools call   --ploidy 1  -mv -Oz -o calls.vcf.gz
tabix calls.vcf.gz
cat $reference | bcftools consensus calls.vcf.gz > $outputfile
rm sorted.rmdup.bam sorted.bam calls.vcf.gz
unlink vcfutils.pl
```

For all the samples:
```sh
cd simplemapping/maurea
#./makeconsensus.sh  Alertalax_blacki4710-16.bam ref.fa Alertalax_blacki4710-16.fa  
#./makeconsensus.sh Alertalax_blacki4710-24.bam ref.fa Alertalax_blacki4710-24.fa
#./makeconsensus.sh Calliostoma_conulus4710-23.bam ref.fa Calliostoma_conulus4710-23.fa
#./makeconsensus.sh Maurea_blacki4710-21.bam ref.fa Maurea_blacki4710-21.fa
#./makeconsensus.sh Maurea_blacki4710-22.bam ref.fa Maurea_blacki4710-22.fa
#./makeconsensus.sh Maurea_foveauxanum4710-20.bam ref.fa Maurea_foveauxanum4710-20.fa
#./makeconsensus.sh Maurea_gibbsorum4710-31.bam ref.fa Maurea_gibbsorum4710-31.fa
#./makeconsensus.sh Maurea_osbornei4710-12.bam ref.fa Maurea_osbornei4710-12.fa
#./makeconsensus.sh Maurea_osbornei4710-28.bam ref.fa Maurea_osbornei4710-28.fa
#./makeconsensus.sh Maurea_pellucida_spirata4710-30.bam ref.fa Maurea_pellucida_spirata4710-30.fa
#./makeconsensus.sh Maurea_tigris4710-13.bam ref.fa Maurea_tigris4710-13.fa
#./makeconsensus.sh Maurea_tigris4710-14.bam ref.fa Maurea_tigris4710-14.fa
#./makeconsensus.sh Maurea_tigris4710-29.bam ref.fa Maurea_tigris4710-29.fa
#./makeconsensus.sh Maurea_turnerarum4710-15.bam ref.fa Maurea_turnerarum4710-15.fa

cd ../tesselatus./makeconsensus.sh  ref.fa

./makeconsensus.sh Micrelenchus_huttonii4710-06.bam ref.fa  Micrelenchus_huttonii4710-06fa
./makeconsensus.sh Micrelenchus_huttonii4710-07.bam ref.fa  Micrelenchus_huttonii4710-07.fa
./makeconsensus.sh Micrelenchus_huttonii4710-08.bam ref.fa  Micrelenchus_huttonii4710-08.fa
./makeconsensus.sh Micrelenchus_huttonii4710-09.bam ref.fa  Micrelenchus_huttonii4710-09.fa


./makeconsensus.sh Micrelenchus_huttonii4710-10.bam ref.fa  Micrelenchus_huttonii4710-10.fa
./makeconsensus.sh Micrelenchus_tenebrosus4710-11.bam ref.fa  Micrelenchus_tenebrosus4710-11.fa
./makeconsensus.sh Micrelenchus_tenebrosus4710-26.bam ref.fa  Micrelenchus_tenebrosus4710-26.fa


./makeconsensus.sh Micrelenchus_tesselatus4710-01.bam ref.fa  Micrelenchus_tesselatus4710-01.fa
./makeconsensus.sh Micrelenchus_tesselatus4710-02.bam ref.fa  Micrelenchus_tesselatus4710-02.fa
./makeconsensus.sh Micrelenchus_tesselatus4710-03.bam ref.fa  Micrelenchus_tesselatus4710-03.fa
./makeconsensus.sh Micrelenchus_tesselatus4710-04.bam ref.fa  Micrelenchus_tesselatus4710-04.fa
./makeconsensus.sh Micrelenchus_tesselatus4710-05.bam ref.fa  Micrelenchus_tesselatus4710-05.fa

./makeconsensus.sh Micrelenchus_sanguineus4710-17.bam ref.fa Micrelenchus_sanguineus4710017.fa
./makeconsensus.sh Micrelenchus_sanguineus4710-18.bam ref.fa Micrelenchus_sanguineus4710-18.fa

```

Explore depth and remove non informative positions


```sh
#Micrelenchus 
samtools sort  Micrelenchus_huttonii4710-06.bam > pouet; samtools depth pouet > Micrelenchus_huttonii4710-06.depth
samtools sort  Micrelenchus_huttonii4710-07.bam > pouet; samtools depth pouet >  Micrelenchus_huttonii4710-07.depth
samtools sort  Micrelenchus_huttonii4710-08.bam > pouet; samtools depth pouet >   Micrelenchus_huttonii4710-08.depth
samtools sort  Micrelenchus_huttonii4710-09.bam > pouet; samtools depth pouet >  Micrelenchus_huttonii4710-09.depth
samtools sort  Micrelenchus_huttonii4710-10.bam > pouet; samtools depth pouet >  Micrelenchus_huttonii4710-10.depth
samtools sort  Micrelenchus_tenebrosus4710-11.bam > pouet; samtools depth pouet >  Micrelenchus_tenebrosus4710-11.depth
samtools sort  Micrelenchus_tenebrosus4710-26.bam > pouet; samtools depth pouet >  Micrelenchus_tenebrosus4710-26.depth
samtools sort  Micrelenchus_tesselatus4710-01.bam > pouet; samtools depth pouet >   Micrelenchus_tesselatus4710-01.depth
samtools sort  Micrelenchus_tesselatus4710-02.bam > pouet; samtools depth pouet >  Micrelenchus_tesselatus4710-02.depth
samtools sort  Micrelenchus_tesselatus4710-03.bam > pouet; samtools depth pouet >  Micrelenchus_tesselatus4710-03.depth
samtools sort  Micrelenchus_tesselatus4710-04.bam > pouet; samtools depth pouet >  Micrelenchus_tesselatus4710-04.depth
samtools sort  Micrelenchus_tesselatus4710-05.bam > pouet; samtools depth pouet >  Micrelenchus_tesselatus4710-05.depth

samtools sort   Micrelenchus_sanguineus4710-17.bam > pouet; samtools depth pouet >  Micrelenchus_sanguineus4710-17.depth
samtools sort   Micrelenchus_sanguineus4710-18.bam > pouet; samtools depth pouet >  Micrelenchus_sanguineus4710-18.depth

wc -l *depth

```

##Replace all non covered positions

```python
from Bio import SeqIO

samples= [
"Micrelenchus_huttonii4710-06",
"Micrelenchus_huttonii4710-07",
"Micrelenchus_huttonii4710-08",
"Micrelenchus_huttonii4710-09",
"Micrelenchus_huttonii4710-10",
"Micrelenchus_tenebrosus4710-11",
"Micrelenchus_tenebrosus4710-26",
"Micrelenchus_tesselatus4710-01",
"Micrelenchus_tesselatus4710-02",
"Micrelenchus_tesselatus4710-03",
"Micrelenchus_tesselatus4710-04",
"Micrelenchus_tesselatus4710-05",
"Micrelenchus_sanguineus4710-17",
"Micrelenchus_sanguineus4710-18"]

for sample in samples:
	covered_bases=0
	#print(sample)
	covered =[[line.strip().split()[0],line.strip().split()[1]] for line in open(sample+".depth")]
	fasta_sequences = SeqIO.parse(open(sample+".fa"),'fasta')
	output=open(sample+"withmissing.fa","w")
	for fasta in fasta_sequences:
		name, sequence = fasta.id, list(fasta.seq)
		for i,base in enumerate(sequence):
			#print (i,base)
			if not [name, str(i+1)] in covered: 
				sequence[i]="N"
			else:
				covered_bases+=1
			#print(sequence[i] )
		if "c2" in name:
			output.write(">"+sample+"\n"+"".join(sequence).upper())
		else:
			output.write("N"*100+"".join(sequence).upper()+"\n")
	print (sample,covered_bases)
output.close()





```


#Group them

```
cat *withmissing* > all_micrelenchus_1contig.fa
 ```

The total sequence is X long but there is Y in the middle

# BONUS quick run mitobim with 

input  from [seed.fa](https://github.com/ldutoit/Micrelenchus/blob/master/output_files/consensusMtesselatusNOEDGESCOMBINED.fa)


```sh
module load Python/2.7.16-gimkl-2018b
mkdir -p mitobim/tigris
cd mitobim/tigris
#https://github.com/ldutoit/Micrelenchus/blob/master/output_files/consensusMtesselatusNOEDGESCOMBINED.fa as seed.fa
module load seqtk
seqtk sample -s12  ../../trimmed_reads/Maurea_tigris4710-13_1.fq 1 >  Maurea_tigris4710-13_R1_sampled05.fq
seqtk sample -s12  ../../trimmed_reads/Maurea_tigris4710-13_2.fq 1 >  Maurea_tigris4710-13_R2_sampled05.fq
paste Maurea_tigris4710-13_R1_sampled05.fq Maurea_tigris4710-13_R2_sampled05.fq  | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}'   > Maurea_tigris4710-13_INTERLEAVED_sampled05.fastq
MITObim.pl --denovo -end 100 --pair -sample Mtesselatus -ref pilotstudytesselatus -readpool  Maurea_tigris4710-13_INTERLEAVED_sampled05.fastq --clean --quick seed.fa -kbait 15

rerunnign smaller proportion but on the output of run 1 as seed after checking that the ong ine is indeed mitochondrion
```



```sh
module load Python/2.7.16-gimkl-2018b
cd /home/ludovic.dutoit/projects/kirsten_28samples
mkdir -p mitobim/tesselatus
cd mitobim/tesselatus
#https://github.com/ldutoit/Micrelenchus/blob/master/output_files/consensusMtesselatusNOEDGESCOMBINED.fa as seed.fa
module load seqtk
seqtk sample -s12  ../../trimmed_reads/*4710-03*_1.fq 1 >  Micrelenchus_tesselatus4710-03_R1_sampled05.fq
seqtk sample -s12  ../../trimmed_reads/*4710-03*_2.fq 1 >  Micrelenchus_tesselatus4710-03_R2_sampled05.fq
paste Micrelenchus_tesselatus4710-03_R1_sampled05.fq Micrelenchus_tesselatus4710-03_R2_sampled05.fq  | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}'   > Micrelenchus_tesselatus4710-03_INTERLEAVED_sampled05.fastq
MITObim.pl --denovo -end 100 --pair -sample Mtesselatus -ref pilotstudytesselatus -readpool  Micrelenchus_tesselatus4710-03_INTERLEAVED_sampled05.fastq --clean --quick seed.fa -kbait 18

#rerunnign smaller proportion but on the output of run 1 as seed after checking that the ong ine is indeed mitochondrion

MITObim.pl --denovo -end 100 --pair -sample Mtesselatus -ref pilotstudytesselatus -readpool  Micrelenchus_tesselatus4710-03_INTERLEAVED_sampled05.fastq --clean --quick seed2.fa -kbait 18
```
