# trim adapters and low quality reads

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --quality-cutoff=15,10 --minimum-length=36 -u 10 -U 10 -o Trimmed_"$xbase" -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log

#fastqc
fastqc --noextract --nogroup "$file"

#alignment
bowtie2 --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x mm10_reference_file -1 "$file" -2 "${file/R1/R2}" | samtools view -u - | samtools sort - > "${xbase%.*}".bam 
samtools index "${xbase%.*}".bam
samtools idxstats $bam > ${file}_idxstats

