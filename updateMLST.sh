echo "" > mlst.out
while IFS=$'\t' read -r genome_id genome_name
do
	echo "${genome_id}:${genome_name}"
	#assign_st_to_genome.pl -m ~/mlstdb -s 1 --org "${genome_name}" < /vol/patric3/downloads/genomes/${genome_id}/${genome_id}.fna >> ${genome_id}.st1
	assign_st_to_genome.pl -m ~/mlstdb-2019-0221 -s 1 --org "${genome_name}" < /homes/mshukla/patric/data/ucla/fasta/${genome_id}.fna >> ${genome_id}.st1
	add_mlstdb_tags.pl < ${genome_id}.st1 > ${genome_id}.st2
 	table_extract_col.pl 0,1,7,8 < ${genome_id}.st2 > ${genome_id}.st 
	grep -H "." ${genome_id}.st >> mlst.out 
	rm ${genome_id}.st*
done < $1
