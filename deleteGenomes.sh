# Purpose: Shell scrip to delete a single genome solr
# Input: genome id, example: 12345.1
# Result: Recods corresponding to the specified genomes will be deleted from folloing solr cores:
#					genome, genome_sequence, genome_feature, pathway, sp_gene, genome_amr
# Usage: deleteGenomeSolr.sh genome_id_list

for genome in `cat $1`;
 
do

	echo ""
	echo "Delete $genome"
	echo ""

	curl "$PATRIC_SOLR/genome/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

	curl "$PATRIC_SOLR/genome_sequence/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

	curl "$PATRIC_SOLR/genome_feature/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

	curl "$PATRIC_SOLR/pathway/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

	curl "$PATRIC_SOLR/subsystem/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

	curl "$PATRIC_SOLR/sp_gene/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

	curl "$PATRIC_SOLR/genome_amr/update?stream.body=<delete><query>genome_id:$genome</query></delete>&commit=false"

done
