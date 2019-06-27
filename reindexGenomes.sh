
###############################################################################
# 
# Purpose: Simple shell script to dump all the solr data for a genome
# 
# Usage: dump_genome.sh 83332.12 
#
###############################################################################

for genome in `cat $1`;

do

  echo ""
  echo "Processing $genome"
  echo ""

	#wget -q "$PATRIC_SOLR/genome/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=100000" -O $genome.genome.json
	#wget -q "$PATRIC_SOLR/genome_amr/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=100000" -O $genome.genome_amr.json
	#wget -q "$PATRIC_SOLR/genome_sequence/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=10000000" -O $genome.genome_sequence.json
	#wget -q "$PATRIC_SOLR/genome_feature/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=10000000" -O $genome.genome_feature.json
	#wget -q "$PATRIC_SOLR/pathway/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=100000" -O $genome.pathway.json
	#wget -q "$PATRIC_SOLR/subsystem/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=100000" -O $genome.subsystem.json
	wget -q "$PATRIC_SOLR/sp_gene/select?q=genome_id%3A$genome&fl=*&wt=json&indent=true&omitHeader=true&rows=100000" -O $genome.sp_gene.json

	perl -pi -e 's/^{\n|^\s*"response".*"docs":|^\s*}}\n$//' $genome.*.json
	perl -pi -e 's/^.*\"(_version_)\":.*},\n$/},/' $genome.*.json
	perl -pi -e 's/^.*\"(_version_)\":.*}]\n$/}]/' $genome.*.json
	perl -pi -e 's/^.*\"(_version_)\":.*,\n$//' $genome.*.json
	#perl -pi -e 's/^.*\"(brc1_cds|sequences|document_type|_version_)\":.*}]\n$/}]/' $genome.*.json

	#curl -H "Content-Type: application/json" "http://localhost:9005/solr/genome/update?commit=false" -d @$genome.genome.json
	#curl -H "Content-Type: application/json" "http://localhost:9005/solr/genome_amr/update?commit=false" -d @$genome.genome_amr.json
	#curl -H "Content-Type: application/json" "http://localhost:9005/solr/genome_sequence/update?commit=false" -d @$genome.genome_sequence.json
	#curl -H "Content-Type: application/json" "http://localhost:9005/solr/genome_feature/update?commit=false" -d @$genome.genome_feature.json
	#curl -H "Content-Type: application/json" "http://localhost:9005/solr/pathway/update?commit=false" -d @$genome.pathway.json
	#curl -H "Content-Type: application/json" "http://localhost:9005/solr/subsystem/update?commit=false" -d @$genome.subsystem.json
	curl -H "Content-Type: application/json" "http://localhost:9005/solr/sp_gene/update?commit=false" -d @$genome.sp_gene.json

	#post.update_dev.sh genome genome.json
	#post.update_dev.sh genome_amr genome_amr.json
	#post.update_dev.sh genome_sequence genome_sequence.json
	#post.update_dev.sh genome_feature genome_feature.json
	#post.update_dev.sh pathway pathway.json
	#post.update_dev.sh subsystem subsystem.json
	#post.update_dev.sh sp_gene sp_gene.json

	rm $genome.*.json

done
