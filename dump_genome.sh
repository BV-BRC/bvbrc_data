
###############################################################################
# 
# Purpose: Simple shell script to dump all the solr data for a genome
# 
# Usage: dump_genome.sh 83332.12 
#
###############################################################################

mkdir $1
pushd $1
wget "$PATRIC_SOLR/genome/select?q=genome_id%3A$1&wt=json&indent=false&omitHeader=true&rows=100000" -O genome.json
wget "$PATRIC_SOLR/genome_sequence/select?q=genome_id%3A$1&wt=json&indent=false&omitHeader=true&rows=100000" -O genome_sequence.json
wget "$PATRIC_SOLR/genome_feature/select?q=genome_id%3A$1+AND+annotation%3APATRIC&wt=json&indent=false&omitHeader=true&rows=100000" -O genome_feature_patric.json
wget "$PATRIC_SOLR/genome_feature/select?q=genome_id%3A$1+AND+annotation%3ARefSeq&wt=json&indent=false&omitHeader=true&rows=100000" -O genome_feature_refseq.json
wget "$PATRIC_SOLR/genome_feature/select?q=genome_id%3A$1+AND+annotation%3ABRC1&wt=json&indent=false&omitHeader=true&rows=100000" -O genome_feature_brc1.json
wget "$PATRIC_SOLR/pathway/select?q=genome_id%3A$1&wt=json&indent=false&omitHeader=true&rows=100000" -O pathway.json
wget "$PATRIC_SOLR/sp_gene/select?q=genome_id%3A$1&wt=json&indent=false&omitHeader=true&rows=100000" -O sp_gene.json
perl -pi -e 's/^{\n|^\s*"response".*"docs":|^\s*}}\n$//' *.json
popd $1

