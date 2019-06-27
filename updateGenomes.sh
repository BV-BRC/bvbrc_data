for genome in `cat $1`;

do

	echo ""
	echo "Processing $genome"
	echo ""

	perl -pi -e 's/commit=true/commit=false/' $genome.delete_pathway.sh
	chmod 755 $genome.delete_pathway.sh
	./$genome.delete_pathway.sh
	post.update.sh genome_feature $genome.genome_feature.json
	post.update.sh pathway $genome.pathway.json
	post.update.sh sp_gene $genome.sp_gene.json

done
