#!/usr/bin/env perl

open IN, $ARGV[0] or die "Can't open genome list file!!\n";

while (my $genome_id = <IN>){

	chomp $genome_id;

	print "Processing $genome_id\n";

	my $solrQuery  = "http://localhost:9006/solr/genome_feature/select?q=genome_id%3A$genome_id".
					"+AND+annotation%3APATRIC+AND+feature_type%3ACDS&fl=feature_id%2Cplfam_id".
					"&wt=csv&csv.separator=%09&rows=100000";

	my @features = `wget -q -O - "$solrQuery" | grep -v "feature_id"`;

	open OUT, ">$genome_id.plfam.json";

	print OUT "[";
	my $count = 0; 
  foreach my $feature (@features){
    $count++;
		my ($feature_id, $plfam_id) = $feature=~/(.*)\t(.*)/;

		next unless $plfam_id;

		if ($plfam_id){	
    	my ($plfam_prefix, $plfam_idx) = $plfam_id=~/(PLF_\d+)_(\d+)/;
    	$plfam_id = $plfam_prefix."_".sprintf("%08d", $plfam_idx);
		}

		print OUT ",\n" unless $count == 1;
		print OUT "{\"feature_id\":\"$feature_id\",\"plfam_id\":{\"set\":\"$plfam_id\"}}";

  } 

	print OUT "]";

	close OUT;

	#`post.update.sh genome_feature $genome_id.plfam.json`;

}

close IN;
