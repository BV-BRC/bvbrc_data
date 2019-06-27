#!/usr/bin/perl

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;&rows=100000";

my $core = "/genome";
my $query = "/select?q=genome_id:*.*";
my $fields = "&fl=genome_id";
my $sort = "";
my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

push @genomes, `wget -q -O - "$solrQuery"`;

foreach my $genome_id (@genomes){

		chomp $genome_id;
		next if $genome_id=~/genome_id/;

		print "Processing genome_id: $genome_id\n";

		my @features;

		my $core = "/genome_feature";
		my $query = "/select?q=genome_id:$genome_id";
		my $fields = "&fl=feature_id";
		my $sort = "";
		my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

		push @features, `wget -q -O - "$solrQuery"`;

		open OUT, ">$genome_id.json";

		my $count = 0;
		print OUT "[\n";

		foreach my $feature_id (@features) {
			
			chomp $feature_id;
			next if $feature_id=~/feature_id/;

			$count++;
			print OUT ",\n" unless $count == 1;
			print OUT "{\"feature_id\":\"$feature_id\",\"document_type\":{\"set\":\"genome_feature\"}}";
		
		}

		print OUT "]";

		close OUT;

		`time post.update.dayhoff.sh $core $genome_id.json`;

		`rm $genome_id.json`;


}
