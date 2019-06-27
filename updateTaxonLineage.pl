#!/usr/bin/perl

open OUT, ">update_taxon_lineage.json";
open TAXON, ">update_taxon_count.json";

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;&rows=1000000";

my %taxon_counts = ();

my $core = "/genome";
my $query = "/select?q=genome_id:1733.* AND public:1";
my $fields = "&fl=genome_id,taxon_id";
my $sort = "";
my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

print "$solrQuery\n";

my @rows = `wget -q -O - "$solrQuery"`;

print OUT "[\n";

foreach my $row (@rows){

		chomp $row;
		next if $row=~/genome_id/;
		$count++;

		my ($genome_id, $taxon_id) = $row=~/(.*)\t(.*)/;

		print "Processing genome_id: $genome_id\t$taxon_id\n";

		my @taxon_lineage;

		my $core = "/taxonomy";
		my $query = "/select?q=taxon_id:$taxon_id";
		my $fields = "&fl=taxon_id,lineage_ids,lineage_names,lineage_ranks";
		my $sort = "";
		my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

		my @results = `wget -q -O - "$solrQuery"`;

		print "\tNot found: $taxon_id\n" if scalar @results < 2; 
		print "\tMultipe records found: $taxon_id\n" if scalar @results > 2; 

		#next;

		foreach my $result (@results) {
			
			chomp $result;
			next if $result=~/taxon_id/;

			my ($id, $lineage_ids, $lineage_names, $lineage_ranks) = $result=~/(.*)\t(.*)\t(.*)\t(.*)/;

			my @lineage_ids = split(/;/, $lineage_ids);
			my @lineage_names = split(/;/, $lineage_names);
			my @lineage_ranks = split(/;/, $lineage_ranks);

			my $taxon_lineage_ids = "\"".join('","', @lineage_ids)."\"";
			my $taxon_lineage_names = "\"".join('","', @lineage_names)."\"";

			my $cnt=0;
			my ($kingdom,$phylum,$class,$order,$family,$genus,$species);
			foreach my $rank (@lineage_ranks){
				$kingdom = $lineage_names[$cnt] if $rank=~/kingdom/i;
				$phylum = $lineage_names[$cnt] if $rank=~/phylum/i;
				$class = $lineage_names[$cnt] if $rank=~/class/i;
				$order = $lineage_names[$cnt] if $rank=~/order/i;
				$family = $lineage_names[$cnt] if $rank=~/family/i;
				$genus = $lineage_names[$cnt] if $rank=~/genus/i;
				$species = $lineage_names[$cnt] if $rank=~/species/i;
				$cnt++;
			}

			print OUT ",\n" unless $count == 1;
			print OUT "{\"genome_id\":\"$genome_id\",".
								"\"taxon_lineage_ids\":{\"set\":[$taxon_lineage_ids]},".
								"\"taxon_lineage_names\":{\"set\":[$taxon_lineage_names]},".
								"\"kingdom\":{\"set\":\"$kingdom\"},".	
								"\"phylum\":{\"set\":\"$phylum\"},".	
								"\"class\":{\"set\":\"$class\"},".	
								"\"order\":{\"set\":\"$order\"},".	
								"\"family\":{\"set\":\"$family\"},".	
								"\"genus\":{\"set\":\"$genus\"},".	
								"\"species\":{\"set\":\"$species\"}".	
								"}";
		

			foreach my $id (@lineage_ids){
				if ($taxon_counts{$id}>=1){
					$taxon_counts{$id}++;
				}else{
					$taxon_counts{$id} = 1;
				} 		
			}
		
	}

#last if $count >= 10; 

}

print OUT "]";
close OUT;

my $count = 0;
print TAXON "[\n";

foreach my $taxon_id (keys %taxon_counts){
	$count++;
	print TAXON ",\n" unless $count == 1;
	print TAXON "{\"taxon_id\":$taxon_id,\"genomes\":{\"set\":$taxon_counts{$taxon_id}}}";
}
print TAXON "]";
		
#`time post.update.dayhoff.sh $core $genome_id.json`;

