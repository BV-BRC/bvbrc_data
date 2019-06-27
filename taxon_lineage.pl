#!/usr/bin/perl

use Bio::Taxon;
use Bio::Tree::Tree;

open IN, $ARGV[0];
open OUT, $ARGV[1];

my $count = 0;
while (my $entry = <IN>){

	$count++;
	
	my ($genome_id, $taxon_id) = $entry=~/(.*)\t(.*) *\n/;

	my $dbh = Bio::DB::Taxonomy->new(-source => 'entrez');
	$taxon = $dbh->get_taxon(-taxonid => '$taxon_id');

	my $tree = Bio::Tree::Tree->new();
	my @lineage = ($tree->get_lineage_nodes($taxon), $taxon);

	#foreach $node (@lineage) { print $node->id, ";"; } print "\n";

	#foreach $node (@lineage) { print $node->scientific_name, ";"; } print "\n";

	@taxon_lineage_id = "";
	@taxon_lineage_name = "";

	foreach $node (@lineage) { 
		push @taxon_lineage_id, $node->id; 
	} 
	print "\n";


	foreach $node (@lineage) { 
		push @taxon_lineage_name, $node->scientific_name; 
	} 

	$taxon_lineage_id = "\"", join('\",\"', @taxon_lineage_id), "\"";
	$taxon_lineage_name = "\"", join('\",\"', @taxon_lineage_name), "\"";
	
	print ",\n" unless $count == 1;
	print OUT "{\"genome_id\":\"$genome_id\",\"taxon_lineage_ids\":{\"set\":[@taxon_lineage_id]},\"taxon_lineage_names\":{\"set\":[@taxon_lineage_name]}}";

}

close IN;
close OUT;
