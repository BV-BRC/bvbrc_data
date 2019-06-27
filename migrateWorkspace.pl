#!/usr/bin/perl

# 
#
# Purpose: Script to migrate P2 user workspaces to P3 Identifier scheme. 
# Replaces P2 genome and feature IDs to corresponding P3 identofiers. 
# If no p3 ID found, it replaces p2 ID with null. 
#
# Input: Single JSON file containing all user workspace data.
# Output: Single JSON file, exactly the same as input file, except 
# 				it contains p3 ids. 
# Usage: migrateWorkspace.pl p2workspace.json p3workspace.json > log  
#
#


use JSON;
use Data::Dumper;

my $solr_server = $ENV{PATRIC_SOLR};

open IN, $ARGV[0];
my @json = <IN>;
my $json_text = join "", @json;
close IN;

my $jsonh = JSON->new->allow_nonref;
my $workspaces = $jsonh->decode($json_text);

foreach my $workspace (@{$workspaces}){

	foreach my $track (@{$workspace->{items}->{data}->{tracks}}){

		my $internal_id = $track->{internalId};
		my $p3ID;

		if ($track->{trackType} eq "Genome"){
			$p3ID = getP3GenomeID($internal_id) if $track->{trackType} eq "Genome";
		}elsif($track->{trackType} eq "Feature"){
			$p3ID = getP3FeatureID($internal_id) if $track->{trackType} eq "Feature";
		}

		if ($track->{trackType} eq "Genome" || $track->{trackType} eq "Feature") {
			print "$workspace->{ownerId}\t$track->{trackId}\t$track->{trackType}\t$internal_id\t$p3ID\n";			
			$p3ID = "" if $p3ID eq "##########";
			$track->{internalId} = $p3ID;
		}
	
	}

}

my $p3workspaces = $jsonh->pretty->encode($workspaces);

open OUT, ">$ARGV[1]"; 
print OUT $p3workspaces;
close OUT;


sub getP3FeatureID {

	my ($p2_feature_id) = @_;

	my $solrQ = $solr_server."/genome_feature/select?q=p2_feature_id%3A$p2_feature_id&fl=feature_id&wt=csv";

	my $result = `wget -q -O - "$solrQ"`;

	my ($p3_feature_id) = $result=~/feature_id\n(.*)/i;

	return $p3_feature_id? $p3_feature_id : "##########";

}

sub getP3GenomeID {

	my ($p2_genome_id) = @_;

	my $solrQ = $solr_server."/genome/select?q=p2_genome_id%3A$p2_genome_id&fl=genome_id&wt=csv";

	my $result = `wget -q -O - "$solrQ"`;

	my ($p3_genome_id) = $result=~/genome_id\n(.*)/;

	return $p3_genome_id? $p3_genome_id : "##########";

}
