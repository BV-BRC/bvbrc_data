#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long::Descriptive;
use JSON;
use Data::Dumper;

$0 =~ m/([^\/]+)$/;
my $self = $1;

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";
my $json = JSON->new->allow_nonref;

my @taxonomy = ();


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["min_genome_count=s", "Minimum number of genomes to compute stats", {default => 5}],
    ["outfile=s", "Output file"],
		["commit=s", "commit, true|false", {default => "false"}],
		[],
		["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;

my $taxa = getSpecies($opt->min_genome_count);

print "Total species with more than ". $opt->min_genome_count ." genomes:".scalar @{$taxa}."\n";

foreach my $species (@$taxa){

	#next unless $species->{name} eq "Erythrobacteraceae bacterium";
	
	print "Species: $species->{name}\t$species->{genome_count}\n";
	
	$species = getSpeciesTaxon($species);
	$species = getSpeciesStats($species);
	
	next unless $species->{taxon_id}; # && $species->{taxon_id} == 2026738;

	my $taxon_update = prepareTaxonUpdate($species);
	push @taxonomy, $taxon_update;

}


my $taxonomy_json = $json->pretty->encode(\@taxonomy);
my $outfile = $opt->outfile ? $opt->outfile : "update_species_stats.json";
open TF, ">$outfile";
print TF "$taxonomy_json";
close TF;


sub getSpecies{

	my ($min_genome_count) = @_;
	my @taxa = ();

	my $core = "/genome";
	my $query = "/select?q=public:1 AND genome_status:(Complete OR WGS)";
	my $facet = "&facet=on&facet.field=species&facet.mincount=$min_genome_count&facet.limit=-1";
	my $rows = "&rows=0"; 
	my $format = "&wt=json&json.nl=map";
	my $solrQuery = $solrServer.$core.$query.$facet.$rows.$format;

	my $result = `wget -q -O - "$solrQuery"`;
	return unless $result;
	my $resultObj = decode_json($result);

	my %hash = %{$resultObj->{facet_counts}->{facet_fields}->{species}};
	foreach my $key (keys %hash){
		my $species;
		$species->{name} = $key;
		$species->{genome_count} = $hash{$key};
		push @taxa, $species if $species->{name};
	}
	
	return \@taxa;	

}

sub getSpeciesTaxon	{

	my ($species) = @_;
	
	my $core = "/taxonomy";
	my $query = "/select?q=taxon_name:$species->{name}";
	my $fields = "&fl=taxon_id";
	my $rows = "&rows=1";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;
	my $taxon_id = `wget -q -O - "$solrQuery" | grep -v taxon_id`;

	chomp $taxon_id;
	$species->{taxon_id} = $taxon_id;

	return $species;
	
}

sub getSpeciesStats {

	my ($species) = @_;

	my @stat_fields = qw(genome_length patric_cds hypothetical_cds_ratio plfam_cds_ratio);

	my $core = "/genome";
	my $query = "/select?q=public:1 AND genome_status:(Complete OR WGS) AND species:\\\"$species->{name}\\\"";
	my $statsQ = "&stats=true&stats.field=".join "&stats.field=", @stat_fields;
	my $rows = "&rows=0"; 
	my $format = "&wt=json&json.nl=map";
	my $solrQuery = $solrServer.$core.$query.$statsQ.$rows.$format;
	
	my $result = `wget -q -O - "$solrQuery"`;
	return unless $result;
	my $resultObj = decode_json($result);
	my $stats = $resultObj->{stats}->{stats_fields};

	$species->{genome_length_mean} = sprintf "%0.2f", $stats->{genome_length}->{mean};
	$species->{genome_length_sd} = sprintf "%0.2f", $stats->{genome_length}->{stddev};
	
	$species->{cds_mean} = sprintf "%0.0f", $stats->{patric_cds}->{mean};
	$species->{cds_sd} = sprintf "%0.0f", $stats->{patric_cds}->{stddev};
	
	$species->{hypothetical_cds_ratio_mean} = $stats->{hypothetical_cds_ratio}->{mean} > 0? sprintf "%0.2f", $stats->{hypothetical_cds_ratio}->{mean} : 0;
	$species->{hypothetical_cds_ratio_sd} = $stats->{hypothetical_cds_ratio}->{stddev} > 0? sprintf "%0.2f", $stats->{hypothetical_cds_ratio}->{stddev} : 0;
	
	$species->{plfam_cds_ratio_mean} = $stats->{plfam_cds_ratio}->{mean} > 0? sprintf "%0.2f", $stats->{plfam_cds_ratio}->{mean} : 0;
	$species->{plfam_cds_ratio_sd} = $stats->{plfam_cds_ratio}->{stddev} > 0? sprintf "%0.2f", $stats->{plfam_cds_ratio}->{stddev} : 0;
	
	return $species;

}	


sub prepareTaxonUpdate {

	my ($taxon) = @_;
	my $taxon_update = {};

	foreach my $key (keys %{$taxon}){
		my $value = $taxon->{$key};
		next if $key=~/name/;
		if ($key=~/taxon_id/i){	
			$taxon_update->{$key} = $value; 
		}else{
			$taxon_update->{$key}->{set} = $value; 
		}	
	}

	return $taxon_update;

}
