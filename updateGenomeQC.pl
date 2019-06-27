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

my @genome_ids = ();
my @genomes = ();
my %species_stats = ();


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["genome_list=s", "File containing list of genome ids"],
    ["genome_id=s", ""],
    ["outfile=s", "Output file"],
		["commit=s", "commit, true|false", {default => "false"}],
		["clean=s", "delete index files, true | false", {default => "false"}],
		[],
		["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless ($opt->genome_list || $opt->genome_id);

if ($opt->genome_list){

	if ($opt->genome_list=~/all/){
		@genome_ids = getGenomes();
	}else{
		open LIST, $opt->genome_list or die "Can't open genome list: $opt->genome_list!";
		@genome_ids = <LIST>;
		close LIST; 
	}
}elsif($opt->genome_id){
	push @genome_ids, $opt->genome_id;
}

foreach my $genome_id (@genome_ids){

	chomp $genome_id;
	next if $genome_id=~/genome_id/;
	
	print "Processing genome: $genome_id\n";
	my $genome = {};
	
	$genome->{genome_id} = $genome_id;
	$genome->{genome_quality_flags} = [];

	getGenomeStats($genome);
	getSequenceStats($genome);
	getFeatureStats($genome);
	#getSpeciesStats($genome) if $genome->{species};

	if (scalar @{$genome->{genome_quality_flags}}){
		$genome->{genome_quality} = "Poor";
	}else{
		$genome->{genome_quality} = "Good";
	}

	my $genome_update = prepareGenomeUpdate($genome);

	push @genomes, $genome_update;
	
}

my $genome_json = $json->pretty->encode(\@genomes);
my $outfile = $opt->outfile ? $opt->outfile : "genome_update.json";
open GF, ">$outfile";
print GF "$genome_json";
close GF;

#my $taxon_json = $json->pretty->encode(\%species_stats);
#open TF, ">taxon.json";
#print TF "$taxon_json";
#close TF;

sub getGenomes{

	my @genome_ids = ();

	my $core = "/genome";
	my $query = "/select?q=public:1 AND genome_status:(Complete OR WGS)";
	my $fields = "&fl=genome_id";
	my $sort = "&sort=genus+asc";
	my $rows = "&rows=1000000"; 
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$rows.$solrFormat;
	
	push @genome_ids, `wget -q -O - "$solrQuery" | grep -v "genome_id"`;

	return @genome_ids;	

}

sub getGenomeStats {

	my ($genome) = @_;

	my $core = "/genome";
	my $query = "/select?q=genome_id:$genome->{genome_id}";
	my $fields = "&fl=genome_name,genus,species,genome_status,genome_length,contigs,patric_cds";
	my $rows = "&rows=1";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;
	
	my $result = `wget -q -O - "$solrQuery" | grep -v genome_name`;
	chomp $result;
	my ($genome_name,$genus,$species,$genome_status,$genome_length,$contigs,$patric_cds) = split /\t/, $result;

	$genome->{genome_name} = $genome_name;
	$genome->{genus} = $genus;
	$genome->{species} = $species;
	$genome->{genome_status} = $genome_status;

	#$genome->{infered_species} = ""; # check with bob

	$genome->{genome_length} = $genome_length;
	$genome->{contigs} = $contigs;
	$genome->{patric_cds} = $patric_cds;

	#$genome->{coarse_consistency} = 0;
	#$genome->{fine_consistency} = 0; # < 85
	#$genome->{checkm_completeness} = 0; # < 80
	#$genome->{checkm_contamination} = 0; # > 10

	$genome->{cds_ratio} = sprintf "%.2f", $genome->{patric_cds}*1000/$genome->{genome_length};

	# Add QC flags
	
	push @{$genome->{genome_quality_flags}}, "Plasmid only" if $genome->{genome_status} =~/plasmid/i; 
	push @{$genome->{genome_quality_flags}}, "Metagenomic bin" if $genome->{genome_status} =~/metagenome bin/i; 
	#push @{$genome->{genome_quality_flags}}, "Misidentified taxon" if $genome->{infered_taxon} && not $genome->{genome_name} =~/$genome->{predicted_species}/i;
	
	push @{$genome->{genome_quality_flags}}, "Genome too long" if $genome->{genome_length} > 15000000;
	push @{$genome->{genome_quality_flags}}, "Genome too short" if $genome->{genome_length} < 300000;

	push @{$genome->{genome_quality_flags}}, "Too many contigs" if $genome->{contigs} > 1000;
	push @{$genome->{genome_quality_flags}}, "Abnormal CDS ratio" if $genome->{cds_ratio} < 0.5 || $genome->{cds_ratio} > 1.5;

	#push @{$genome->{genome_quality_flags}}, "Low CheckM completeness score" if $genome->{checkm_completeness} && $genome->{checkm_completeness} < 80;
	#push @{$genome->{genome_quality_flags}}, "High CheckM contaminationi score" if $genome->{checkm_contamination} && $genome->{checkm_contamination} > 10;
	#push @{$genome->{genome_quality_flags}}, "Low Fine consistency score" if $genome->{fine_consistency} && $genome->{fine_consistency} < 85;


}	


sub getSequenceStats {

	my ($genome) = @_;

	my $core = "/genome_sequence";
	my $query = "/select?q=genome_id:$genome->{genome_id}";
	my $fields = "&fl=accession,length";
	my $rows = "&rows=1000000";
	my $sort = "&sort=length desc";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;
	my @results = `wget -q -O - "$solrQuery" | grep -v accession`;

	my $genome_length = 0;
	my $contigs = 0;

	foreach my $row (@results){
		chomp $row;
		my ($accession, $length) = split /\t/, $row;
		$genome_length += $length;
		$contigs++;

		$genome->{contig_l50} = $contigs if $genome_length >= $genome->{genome_length}/2 && !$genome->{contig_l50};
		$genome->{contig_n50} = $length if $genome_length >= $genome->{genome_length}/2 && !$genome->{contig_n50}; 		
	}

	push @{$genome->{genome_quality_flags}}, "High contig L50" if $genome->{contig_l50} > 500;
	push @{$genome->{genome_quality_flags}}, "Low contig N50" if $genome->{contig_n50} < 5000;

	#push @{$genome->{genome_quality_flags}}, "Inconsistent data: number of contigs" unless $contigs == $genome->{contigs};
	#push @{$genome->{genome_quality_flags}}, "Inconsistent data: genome length" unless $genome_length == $genome->{genome_length};
	
	$genome->{contigs} = $contigs;
	$genome->{genome_length} = $genome_length;

}	


sub getFeatureStats {

	my ($genome) = @_;
	
	my $core = "/genome_feature";
	my $query = "/select?q=genome_id:$genome->{genome_id} AND annotation:PATRIC";
	my $fields = "&fl=patric_id,feature_type,product,plfam_id,pgfam_id,aa_sequence_md5";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;
	my @results = `wget -q -O - "$solrQuery" | grep -v genome_name`;

	$genome->{cds} = 0;
	$genome->{trna} = 0;
	$genome->{rrna} = 0;
	#$genome->{pseudogene} = 0;
	$genome->{partial_cds} = 0;
	$genome->{hypothetical_cds} = 0;
	$genome->{plfam_cds} = 0;

	foreach my $row (@results){
		my ($patric_id,$feature_type,$product,$plfam_id,$pgfam_id,$aa_sequence_md5) = split /\t/, $row;
	
		$genome->{cds}++ if $feature_type=~/CDS/;
		$genome->{trna}++ if $feature_type=~/tRNA/;
		$genome->{rrna}++ if $feature_type=~/rRNA/;
		#$genome->{pseudogene}++ if $feature_type=~/pseudogene/;
		$genome->{partial_cds}++ if $feature_type=~/CDS/ && !$aa_sequence_md5;
		$genome->{hypothetical_cds}++ if $feature_type=~/CDS/ && $product=~/hypothetical protein/i;
		$genome->{plfam_cds}++ if $feature_type=~/CDS/ && $plfam_id=~/PLF/;
	
	}

	if ($genome->{cds} == 0){
		push @{$genome->{genome_quality_flags}}, "No CDS";
		return;	
	}

	$genome->{patric_cds} = $genome->{cds};	
	$genome->{cds_ratio} = sprintf "%.2f", $genome->{cds} * 1000 / $genome->{genome_length};
	$genome->{hypothetical_cds_ratio} = sprintf "%.2f", $genome->{hypothetical_cds} / $genome->{cds};
	$genome->{partial_cds_ratio} = sprintf "%.2f", $genome->{partial_cds} / $genome->{cds};
	$genome->{plfam_cds_ratio} = sprintf "%.2f", $genome->{plfam_cds} / $genome->{cds};

	#push @{$genome->{genome_quality_flags}}, "Inconsistent data: number of CDS" unless $genome->{cds} == $genome->{patric_cds};

	push @{$genome->{genome_quality_flags}}, "Abnormal CDS ratio" if $genome->{cds_ratio} < 0.5 || $genome->{cds_ratio} > 1.5;
	push @{$genome->{genome_quality_flags}}, "Too many hypothetical CDS" if $genome->{hypothetical_cds_ratio} > 0.7;
	push @{$genome->{genome_quality_flags}}, "Too many partial CDS" if $genome->{partial_cds_ratio} > 0.3;

}

sub getSpeciesStats {

	my ($genome) = @_;

	my $species = $genome->{species};

	if ($species_stats{$species}->{genome_count}){

	}else{

  my $core = "/genome";
  my $query = "/select?q=public:1 AND species:%22$species%22 AND genome_status:(Complete OR WGS)";
  my $fields = "";
  my $rows = "&rows=0";
  my $stat = "&stats=true&stats.field=genome_length&stats.field=patric_cds";
	my $format = "&wt=json";
  my $solrQuery = $solrServer.$core.$query.$fields.$rows.$stat.$format;
	
  my $result = `wget -q -O - "$solrQuery"`;
	return unless $result;
	my $resultObj = decode_json($result);

	$species_stats{$species}->{genome_count} = $resultObj->{stats}->{stats_fields}->{genome_length}->{count};

	$species_stats{$species}->{genome_length_mean} = $resultObj->{stats}->{stats_fields}->{genome_length}->{mean};
	$species_stats{$species}->{genome_length_sd} = $resultObj->{stats}->{stats_fields}->{genome_length}->{stddev};
	
	$species_stats{$species}->{cds_mean} = $resultObj->{stats}->{stats_fields}->{patric_cds}->{mean};
	$species_stats{$species}->{cds_sd} = $resultObj->{stats}->{stats_fields}->{patric_cds}->{stddev};
	
	#$species_stats{$species}->{hypothetical_cds_mean} = $resultObj->{stats}->{stats_fields}->{hypothetical_cds}->{mean};
	#$species_stats{$species}->{hypothetical_cds_sd} = $resultObj->{stats}->{stats_fields}->{hypothetical_cds}->{stddev};
	
	#print Dumper $result;	

	}

	return unless $species_stats{$species}->{genome_count} > 5;

	push @{$genome->{genome_quality_flags}}, "Genome too long" if 
		$genome->{genome_length} > $species_stats{$species}->{genome_length_mean} * 1.5 ||
		$genome->{genome_length} > $species_stats{$species}->{genome_length_mean} + ( 5 * $species_stats{$species}->{genome_length_sd});

	push @{$genome->{genome_quality_flags}}, "Genome too short" if 
		$genome->{genome_length} < $species_stats{$species}->{genome_length_mean} * 0.5 ||
		$genome->{genome_length} < $species_stats{$species}->{genome_length_mean} - ( 5 * $species_stats{$species}->{genome_length_sd});
	
}

sub prepareGenomeUpdate {

	my ($genome) = @_;
	my $genome_update = {};


	foreach my $key (keys %{$genome}){
		next if $key=~/^(genome_name\|genus\|species\|genome_status)$/;
		my $value = $genome->{$key};

		if ($key=~/genome_id/i){	
			$genome_update->{$key} = $value; 
		}else{
			$genome_update->{$key}->{set} = $value; 
		}	
	}

	return $genome_update;

}
