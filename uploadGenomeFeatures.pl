#!/usr/bin/env perl 

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use JSON;
use Data::Dumper;

use lib "$Bin";
use SolrAPI;

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my $solrh = SolrAPI->new();
my $json = JSON->new->allow_nonref;


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["file=s", "File containing list of annotation files"],
		["annotation=s", "annotation: PATRIC|RefSeq|TPA"],
    ["owner=s", "owner", { default => "PATRIC\@patricbrc.org"}],
    ["public=s", "Public:true|false", { default => "false"}],
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->file;

my $public = ($opt->public=~/true/i) ? 1 : 0;

open LIST, $opt->file or die "Can't open file: $opt->file!!\n";


# global arrays to hold structured assertions
my @features;
my @attribs;

while (my $row = <LIST>){

	chomp $row;

	next if $row=~/genome.id/i;

	# genome_id	accession	source	feature	locus_tag	start	end	strand	product
	my ($genome_id,$accession,$source,$feature_type,$locus_tag,$start,$end,$strand,$product) = split /\t/, $row;

	next unless $genome_id && $accession;

	print "$row\n";

	# new genome, get primary identifiers and existing features
	my $core = "/genome_sequence";
	my $query = "/select?q=genome_id:$genome_id AND accession:$accession";
	my $fields = "&fl=genome_name,taxon_id,sequence_id,sequence";
	my $rows = "&rows=1";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$solrFormat;

	my $result = `wget -q -O - "$solrQuery" | grep -v "taxon_id"`;
	chomp $result;

	my ($genome_name, $taxon_id, $sequence_id, $sequence) = split /\t/, $result;

	next unless $result;

	print "$row\n";
	
	my $feature;
	
	$feature->{genome_id} = $genome_id;
	$feature->{genome_name} = $genome_name;
	$feature->{taxon_id} = $taxon_id;
	
	$feature->{sequence_id} = $sequence_id;
	$feature->{accession} = $accession;
	
	$feature->{annotation} = $source;
	$feature->{feature_type} = $feature_type;
	
	my $strand_wd = ($strand=~/\+/)? 'fwd':'rev';
	$feature->{feature_id} = "$feature->{annotation}.$feature->{genome_id}.$feature->{accession}."
														."$feature->{feature_type}.$start.$end.$strand_wd";

	#$feature->{patric_id} = "";
	$feature->{refseq_locus_tag} = $locus_tag;
	$feature->{alt_locus_tag} = $locus_tag;
	
	#$feature->{gene} = "";
	$feature->{product} = $product;

	$feature->{start} = $start;
	$feature->{end} = $end;
	$feature->{strand} = $strand;
	$feature->{location} = $strand=~/\+/ ? "$start..$end": "complement($start..$end)";
	$feature->{segments} = ["$start..$end"];

	$feature->{pos_group} = "$feature->{sequence_id}:$feature->{end}:+" if $feature->{strand} eq '+';
	$feature->{pos_group} = "$feature->{sequence_id}:$feature->{start}:-" if $feature->{strand} eq '-';
	
	$feature->{na_length} = $end-$start+1;
	$feature->{na_sequence} = substr($sequence, $start-1, $feature->{na_length});
	$feature->{na_sequence} =~tr/ACGTacgt/TGCAtgca/ if $feature->{strand} eq "-";
	$feature->{na_sequence} = reverse($feature->{na_sequence}) if $feature->{strand} eq "-";

	#$feature->{aa_length} = "";
	#$feature->{aa_sequence} = "";
	#$feature->{aa_sequence_md5} = "";

	$feature->{owner} = $opt->owner;
	$feature->{public} = $opt->public;

	push @features, $feature;

}

close LIST;

my $features_json = $json->pretty->encode(\@features);
my $outfile = $opt->file;
$outfile=~s/\.\w*$//;
$outfile=$outfile.".json";

open OUT, ">$outfile";
print OUT "$features_json";
close OUT;

`post.update.sh genome_feature $outfile` if $opt->commit eq "true";

