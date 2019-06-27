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


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["genome_list=s", "File containing list of genome ids"],
    ["genome_id=s", ""],
    ["outfile=s", "Output file"],
		["mlst=s", "MLST database directory", {default => "/homes/mshukla/mlstdb/blastdb"}],
		["fasta=s", "FNA file directory", {default => "/homes/mshukla/mlstdb/blastdb"}],
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

	
	
}

my $genome_json = $json->pretty->encode(\@genomes);
my $outfile = $opt->outfile ? $opt->outfile : "genome_update.json";
open GF, ">$outfile";
print GF "$genome_json";
close GF;

sub getGenomes{

	my @genome_ids = ();

	my $core = "/genome";
	my $query = "/select?q=public:1 AND genome_status:(Complete OR WGS)";
	my $fields = "&fl=genome_id,species";
	my $sort = "&sort=genus+asc";
	my $rows = "&rows=1000000"; 
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$rows.$solrFormat;
	
	push @genome_ids, `wget -q -O - "$solrQuery" | grep -v "genome_id"`;

	return @genome_ids;	

}

