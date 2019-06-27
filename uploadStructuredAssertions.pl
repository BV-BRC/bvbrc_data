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
		["property=s", "Property: Gene name|Function|GO Term|Interaction|Essentiality|Citation"],
		["source=s", "Source: external database, project, or publication"],
		["value=s", "Value"],
		["comment=s", "Comment or assertion"],
		["evidence_code=s", "GO eviodence code"], 
		["pmid=s", "PMID: multi-value, separated by comma"], 
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
    ["public=s", "Public:true|false", { default => "false"}],
    ["owner=s", "owner", { default => "PATRIC\@patricbrc.org"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->file;

my @attribs_all = qw (feature_id patric_id refseq_locus_tag property source value comment evidence_code pmid public owner);

my $property = $opt->property;
my $source = $opt->source;
my $value = $opt->value;
my $comment = $opt->comment;
my $evidence_code = $opt->evidence_code;
my $pmid = $opt->pmid;
my $commit = $opt->commit;
my $public = ($opt->public=~/true/i) ? 1 : 0;
my $owner = $opt->owner;

open LIST, $opt->file or die "Can't open file: $opt->file!!\n";


# global arrays to hold structured assertions
my @assertions;
my @attribs;

while (my $row = <LIST>){

	chomp $row;
	my @values;

	if ($row=~/gene.id/i){
		@attribs = split /\t/, $row; 
	}else{
		@values = split /\t/, $row;
	}

	my $gene_id = $values[0] if $attribs[0]=~/gene.id/i;
	next unless $gene_id;

	print "$gene_id\n";

	# new genome, get primary identifiers and existing features
	my $core = "/genome_feature";
	my $query = "/select?q=annotation:PATRIC AND feature_type:CDS AND (patric_id:$gene_id OR refseq_locus_tag:$gene_id)";
	my $fields = "&fl=feature_id,patric_id,refseq_locus_tag";
	my $rows = "&rows=1";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$solrFormat;
	
	my $result = `wget -q -O - "$solrQuery" | grep -v "feature_id"`;
	chomp $result;

	my ($feature_id, $patric_id, $refseq_locus_tag) = split /\t/, $result;

	next unless $feature_id && $patric_id;

	# build structured assertion record
	my $assertion;
	
	$assertion->{feature_id} = $feature_id;
	$assertion->{patric_id} = $patric_id;
	$assertion->{refseq_locus_tag} = $refseq_locus_tag;

	for (my $i=1; $i< scalar @values; $i++){
		my $attrib = $attribs[$i];
		my $value = $values[$i];
		$assertion->{$attrib} = $value? $value : $opt->{$attrib};
	}

	foreach my $attrib (@attribs_all){
		next if $assertion->{$attrib};
		$assertion->{$attrib} = $opt->{$attrib};
	}

	push @assertions, $assertion;


}

close LIST;

print "\tPrepare structured_assertions.json\n";

my $assertions_json = $json->pretty->encode(\@assertions);
my $outfile = $opt->file;
$outfile=~s/\.\w*$//;
$outfile=$outfile.".json";

open OUT, ">$outfile";
print OUT "$assertions_json";
close OUT;

`post.update.sh structured_assertion $outfile` if $opt->commit eq "true";

