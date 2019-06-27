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
    ["file=s", "File containing list of pairwise interactions: gene1\tgene2"],
		["domain=s", "Domain:Bacteria|Archaea|Eukaryota", { default => "Bacteria"}],
		["interactor_type=s", "interactor_type", { default => "Protein"}],
		["category=s", "Category:PPI|HPI", { default => "PPI"}],
		["interaction_type=s", "Interaction type: genetic interaction | physical interaction", { default => "genetic interaction"}],
		["detection_method=s", "Detection method", { default => "Tn-seq experiment"}],
		["evidence=s", "Evidence", { default => "experimental"}], 
		["pmid=s", "PMID: multi-value, separated by comma"], 
    ["source_db=s", "Source DB"],
		["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->file;


open LIST, $opt->file or die "Can't open file: $opt->file!!\n";


# global arrays to hold structured assertions
my @interactions;

while (my $row = <LIST>){

	chomp $row;
	
	my ($gene1, $gene2) = $row=~/(.*)\t(.*)/;
	next unless $gene1 && $gene2;
	next if $row=~/gene/i;

	print "$gene1\t$gene2\n";

	my ($core, $query, $fields, $rows, $solrQuery, $result);
	# Get info for gene1
	$core = "/genome_feature";
	$query = "/select?q=annotation:PATRIC AND feature_type:CDS AND (patric_id:$gene1 OR refseq_locus_tag:$gene1)";
	$fields = "&fl=taxon_id,genome_id,genome_name,feature_id,patric_id,refseq_locus_tag,gene,product";
	$rows = "&rows=1";
	$solrQuery = $solrServer.$core.$query.$fields.$rows.$solrFormat;
	
	$result = `wget -q -O - "$solrQuery" | grep -v "feature_id"`;
	chomp $result;
	my ($taxon_id_a,$genome_id_a,$genome_name_a,$feature_id_a,$patric_id_a,$refseq_locus_tag_a,$gene_a,$product_a) = split /\t/, $result;


	# Get info for gene2
	$core = "/genome_feature";
	$query = "/select?q=annotation:PATRIC AND feature_type:CDS AND (patric_id:$gene2 OR refseq_locus_tag:$gene2)";
	$fields = "&fl=taxon_id,genome_id,genome_name,feature_id,patric_id,refseq_locus_tag,gene,product";
	$rows = "&rows=1";
	$solrQuery = $solrServer.$core.$query.$fields.$rows.$solrFormat;
	
	$result = `wget -q -O - "$solrQuery" | grep -v "feature_id"`;
	chomp $result;
	my ($taxon_id_b,$genome_id_b,$genome_name_b,$feature_id_b,$patric_id_b,$refseq_locus_tag_b,$gene_b,$product_b) = split /\t/, $result;

	next unless $patric_id_a && $patric_id_b;


	# build interaction record
	my $interaction;

	$interaction->{taxon_id_a} = $taxon_id_a;
	$interaction->{genome_id_a} = $genome_id_a;
	$interaction->{genome_name_a} = $genome_name_a;
	$interaction->{feature_id_a} = $feature_id_a;
	$interaction->{interactor_a} = $patric_id_a;
	$interaction->{refseq_locus_tag_a} = $refseq_locus_tag_a;
	$interaction->{gene_a} = $gene_a;
	$interaction->{interactor_desc_a} = $product_a;
	$interaction->{interactor_type_a} = $opt->interactor_type;
	$interaction->{domain_a} = $opt->domain;
	
	$interaction->{taxon_id_b} = $taxon_id_b;
	$interaction->{genome_id_b} = $genome_id_b;
	$interaction->{genome_name_b} = $genome_name_b;
	$interaction->{feature_id_b} = $feature_id_b;
	$interaction->{interactor_b} = $patric_id_b;
	$interaction->{refseq_locus_tag_b} = $refseq_locus_tag_b;
	$interaction->{gene_b} = $gene_b;
	$interaction->{interactor_desc_b} = $product_b;
	$interaction->{interactor_type_b} = $opt->interactor_type;
	$interaction->{domain_b} = $opt->domain;

	$interaction->{category} = $opt->category;
	$interaction->{interaction_type} = $opt->interaction_type;
	$interaction->{detection_method} = $opt->detection_method;
	$interaction->{evidence} = $opt->evidence;
	$interaction->{pmid} = $opt->pmid;
	$interaction->{source_db} = $opt->source_db;
	
	push @interactions, $interaction;

}

close LIST;


my $interactions_json = $json->pretty->encode(\@interactions);
my $outfile = $opt->file;
$outfile=~s/\.\w*$//;
$outfile=$outfile.".json";

print "\tPrepare output file: $outfile\n";

open OUT, ">$outfile";
print OUT "$interactions_json";
close OUT;

`post.update.sh ppi $outfile` if $opt->commit eq "true";

