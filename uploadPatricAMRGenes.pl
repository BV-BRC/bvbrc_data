#!/usr/bin/env perl 

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use JSON;
use Data::Dumper;
use URI::Escape;

use lib "$Bin";
use SolrAPI;

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my $solrh = SolrAPI->new();
my $json = JSON->new->allow_nonref;


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["amr_roles=s", "File containing list of curated AMR roles and related metadata"],
    ["json_file=s", "Output JSON files for indexing new sp genes in Solr"],
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->amr_roles && $opt->json_file;

my $json_file = $opt->json_file;
my $property = "Antibiotic Resistance";
my $source = "PATRIC";
my $evidence = "K-mer Search";
	
# global arrays to hold spgene records
my @spgenes;

open LIST, $opt->amr_roles or die "Can't open AMR roles file: $opt->amr_roles!!\n";

while (my $row = <LIST>) {
	chomp $row;

	#my ($amr_role, $role_abbr, $pubmed_ids, $classification) = split /\t/, $row;
	my ($amr_role, $classification, $antibiotics_class, $antibiotics, $pubmed_ids) = split /\t/, $row;

	next if $amr_role=~/functional role/i;

	$pubmed_ids=~s/[,;] */;/g;
	$pubmed_ids=~s/(conjecture|uncharacterized|not found)(,*)//gi;
	$classification=~s/_/ /g;

	#$antibiotics = lc ($antibiotics_class.",".$antibiotics);
	$antibiotics = lc ($antibiotics);
	$antibiotics=~s/_/ /g;
	$antibiotics=~s/[,;] */;/g;
	
	$antibiotics_class= lc ($antibiotics_class);
	$antibiotics_class=~s/_/ /g;
	$antibiotics_class=~s/, */,/g;

	print "Processing AMR Role: $amr_role";

	# new genome, get primary identifiers and existing features
	my $features = getFeatures($amr_role);

	print "\tFeatures found = ". scalar @{$features}. "\n";

	foreach my $feature (@{$features} ){
		push @spgenes, prepareSPGenes($feature, $pubmed_ids, $classification, $antibiotics_class, $antibiotics);
  }

}

close LIST;

prepareJsonFile($json_file, \@spgenes);
postJsonFile($json_file) if $opt->commit eq "true";


sub prepareJsonFile {

	my($json_file, $spgenes) = @_;

	print "\tPrepare JSON file $json_file\n";
	my $spgene_json = $json->pretty->encode($spgenes);
	open SPG, ">$json_file";
	print SPG "$spgene_json";
	close SPG;

}


sub postJsonFile {

	my($json_file) = @_;

	print "\tPost $json_file\n";
	`post.update.sh sp_gene $json_file`;
	#`rm $genome_id.genome_feature.json`;
	
}


sub getFeatures {

	my ($amr_role) = @_;
	my @features = ();

	my $amr_role_encoded = uri_escape($amr_role);
	
	my $core = "/genome_feature";
	my $query = "/select?q=public:1 AND annotation:PATRIC AND feature_type:CDS AND product:\"$amr_role_encoded\"";
	my $fields = "&fl=genome_id,genome_name,taxon_id,sequence_id,accession,annotation,".
		"feature_id,patric_id,alt_locus_tag,refseq_locus_tag,gene,product,owner,public";
	my $rows = "&rows=10000000";
	my $sort = "";
	my $solrFormat="&wt=json";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;
	$solrQuery=~s/"/\\"/g;	

	my $result = `wget -q -O - "$solrQuery"`;

	return \@features unless $result;

	my $resultObj = decode_json($result);

	foreach my $feature (@{$resultObj->{'response'}->{'docs'}}){
		next unless $feature->{product} eq $amr_role;	
		push @features, $feature;
	}

	return \@features;

}

sub prepareSPGenes {

	my ($feature, $pubmed_ids, $classification, $antibiotics_class, $antibiotics) = @_;

	my $spgene;
	
	$spgene->{genome_id} = $feature->{genome_id};	
	$spgene->{genome_name} = $feature->{genome_name};	
	$spgene->{taxon_id} = $feature->{taxon_id};	
		
	$spgene->{feature_id} = $feature->{feature_id};
	$spgene->{patric_id} = $feature->{patric_id};	
	$spgene->{alt_locus_tag} = $feature->{alt_locus_tag};
	$spgene->{refseq_locus_tag} = $feature->{refseq_locus_tag};
		
	$spgene->{gene} = $feature->{gene};
	$spgene->{product} = $feature->{product};

	$spgene->{property} = $property;
	$spgene->{source} = $source;
	$spgene->{property_source} = $property.': '.$source;
	$spgene->{evidence} = $evidence;
		
	$spgene->{source_id} = "";
	$spgene->{organism} = "";
	$spgene->{classification} = $classification if $classification; 
	#$spgene->{antibiotics_class} = [split /,/, $antibiotics_class] if $antibiotics_class;
	$spgene->{antibiotics_class} = $antibiotics_class if $antibiotics_class;
	$spgene->{antibiotics} = [split /[,;]/, $antibiotics] if $antibiotics;
	$spgene->{pmid} = $pubmed_ids if $pubmed_ids;

	$spgene->{owner} = $feature->{owner};
	$spgene->{public} = $feature->{public};

	return $spgene;

}
