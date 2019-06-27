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
    ["genome_list=s", "File containing list of annotation files"],
		["property=s", "Property: Antibiotic Resistance|Essential Gene"],
		["source=s", "Source database, PATRIC if generated using prediction tool"],
		["evidence=s", "Evidence: Literature|BLASTP|K-mer Search|FBA"], 
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->genome_list && $opt->property && $opt->source && $opt->evidence;

open LIST, $opt->genome_list or die "Can't open genome_list file: $opt->genome_list!!\n";

while (my $file_name = <LIST>) {
	chomp $file_name;

	my $genome_id = $file_name;
	$genome_id =~s/.*\///;
	$genome_id=~s/.essentials//;
	print "Processing $genome_id\n";

	# process each genome files
	if (-f $file_name){ 
		open GENOME, "$file_name" or next;
	}else {
		print "\t$file_name doesn't exist!!\n";
	}  

	# global arrays to hold spgene records
	my @spgenes;

	# new genome, get primary identifiers and existing features
	my $features = getFeatures($genome_id);

	while (my $patric_id = <GENOME>){
		chomp $patric_id;
		next unless $patric_id;
		push @spgenes, prepareSPGenes($features->{$patric_id}) if $features->{$patric_id};
  }
	close GENOME;

	prepareJsonFile($genome_id, \@spgenes);
	postJsonFile($genome_id) if $opt->commit eq "true";

}

close LIST;



sub prepareJsonFile {

	my($genome_id, $spgenes) = @_;

	print "\tPrepare $genome_id.sp_gene.json\n";
	my $spgene_json = $json->pretty->encode($spgenes);
	open SPG, ">$genome_id.sp_gene.json";
	print SPG "$spgene_json";
	close SPG;

}


sub postJsonFile {

	my($genome_id) = @_;

	print "\tPost $genome_id.sp_gene.json\n";
	`post.update.sh sp_gene $genome_id.sp_gene.json`;
	#`rm $genome_id.genome_feature.json`;
	
}


sub getFeatures {

	my ($genome_id) = @_;

	my $core = "/genome_feature";
	my $query = "/select?q=annotation:PATRIC AND feature_type:CDS AND genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,taxon_id,sequence_id,accession,annotation,feature_id,patric_id,alt_locus_tag,refseq_locus_tag,gene,product,owner,public";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrFormat="&wt=json";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my $result = `wget -q -O - "$solrQuery"`;

	my $resultObj = decode_json($result);
	my $features;

	foreach my $record(@{$resultObj->{'response'}->{'docs'}}){
		$features->{$record->{patric_id}} = $record;

	}

	return $features;

}

sub prepareSPGenes {

	my ($feature) = @_;

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

	$spgene->{property} = $opt->property;
	$spgene->{source} = $opt->source;
	$spgene->{property_source} = $opt->property.': '.$opt->source;
	$spgene->{evidence} = $opt->evidence;
		
	$spgene->{source_id} = "";
	$spgene->{organism} = "";
	$spgene->{function} = "Essential gene, predicted using metabolic flux-balance analysis in LB media" if $opt->property=~/Essential Gene/ && $opt->evidence eq "FBA";
	$spgene->{classification} = ""; 

	$spgene->{owner} = $feature->{owner};
	$spgene->{public} = $feature->{public};

	return $spgene;

}
