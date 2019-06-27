#!/usr/bin/env perl

###########################################################
#
# expression2solr.pl: 
#
# Script to parse Differential Expression Servivce output and
# and prepare separate JSON files for each of the Solr cores.  
#
# Input: 
#	- Diff expression output directory 
# 
# Output: Six JSON files each correspodning to a Solr core
#	- transcriptomics_experiment.json
#	- transcriptomics_sample.json
#	- transcriptomics_gene.json
#
###########################################################

use strict;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use POSIX;
use JSON;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Date::Parse;
use XML::Simple;
#use Bio::KBase::AppService::AppConfig;

use lib "$Bin";
use SolrAPI;

#my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
#my $solrh = SolrAPI->new($data_api);

my $solrh = SolrAPI->new();
my $json = JSON->new->allow_nonref;

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my ($opt, $usage) = describe_options( 
				"%c %o",
				[],
				["exp_dir=s", "Differential Expression output directory"],
				["owner=s", "owner", { default => "PATRIC\@patricbrc.org"}],
				["public=s", "Public:true|false", { default => "false"}],
				["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
				[],
				["help|h", "Print usage message and exit"] );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->exp_dir;


my $exp_dir = $opt->exp_dir;
my $owner = $opt->owner;
my $public = ($opt->public=~/true/i) ? 1 : 0;
my $commit = ($opt->commit=~/true/i) ? 1 : 0;

print "\n\nProcessing $exp_dir\n";

# Initialize global arrays to hold solr data
my $experiment;
my @samples = ();
my @genes = ();

my $eid = getEid();
my $pid = getPid();

my %pids = ();

# Process differential expression result files
processExperiment();
processSamples();
processGenes();

# write to json files
writeJson();


sub writeJson {

	print "Preparing JSON files ...\n";

	my $experiment_json = $json->pretty->encode($experiment);
	my $sample_json = $json->pretty->encode(\@samples);
	my $gene_json = $json->pretty->encode(\@genes);
	
	open FH, ">$exp_dir/transcriptomics_experiment.json" or die "Cannot write transcriptomics_experiment.json: $!"; 
	print FH "[".$experiment_json."]";
	close FH;
	`post.update.sh transcriptomics_experiment $exp_dir/transcriptomics_experiment.json` if $commit==1;

	open FH, ">$exp_dir/transcriptomics_sample.json" or die "Cannot write transcriptomics_sample.json: $!"; 
	print FH $sample_json;
	close FH;
	`post.update.sh transcriptomics_sample $exp_dir/transcriptomics_sample.json` if $commit==1;

	open FH, ">$exp_dir/transcriptomics_gene.json" or die "Cannot write transcriptomics_gene.json: $!"; 
	print FH $gene_json;
	close FH;
	`post.update.sh transcriptomics_gene $exp_dir/transcriptomics_gene.json` if $commit==1;

}

sub processExperiment {

	print "Getting experiment data ...\n";

	print "$exp_dir/experiment.json\n";

	open FH, "$exp_dir/experiment.json" or die "Can't open $exp_dir/experiment.json file!\n";
	my $expObj = decode_json(join "", <FH>);
	close FH;

	#print Dumper $expObj;
	
	$experiment->{eid} = $eid;	# how to generate eid? 
	#$experiment->{expid} = $expObj->{expid};	# ws exp id
	$experiment->{accession} = $expObj->{accession}; # add to submission form
	
	$experiment->{title} = $expObj->{title};	
	$experiment->{description} = $expObj->{desc};
	
	$experiment->{organism} = $expObj->{organism};	
	$experiment->{genome_ids} = $expObj->{genome_id};	# get from organism or id mapping? 

	$experiment->{samples} = $expObj->{samples};	
	$experiment->{genes} = $expObj->{geneMapped};
	
	$experiment->{pmid} = $expObj->{pmid};
	
	#$experiment->{owner} = $owner;	
	#$experiment->{public} = $public;

	#print Dumper $experiment;	

}


sub processSamples {

	print "Getting sample data ...\n";

	open FH, "$exp_dir/sample.json" or die "Can't open $exp_dir/sample.json file!\n";
	my $sampleObjs = decode_json(join "", <FH>);
	close FH;
	
	#print Dumper $sampleObjs;

	foreach my $sampleObj (@{$sampleObjs->{sample}}){

		my $sample;
		$sample->{eid} = $experiment->{eid};
		#$sample->{expid} = $experiment->{expid};
		$sample->{accession} = $experiment->{accession};

		$pid++;
		$sample->{pid} = $pid;
		$sample->{expname} = $sampleObj->{expname};
		
		$pids{$sampleObj->{pid}} = "$pid\t$sampleObj->{expname}";
		
		$sample->{organism} = $experiment->{organism};
		$sample->{genome_ids} = ($experiment->{genome_id});

		$sample->{expmean} = $sampleObj->{expmean};
		$sample->{expstddev} = $sampleObj->{expstddev};
		$sample->{genes} = $sampleObj->{genes};
	
		$sample->{sig_log_ratio} = $sampleObj->{sig_log_ratio};
		$sample->{sig_z_score} = $sampleObj->{sig_z_score};

		#$sample->{owner} = $owner;
  	#$sample->{public} = $public;

		push @samples, $sample;

	}

}


sub processGenes {

	print "Getting gene expression data ...\n";

	open FH, "$exp_dir/expression.json" or die "Can't open $exp_dir/expression.json file!\n";
	my $geneObjs = decode_json(join "", <FH>);
	close FH;

	#print Dumper $geneObjs;

	foreach my $geneObj (@{$geneObjs->{expression}}){

		my $gene;

		$gene->{eid} = $experiment->{eid};			
		$gene->{accession} = $experiment->{accession};
		
		my ($pid, $expname) = $pids{$geneObj->{pid}}=~/(\d+)\t(.*)/;
		$gene->{pid} = $pid;			
		$gene->{expname} = $expname;			

		$gene->{genome_id} = $2 if $geneObj->{feature_id}=~/^(PATRIC|RefSeq)\.(\d+\.\d+)\./;
		$gene->{feature_id} = $geneObj->{feature_id};
		#$gene->{patric_id} = "";
		#$gene->{refseq_locus_tag} = "";
		$gene->{alt_locus_tag} = $geneObj->{exp_locus_tag};
	
		$gene->{organism} = $experiment->{organism};
		
		#$gene->{avg_intensity} = "";
		$gene->{log_ratio} = $geneObj->{log_ratio};
		$gene->{z_score} = $geneObj->{z_score};
		
		#$gene->{owner} = $owner;
  	#$gene->{public} = $public;

		push @genes, $gene;

	}		

}


sub getEid {

	my $solrQuery = $solrServer."/transcriptomics_experiment/select?";
  $solrQuery .= "q=*:*";
  $solrQuery .= "&fl=eid&sort=eid+desc&rows=1";
  $solrQuery .= $solrFormat;
 
	my $eid = `wget -q -O - "$solrQuery" | tail -n 1`;
	$eid++;

	return $eid;

}


sub getPid {

	my $solrQuery = $solrServer."/transcriptomics_sample/select?";
  $solrQuery .= "q=*:*";
  $solrQuery .= "&fl=pid&sort=pid+desc&rows=1";
  $solrQuery .= $solrFormat;
  
	my $pid = `wget -q -O - "$solrQuery" | tail -n 1`;	
	$pid++;

	return $pid;
}
