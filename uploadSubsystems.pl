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

my $solrh = SolrAPI->new($ENV{PATRIC_DATA_API}, $ENV{PATRIC_REFERENCE_DATA});
my $json = JSON->new->allow_nonref;


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["genome_list=s", "File containing list of annotation files"],
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->genome_list;

open LIST, $opt->genome_list or die "Can't open genome_list file: $opt->genome_list!!\n";

while (my $file_name = <LIST>) {
	chomp $file_name;

	my $genome_id = $file_name;
	$genome_id =~s/.*\///;
	$genome_id=~s/.tbl//;
	print "Processing $genome_id\n";

	# process each genome files
	if (-f $file_name){ 
		open GENOME, "$file_name" or next;
	}else {
		print "\t$file_name doesn't exist!!\n";
	}  

	# global arrays to hold spgene records
	my @subsystems;

	# new genome, get primary identifiers and existing features
	my $features = getFeatures($genome_id);
	my $update_features;

	while (my $subsystem = <GENOME>){
		chomp $subsystem;
		my ($genome_id, $patric_id, $ss_id, $ss_name, $ss_superclass, $ss_class, $ss_subclass, $role_id, $role_name, $active) = split /\t/, $subsystem;

		push @subsystems, prepareSubsystems($features->{$patric_id}, $subsystem) if $features->{$patric_id};
  	push @{$update_features->{$features->{$patric_id}->{feature_id}}->{set}}, $ss_name 
			unless (grep {$_ eq $ss_name} @{$update_features->{$features->{$patric_id}->{feature_id}}->{set}}); 
	}
	close GENOME;

	prepareJsonFile($genome_id, \@subsystems, $update_features);
	postJsonFile($genome_id) if $opt->commit eq "true";

}

close LIST;



sub prepareJsonFile {

	my($genome_id, $subsystems, $update_features) = @_;

	print "\tPrepare $genome_id.genome_feature.json\n";
	my $feature_json = $json->pretty->encode($update_features);
	open GF, ">$genome_id.genome_features.json";
	print GF "$feature_json";
	close GF;
	
	print "\tPrepare $genome_id.subsystem.json\n";
	my $subsystem_json = $json->pretty->encode($subsystems);
	open SS, ">$genome_id.subsystem.json";
	print SS "$subsystem_json";
	close SS;

}


sub postJsonFile {

	my($genome_id) = @_;

	print "\tPost $genome_id.subsystem.json\n";
	`post.update.sh subsystem $genome_id.subsystem.json`;
	#`rm $genome_id.genome_feature.json`;
	
}


sub getFeatures {

	my ($genome_id) = @_;

	my $core = "/genome_feature";
	my $query = "/select?q=annotation:PATRIC AND feature_type:CDS AND genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,taxon_id,feature_id,patric_id,refseq_locus_tag,gene,product,owner,public";
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

sub prepareSubsystems {

	my ($feature,$ss) = @_;

	# Genome, Feature, SubsystemId, SubsystemName, Classification 1, Classification 2, RoleId, RoleName, ActiveorLikely
	my ($genome_id, $patric_id, $ss_id, $ss_name, $ss_superclass, $ss_class, $ss_subclass, $role_id, $role_name, $active) = split /\t/, $ss;

	my $subsystem;
	
	$subsystem->{genome_id} = $feature->{genome_id};	
	$subsystem->{genome_name} = $feature->{genome_name};	
	$subsystem->{taxon_id} = $feature->{taxon_id};	
		
	$subsystem->{feature_id} = $feature->{feature_id};
	$subsystem->{patric_id} = $feature->{patric_id};	
	$subsystem->{refseq_locus_tag} = $feature->{refseq_locus_tag};
	#$subsystem->{alt_locus_tag} = $feature->{alt_locus_tag};
		
	$subsystem->{gene} = $feature->{gene};
	$subsystem->{product} = $feature->{product};

	$subsystem->{role_id} = $role_id;
	$subsystem->{role_name} = $role_name;

	$subsystem->{subsystem_id} = $ss_id;
	$subsystem->{subsystem_name} = $ss_name;
	$subsystem->{superclass} = $ss_superclass;
	$subsystem->{class} = $ss_class;
	$subsystem->{subclass} = $ss_subclass;
	$subsystem->{active} = $active;

	$subsystem->{owner} = $feature->{owner};
	$subsystem->{public} = $feature->{public};

	return $subsystem;

}
