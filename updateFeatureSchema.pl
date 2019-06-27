#!/usr/bin/env perl

# Script to migrate data for genome_feature core optimization
# Remove attribites:
# Add attribites: 
# Move to new core:  

use strict;
use warnings;
use Getopt::Long;
use JSON;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Data::Dumper;

$0 =~ m/([^\/]+)$/;
my $self = $1;
my $usage = $self. qq( --genome_list genome_list.tab --genome_id 83332.12 --help);


my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=json&rows=1000000";
#my $json = JSON->new->allow_nonref;

my ($help, $genome_list, $genome_id);
my @genome_ids = ();
my %md5_hash = ();


my $opts = GetOptions(
	"help" => \$help,
	"genome_list=s" => \$genome_list,
	"genome_id=s" => \$genome_id
);

if (!$opts || $help || !($genome_list || $genome_id)){
	warn qq(\n   usage: $usage\n\n);
	exit(0);
}

if ($genome_list){
	open LIST, $genome_list;
	@genome_ids = <LIST>;
	close LIST;
}elsif ($genome_id){
	push @genome_ids, $genome_id;
}else{

}

foreach my $genome_id (@genome_ids){

		chomp $genome_id;
		next if $genome_id=~/genome_id/;

		print "Processing genome_id: $genome_id\n";


		my $core = "/genome_feature";
		my $query = "/select?q=genome_id:$genome_id";
		my $fields = "&fl=*";
		my $sort = "";
		my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

		my $result =  `wget -q -O - "$solrQuery"`;
		my $resultObj = decode_json($result);
		
		my @features_orig = ();
		my @features = ();
		my @sequences = ();
		#my %md5_hash = ();

		foreach my $feature (@{$resultObj->{'response'}->{'docs'}}) {

			#my %feature_orig = %{$feature};
			#my $feature_orig_ref = \%feature_orig; 		
			#delete $feature_orig_ref->{_version_};
			#push @features_orig, $feature_orig_ref;
			
			if ($feature->{na_sequence}){
				my $md5 = md5_hex(lc($feature->{na_sequence}));
				$feature->{na_sequence_md5} = $md5;
				my $na_seq = { "md5" => $md5, "sequence" => lc($feature->{na_sequence}), "sequence_type" => "NA"};	
				if (!$md5_hash{$md5}){
					$md5_hash{$md5} = 1;
					push @sequences, $na_seq;
				}
			}

			if ($feature->{aa_sequence}){
				my $aa_seq = { "md5" => $feature->{aa_sequence_md5}, "sequence" => $feature->{aa_sequence}, "sequence_type" => "AA"};	
				if (!$md5_hash{$feature->{aa_sequence_md5}}){
					$md5_hash{$feature->{aa_sequence_md5}} = 1;
					push @sequences, $aa_seq;
				}
			}
			
			my @property = ();

			if ($feature->{ec}){delete $feature->{ec}; push @property, "EC number"; }
			if ($feature->{pathway}){delete $feature->{pathway}; push @property, "Pathway"; } 
			if ($feature->{subsystem}){delete $feature->{subsystem}; push @property, "Subsystem"; }

			$feature->{property} = \@property if scalar @property;

			delete $feature->{na_sequence};
			delete $feature->{aa_sequence};
			delete $feature->{uniprotkb_accession};
			delete $feature->{gi};
			delete $feature->{pos_group};
			delete $feature->{_version_};
			delete $feature->{document_type};
			delete $feature->{document_id};

			push @features, $feature;

		}
		
		#my $feature_orig_json = to_json(\@features_orig, {pretty => 1});
		#open JSON, ">$genome_id.feature_orig.json";
		#print JSON "$feature_orig_json";
		#close JSON;	

		my $feature_json = to_json(\@features, {pretty => 1});
		open JSON, ">$genome_id.feature.json";
		print JSON "$feature_json";
		close JSON;	

		my $sequence_json = to_json(\@sequences, {pretty => 1});
		open JSON, ">$genome_id.sequence.json";
		print JSON "$sequence_json";
		close JSON;

		`post.update_dev.sh genome_feature $genome_id.feature.json`;	
		`post.update_dev.sh feature_sequence $genome_id.sequence.json`;	

		#`rm $genome_id.feature_orig.json`;
		`rm $genome_id.feature.json`;
		`rm $genome_id.sequence.json`;

}
