#!/usr/bin/env perl

use strict; 
use warnings;
use FindBin qw($Bin);
use Getopt::Long::Descriptive; 
use Data::Dumper; 

my @annotations = ('PATRIC', 'RefSeq');
my %uid_names = (		'genome' => 'genome_id', 
										'genome_sequence' => 'sequence_id',
										'genome_feature' => 'feature_id', 
										'pathway' => 'id',
										'sp_gene' => 'id',
										'genome_amr' => 'id'
									);

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my ($opt, $usage) =  
  describe_options( 
    "%c %o", 
    ["core=s", "Solr core name"], 
    ["attribute=s", "Solr attribute name"], 
    ["annotation=s", "Annotation: PATRIC | RefSeq", {default => "PATRIC"}], 
    ["owner=s", "Owner: All | PATRIC", { default => "PATRIC"}], 
    ["public=s", "Pubic: true | false", { default => "true"}], 
    ["replace=s", "Tab-delimited file containing current and new attribute values to be replaced"], 
    [], 
    ["help|h", "Print usage message and exit"] ); 
 
print($usage->text), exit 0 if $opt->help; 
print($usage->text), exit 1 unless ($opt->core && $opt->attribute); 
print($usage->text), exit 1 unless $opt->replace; 

my $core = $opt->core;
my $attribute = $opt->attribute;
my $annotation = $opt->annotation;
my $owner = $opt->owner;
my $public = $opt->public;
my $replace = $opt->replace;
my $uid_name = $uid_names{$core};

my $jsonFile = $core.".json";
open JSON, ">$jsonFile" or die "Can't open json file for writing: $jsonFile";

open IN, $replace or die "Can't open replace attribute file: $replace\n!";

my $count = 0;
print JSON "[";

foreach my $entry (<IN>){

	my ($value, $new_value) = $entry=~/(.*)\t(.*)/;

	my $core = "/$core";
	my $query = "/select?q=$attribute:\\\"$value\\\"";
	$query .= " AND annotation:$annotation" unless $core=~/(genome|genome_amr|genome_sequence)$/i ;
	$query .= " AND owner:(PATRIC OR PATRIC\@patricbrc.org)" if $owner=~/PATRIC/i;
	$query .= " AND public:$public" if $public;
	my $fields = "&fl=$uid_name,$attribute";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	#print "\n$solrQuery\n";

	my @results = `wget -q -O - "$solrQuery"`;

	next unless scalar @results > 1;

	my $no_of_results = (scalar @results) - 1;

	#print "$no_of_results\t$value => $new_value\n";


	foreach my $result (@results) {
			
		my ($uid, $attribute_value) = $result=~/(.*)\t(.*)/;
		next if $uid=~/$uid_name/;
		next if $attribute_value eq $new_value;
		next unless $attribute_value eq $value || ucfirst $attribute_value eq $value;

		print "$uid\t$attribute_value\t$value\t$new_value\n";
		

		$count++;
		print JSON ",\n" unless $count == 1;
		print JSON "{\"$uid_name\":\"$uid\",\"$attribute\":{\"set\":\"$new_value\"}}" unless $new_value=~/^null$|^\s*$/i;
		print JSON "{\"$uid_name\":\"$uid\",\"$attribute\":{\"set\":null}}" if $new_value=~/^null$|^\s*$/;
		
	}

}
	
print JSON "]";
close JSON;

#`post.update.dayhoff.sh $core_name $jsonFile`;

#`rm $jsonFile`;
	

