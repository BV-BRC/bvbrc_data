#!/usr/bin/env perl

use strict; 
use warnings;
use FindBin qw($Bin);
use Getopt::Long::Descriptive; 
use Data::Dumper; 
use JSON;

my %uid_names = (		'genome' => 'genome_id', 
										'genome_sequence' => 'sequence_id',
										'genome_feature' => 'feature_id', 
										'pathway' => 'id',
										'subsystem' => 'id',
										'sp_gene' => 'id',
										'genome_amr' => 'id',
										'feature_sequence' => 'md5'
									);

#my $solrServer = $ENV{PATRIC_SOLR_DEV};
my $solrServer = $ENV{PATRIC_SOLR};

my ($opt, $usage) =  
  describe_options( 
    "%c %o", 
    ["core=s", "Solr core name"], 
    [], 
    ["help|h", "Print usage message and exit"] ); 
 
print($usage->text), exit 0 if $opt->help; 
print($usage->text), exit 1 unless ($opt->core); 

my $core = $opt->core;
my $uid_name = $uid_names{$core}? $uid_names{$core} : "id";
my $cursor = "*";
my $next_cursor = "";
my $rows = 100000;
my $counter = 0;
my $format="&wt=json&indent=on";


while ($cursor ne $next_cursor){
	
	$counter++;

	$cursor = $next_cursor unless $counter == 1;

	print "Processing $core\t$counter\t$cursor\n";

	my $outfile = $core."_".$counter.".json";

	my $core = "/$core";
	my $query = "/select?q=date_inserted:[2019-04-15T00:00:00Z TO 2019-04-30T00:00:00Z]";
	#my $query = "/select?q=*:*";
	my $fields = "&fl=*";
	my $sort = "&sort=$uid_name+asc";
	my $rows = "&rows=$rows";
	my $cursorMark = "&cursorMark=$cursor";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$rows.$format.$cursorMark;

	print "\t$solrQuery\n";

	#`wget -q -O $outfile "$solrQuery"`;

	#$next_cursor = `grep nextCursor $outfile`;
	#$next_cursor=~s/^\s*\"nextCursorMark\"\s*:\s*\"|\"\s*}\s*$//g;
	#$next_cursor=~s/\//%2F/g;

	my $result = `wget -q -O - "$solrQuery"`;
	my $resultObj = decode_json($result);

	$next_cursor = $resultObj->{nextCursorMark};
	$next_cursor=~s/\//%2F/g;
	
	foreach my $doc (@{$resultObj->{response}->{docs}}){
		delete $doc->{_version_};
	}

	my $docsJson = to_json(\@{$resultObj->{response}->{docs}}, {pretty => 1});
	open FH, ">$outfile";
	print FH "$docsJson";
	close FH;

}
	

