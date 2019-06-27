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
    ["infile=s", "Tab-delimited file"],
    #["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->infile;

open IN, $opt->infile or die "Can't open infile: $opt->infile!!\n";

my $outfile = $opt->infile;
$outfile=~s/\..*$/.json/;

my $count = 0;
my @attribs = ();
my @records = ();

while (my $row = <IN>) {
	chomp $row;
	$count++;

	if ($count == 1){
		@attribs = split /\t/, $row; 
	}else{
		my @values = split /\t/, $row;

		my $record;
		for (my $i=0; $i < scalar @values; $i++){
			my $attrib = lc $attribs[$i];
			$attrib =~s/  */_/g;
			
			if ($attrib=~/^(classification|antibiotics|pmid)$/i){
				@{$record->{$attrib}} = (split /,\s*/, $values[$i]);
			}else{
				$record->{$attrib} = $values[$i];
			}
	
		}
		push @records, $record;
		 
	}

}

close IN;

my $records_json = $json->pretty->encode(\@records);
open OUT, ">$outfile";
print OUT "$records_json";
close OUT;

