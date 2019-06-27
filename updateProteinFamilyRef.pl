#!/usr/bin/env perl 

#use strict;
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

my %plfams = ();
my %pgfams = ();

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


# process functions and prepare update files
processAnnotationFiles($opt->genome_list);

open PGF, ">pgfam.ref" or die "Can't open pgfam.ref for writing!\n";
open PLF, ">plfam.ref" or die "Can't open plfam.ref for writing!\n";

sub processAnnotationFiles {

 my ($genome_list) = @_;
 open LIST, "$genome_list" or die "Can't open genome_list file: $genome_list!!\n";

 while (my $file_name = <LIST>) {

	chomp $file_name;
 	my $genome_id = $file_name;
	$genome_id =~s/.*\///;		
	print "Processing $genome_id\n";

	# process each genome files
	if (-f $file_name){ 
		open GENOME, "$file_name" or next;
	}else {
		print "\t$file_name doesn't exist!!\n";
	}  

	while (my $row = <GENOME>){

		chomp $row;

		my ($patric_id, $update_code, $old_product, $old_pgfam_id, $old_plfam_id, 
				$product, $pgfam_id, $plfam_id) = split /\t/, $row;


		my ($plfam, $pgfam);

		if ($plfam_id){	
			my ($plfam_prefix, $plfam_idx) = $plfam_id=~/(PLF_\d+)_(\d+)/;
			$plfam = $plfam_prefix."_".sprintf("%08d", $plfam_idx);
		}

		$pgfam = $pgfam_id if $pgfam_id; 
		
		#print "$patric_id\t$pgfam\t$plfam\t$product\n";

		if ($pgfam && $product){
			#$pgfams{$pgfam}{$product} = 0 unless $pgfams{$pgfam}{$product};
			$pgfams{$pgfam}{$product}++;
			#print "$pgfam\t$product\t$pgfams{$pgfam}{$product}\n";
		}

		if ($plfam && $product){
			#$plfams{$plfam}{$product} = 0 unless $plfams{$plfam}{$product};
			$plfams{$plfam}{$product}++;
			#print "$plfam\t$product\t$plfams{$plfam}{$product}\n";
		}

  }

	close GENOME;

 }

  close LIST;

}


foreach my $pgfam (keys %pgfams){
	my $max = 0;
	foreach my $product (keys %{$pgfams{$pgfam}}){
		next unless $product && $pgfams{$pgfam}{$product};
		print "$pgfam\t$product\t$pgfams{$pgfam}{$product}\n";
		 if ($pgfams{$pgfam}{$product} > $max){
			$pgfams{$pgfam} = $product;
			$max = $pgfams{$pgfam}{$product} if $pgfams{$pgfam}{$product};
		}
	}	
}

foreach my $pgfam (keys %pgfams){
	print PGF "$pgfam\t$pgfams{$pgfam}\n";
}


foreach my $plfam (keys %plfams){
	my $max = 0;
	foreach my $product (keys %{$plfams{$plfam}}){
		next unless $product && $plfams{$plfam}{$product};	
		print "$plfam\t$product\t$plfams{$plfam}{$product}\n";
		if ($plfams{$plfam}{$product} > $max){
			$plfams{$plfam} = $product;
			$max = $plfams{$plfam}{$product} if $plfams{$plfam}{$product};
		}
	}	
}

foreach my $plfam (keys %plfams){
	print PLF "$plfam\t$plfams{$plfam}\n";
}


close PGF;
close PLF;

