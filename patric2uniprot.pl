#!/usr/bin/perl

###############################################################################
# Script to update genome_name. Input is a tab delimited file containing 
# genome_id and new_genome_name for each genome on a separate line. The 
# script queries all relevant solr cores containing genome data and updates 
# genome names using atomic updates.    
# 
# Usage: updateGenomeName.pl new_names.tab
# 
# #############################################################################

use strict;
use warnings;
use Getopt::Long;

$0 =~ m/([^\/]+)$/;
my $self = $1;
my $usage = $self. qq( --genome_list genome_list.tab --help);

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my ($help, $genome_list, $genome_id, $outfile);
my @genome_ids = ();

my $opts = GetOptions(
    "help" => \$help,
    "genome_list=s" => \$genome_list,
    "genome_id=s" => \$genome_id,
		"outfile=s" => \$outfile
);

if (!$opts || $help || !($genome_list || $genome_id)){
  warn qq(\n   usage: $usage\n\n);
  exit(0);
}

if ($genome_list){
	open LIST, $genome_list;
	@genome_ids = <LIST>;
	close LIST; 
}

if ($genome_id){
	push @genome_ids, $genome_id;
} 

foreach my $genome_id (@genome_ids){

	chomp $genome_id;
	next if $genome_id=~/genome_id/;

	my %gis;
	my %gene_ids;
	my %protein_ids;

	my $core = "/genome_feature";
	my $query = "/select?q=genome_id:$genome_id AND annotation:PATRIC AND feature_type:CDS AND refseq_locus_tag:*";
	my $fields = "&fl=patric_id,refseq_locus_tag,gi,gene_id,protein_id";
	my $rows = "&rows=100000";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my @results = `wget -q -O - "$solrQuery" | grep -v patric_id`;

 	next unless scalar @results > 0;

	print "Genome_id:$genome_id\tFeatures:". scalar @results. "\n";


	foreach my $result (@results){
		chomp $result;
		my ($patric_id, $refseq_locus_tag, $gi, $gene_id, $protein_id) = split /\t/, $result;

		$gis{$gi} = $patric_id;
		$gene_ids{$gene_id} = $patric_id;
		$protein_ids{$protein_id} = $patric_id;

	}

	my %patric2uniprot = ();

	my $url1 = $solrServer."/id_ref/select/";
	my $query1 = "q=id_type:GeneID AND id_value:(". join (" OR ", keys %gene_ids).")";
	$query1 .= "&fl=uniprotkb_accession,id_value&rows=100000";
	$query1 .= $solrFormat;
	my @rows1 = `wget -q -O - "$url1" --post-data="$query1" | grep -v "uniprotkb_accession"`;

	foreach my $row1 (@rows1){
		my ($uniprotkb_accession, $gene_id) = $row1=~/(.*)\t(.*)/;
		my $patric_id = $gene_ids{$gene_id};
		$patric2uniprot{$patric_id} = $uniprotkb_accession;
		print "$uniprotkb_accession\t$patric_id\t$gene_id\n";
	}
	
	my $url2 = $solrServer."/id_ref/select/";
	my $query2 = "q=id_type:RefSeq AND id_value:(". join (" OR ", keys %protein_ids).")";
	$query2 .= "&fl=uniprotkb_accession,id_value&rows=100000";
	$query2 .= $solrFormat;
	my @rows2 = `wget -q -O - "$url2" --post-data="$query2" | grep -v "uniprotkb_accession"`;

	foreach my $row2 (@rows2){
		my ($uniprotkb_accession, $protein_id) = $row2=~/(.*)\t(.*)/;
		my $patric_id = $protein_ids{$protein_id};
		$patric2uniprot{$patric_id} = $uniprotkb_accession;
		print "$uniprotkb_accession\t$patric_id\t$protein_id\n";
	}

	my $url3 = $solrServer."/id_ref/select/";
	my $query3 = "q=id_type:GI AND id_value:(". join (" OR ", keys %gis).")";
	$query3 .= "&fl=uniprotkb_accession,id_value&rows=100000";
	$query3 .= $solrFormat;
	my @rows3 = `wget -q -O - "$url3" --post-data="$query3" | grep -v "uniprotkb_accession"`;

	foreach my $row3 (@rows3){
		my ($uniprotkb_accession, $gi) = $row3=~/(.*)\t(.*)/;
		my $patric_id = $gis{$gi};
		$patric2uniprot{$patric_id} = $uniprotkb_accession;
		print "$uniprotkb_accession\t$patric_id\t$gi\n";
	}

	next unless scalar (keys %patric2uniprot) > 0;

	open MAP, ">$genome_id.map"; 
	
	foreach my $patric_id (keys %patric2uniprot){
		print MAP "$patric2uniprot{$patric_id}\t$patric_id\n";
	}
	
	close MAP;

}
