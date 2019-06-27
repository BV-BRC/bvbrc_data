#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use Data::Dumper;
use JSON;
use POSIX;
use Date::Parse;

my $json = JSON->new->allow_nonref;

my $public = 1;
my $owner = "PATRIC\@patricbrc.org";


# Get valid genome metadata attributes

open ATTR, "$Bin/attrib_solr";
my %solr=();
while (my $entry = <ATTR>){
	chomp $entry;
	my ($attrib, $type) = $entry=~/(.*)\t(.*)/;
	$solr{$attrib}=$type;
}
close ATTR;


# Get list of public genomes 

my $solrServer = $ENV{PATRIC_SOLR};
my $url = "$solrServer/genome/select?q=public%3Atrue AND genome_status:(Complete OR WGS)&rows=1000000&fl=genome_id%2Cgenome_name%2Ctaxon_id&wt=csv&csv.separator=%09";

`wget -q "$url" -O patric_genomes` unless -e "patric_genomes";

open GENOME, "patric_genomes";
my %genomeName=();
my %genomeID=();
my %taxon=();

while (my $entry = <GENOME>){
	chomp $entry;
	my ($genome_id, $genome_name, $taxon_id) = $entry=~/(.*)\t(.*)\t(.*)/;
	$genomeID{$genome_name}=$genome_id;
	$genomeName{$genome_id}=$genome_name;
	$taxon{$genome_id}=$taxon_id;
}
close GENOME;


# Parse metadata
my @metadata_records = ();
my @amr_records = ();
my @attribs=();
my @amr=();
my ($prev_genome_id, $metadata_record);

my $infile = $ARGV[0];
open IN, $infile or die "Can't open metadata file: $ARGV[0]!!";

my ($outfile)= $infile=~/(.*).txt/;
open LOG, "> $outfile.metadata.log";
open TAB, "> $outfile.metadata.tab";
open AMR, "> $outfile.amr.tab";


my $count=0;
while (my $line = <IN>){
 
	chomp $line;
	$count++;

	my ($genome_id, $genome_name, $taxon_id, $amr_record);
	my @values=();

	if ($count == 1) {	# parse first / header line
	
		@attribs = split (/\t/, $line);
	
	}else{	# parse all other lines

		@values = split (/\t/, $line);

		$genome_name = $values[0] if $attribs[0]=~/genome.name/i;
		$genome_id = ($attribs[0]=~/genome.id/i) ? $values[0] : $genomeID{$genome_name};
		$taxon_id = $taxon{$genome_id};
		$genome_name = $genomeName{$genome_id} unless $genome_name;

		if ($genome_id && $genome_name){
			# process the record
		}else{
			print LOG "$line\n";
			next;
		}
		
		if (scalar @amr > 0 && $genome_id ne $prev_genome_id){
			$metadata_record->{genome_id} = $prev_genome_id;
			$metadata_record->{antimicrobial_resistance}->{set} = \@amr;
			$metadata_record->{antimicrobial_resistance_evidence}->{set} = "AMR Panel";
			push @metadata_records, $metadata_record;	
			@amr=();
		}

		my $i = 0;
		foreach (@values) { 
			
			$i++;

			my $attrib = lc($attribs[$i]);
			$attrib=~s/ +/_/g;
			my $value = $values[$i];
			$attrib=~s/^\s+|\s+$//g;
			$value=~s/^\s+|\s+$//g;
			$value=~s/ *" *//g;
			
			next unless ($attrib && $value);

			next if $attrib=~/reagent|dst.media/i;
			next if $value=~/missing|not defined/i;

			if($attrib=~/critical.concentration/){
				$amr_record->{measurement} = $value;
				$amr_record->{measurement_unit} = "mg/L";
				next;
			}

			$value = ucfirst $value if $value=~/resistant|susceptible|intermediate|non-susceptible/;


			$amr_record->{genome_id} = $genome_id;
			$amr_record->{genome_name} = $genome_name;
			$amr_record->{taxon_id} = $taxon_id;
			
			$amr_record->{$attrib} = $value;

			$amr_record->{owner} = $owner;
			$amr_record->{public} = $public;

			if ($attrib=~/resistant phenotype/i && $value=~/resistant|susceptible|intermediate/i){
				push @amr, $value unless (grep {$_ eq $value} @amr); 
			}
	
			$prev_genome_id = $genome_id;

		}

		push @amr_records, $amr_record if $amr_record->{measurement} || $amr_record->{resistant_phenotype};
		
			
	}

}

if (scalar @amr > 0){
			$metadata_record->{genome_id} = $prev_genome_id;
			$metadata_record->{antimicrobial_resistance}->{set} = \@amr;
			$metadata_record->{antimicrobial_resistance_evidence}->{set} = "AMR Panel";
			push @metadata_records;	
}


close IN;
close TAB;
close LOG;

my $genome_metadata = $json->pretty->encode(\@metadata_records);
#print "$genome_metadata\n";
open JS, ">$outfile.metadata.json";
print JS "$genome_metadata";
close JS;

my $genome_amr = $json->pretty->encode(\@amr_records);
#print "$genome_amr\n";
open JS, ">$outfile.amr.json";
print JS "$genome_amr";
close JS;

