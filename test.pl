#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use POSIX;
use JSON;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Date::Parse;
use XML::Simple;

getMetadataFromBioSample($ARGV[0]);


sub getMetadataFromBioSample {

	my($biosample_accn) = @_;

	print "Getting genome metadata from BioSample: $biosample_accn ...\n";
  
	my $xml = `wget -q -O - "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term=$biosample_accn"`;
  $xml=~s/\n//;
  my ($biosample_id) = $xml=~/<Id>(\d+)<\/Id>/;

	`wget -q -O "biosample.xml" "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=$biosample_id"`;

	return unless -f "biosample.xml";
	
	#my $xml = XMLin("biosample.xml", ForceArray => ["Row"]);
	my $xml = XMLin("biosample.xml");

	return unless ref $xml->{BioSample}->{Attributes}->{Attribute} eq 'ARRAY';

	foreach my $attrib (@{$xml->{BioSample}->{Attributes}->{Attribute}}){
	
		my $attrib_name = $attrib->{harmonized_name};
		my $attrib_value = $attrib->{content};


	# parse AMR metadata

	return unless $xml->{BioSample}->{Description}->{Comment}->{Table}->{class}=~/Antibiogram/i;

	foreach my $row (@{$xml->{BioSample}->{Description}->{Comment}->{Table}->{Body}->{Row}}){
	
		my @amr1 = @{$row->{Cell}};

		print Dumper(@amr1);

	}

	}

}

