#!/usr/bin/env perl

use XML::Simple;
use Data::Dumper;

my $xmlfile = $ARGV[0];

open $fh_mt, ">metadata.tab";
open $fh_dst, ">genome_dst.tab";
open $fh_amr, ">genome_amr.tab";

#my $xml = XMLin($xmlfile);
my $xml = XMLin("$xmlfile", ForceArray => ["Row"]);

#XMLout($xml, OutputFile => "$xmlfile.clean");
#exit;

foreach my $biosample_id (keys %{$xml->{BioSample}}){

	my $genome;
	my $biosample = $xml->{BioSample}->{$biosample_id};
	my $accession = $biosample->{accession};
	my $taxonomy_id = $biosample->{Description}->{Organism}->{taxonomy_id};
	my $taxonomy_name = $biosample->{Description}->{Organism}->{taxonomy_name};
	
	my $biosample_label;
	foreach my $id (@{$biosample->{Ids}->{Id}}){
			$biosample_label = $id->{content} if $id->{db_label}=~/Sample name/;
	}

	$genome->{genome_name} = $taxonomy_name;

	next unless ref $biosample->{Attributes}->{Attribute} eq 'ARRAY';

	foreach my $attrib (@{$biosample->{Attributes}->{Attribute}}){
	
		my $attrib_name = $attrib->{harmonized_name};
		my $attrib_value = $attrib->{content};

		next unless $attrib_value;

		if ($attrib_name=~/lat_lon/i){ 
			my ($latitude, $longitude) = $attrib=~/(\d+\.*\d* [NS])\s+(\d+\.*\d* [EW])/;
			$genome->{latitude} = $latitude if $latitude;
			$genome->{longitude} = $longitude if $longitude;	
		}else{
			my $patric_attrib_name = biosample2patricAttrib($attrib_name);
			next unless ($patric_attrib_name && $attrib_value && not $attrib_value=~/^ *(-|missing|na|n\/a|not available|not provided|not determined|nd|unknown) *$/i);
		
			if ($patric_attrib_name=~/other|additional/i){
				my ($attrib1, $attrib2) = $patric_attrib_name=~/^([^:]*):(.*)/;
				push @{$genome->{$attrib1}}, "$attrib2:$attrib_value";
			}elsif ($patric_attrib_name=~/comments/i){
				push @{$genome->{$patric_attrib_name}}, $attrib_value;
			}else{
				$genome->{$patric_attrib_name} = $attrib_value;
			} 
		}

	}

	$genome->{isolation_country} = $1 if $genome->{geographic_location}=~/^([^:]*):.*/;

	$genome->{strain}=~s/Primary culture//;
	$genome->{strain} = $biosample_label unless $genome->{strain};
	$genome->{strain}=~s/$taxonomy_name //;

	if ($genome->{strain} && (not $genome->{genome_name}=~/$genome->{strain}/)){
		$genome->{genome_name} .= " strain $genome->{strain}";
	}else{

	}
	#print "$accession\t$taxonomy_id\t$genome->{genome_name}\n";

	foreach my $attrib (keys %{$genome}){
		next if $attrib=~/genome_name|taxon_id/;
		my $value;
		if (ref $genome->{$attrib} eq "ARRAY"){
			$value = "[\"". join("\",\"",@{$genome->{$attrib}}). "\"]";
		}else {
			$value = $genome->{$attrib};
		}
		print "$accession\t$attrib\t$value\n";
	} 

	
	# parse AMR metadata



	
	#print "### $accession\n";	
	my $header, $rows, $format;

	if ( ref($biosample->{Description}->{Comment}->{Table}) eq "HASH" && $biosample->{Description}->{Comment}->{Table}->{class}=~/Antibiogram/i){	
		$header = $biosample->{Description}->{Comment}->{Table}->{Header};
		$rows = $biosample->{Description}->{Comment}->{Table}->{Body}->{Row};
	}elsif( ref($biosample->{Description}->{Comment}) && $biosample->{Description}->{Comment}->{class}=~/Antibiogram/i){
		$header = $biosample->{Description}->{Comment}->{Header};
		$rows = $biosample->{Description}->{Comment}->{Body}->{Row};
	}
	
	#print "###Accession=$accession\n";
	
	$fh = join(",", @{$header->{Cell}})=~/DST/ ? $fh_dst : $fh_amr;

	print $fh "#Accession";
	foreach my $cell (@{$header->{Cell}}){
		$cell = "" if ref $cell eq ref {};
		print $fh "\t$cell";
	}
	print $fh "\n";

	foreach my $row (@{$rows}){
		print $fh "$accession";
		foreach my $cell (@{$row->{Cell}}){
			$cell = "" if ref $cell eq ref {};
			print $fh "\t$cell";
		}
		print $fh "\n";
	}

}

close $fh_mt, $fh_dst, $fh_amr;


sub biosample2patricAttrib{

	my ($attribute) = @_;

	my %biosample2patric = (

		"altitude" => "altitude",
		"biomaterial_provider" => "additional_metadata:biomaterial_provider",
		"collected_by" => "additional_metadata:collected_by",
		"collection_date" => "collection_date",			
		"culture_collection" => "culture_collection", 
		"depth" => "depth", 
		"description" => "comments", 
		"env_biome" => "other_environmental:env_biome", 
		"genotype" => "other_typing:genotype", 
		"geo_loc_name" => "geographic_location", 
		"host" => "host_name", 
		"host_age" => "host_age", 
		"host_description" => "other_clinical:host_description", 
		"host_disease" => "host_health",
		"host_disease_outcome" => "other_clinical:host_disease_outcome",
		"host_disease_stage" => "other_clinical:host_disease_stage",
		"host_health_state" => "other_clinical:host_health_state",
		"host_sex" => "host_gender",
		"host_subject_id" => "other_clinical:host_subject_id", 
		"host_tissue_sampled" => "isolation_source",
		"identified_by" => "additional_metadata:identified_by",
		"isolate" => "additional_metadata:isolate", 
		"isolation_source" => "isolation_source", 
		"lab_host" => "additional_metadata:lab_host", 
		"lat_lon" => "other_environmental:lat_lon", 
		"mating_type" => "additional_metadata:mating_type", 
		"organism" => "",
		"passage_history" => "additional_metadata:passage_history",
		"pathotype" => "pathovar",
		"sample_name" => "",
		"sample_title" => "",
		"sample_type" => "additional_metadata:sample_type",
		"samp_size" => "",
		"serotype" => "serovar",
		"serovar" => "serovar",
		"specimen_voucher" => "additional_metadata:specimen_voucher", 
		"strain" => "strain",
		"isolate" => "strain",
		"subgroup" => "",
		"subtype" => "",
		"temp" => "other_environmental:temperature"
	);

	return $biosample2patric{$attribute} if $attribute;	

}
