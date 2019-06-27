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

my $infile = $ARGV[0];
open IN, $infile or die "Can't open metadata file: $ARGV[0]!!";

my ($outfile)= $infile=~/(.*).txt/;
open LOG, "> $outfile.metadata.log";
open TAB, "> $outfile.metadata.tab";
open AMR, "> $outfile.amr.tab";


my @attribs=();
my $count=0;
while (my $line = <IN>){
 
	chomp $line;
	$count++;

	my ($genome_id, $genome_name, $taxon_id, $metadata_record, $amr_info, $amr_record, $prev_attrib, $pmid);
	my @amr=();
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

		my $i = 0;
		foreach (@values) { 
			
			$i++;

			my $attrib = lc($attribs[$i]);

			my $value = $values[$i];
			$attrib=~s/^\s+|\s+$//g;
			$value=~s/^\s+|\s+$//g;
			$value=~s/ *" *//g;
			
			next unless ($attrib && $value);

			$pmid = $value if $attrib=~/publication|pmid/i; 	

			if ($attrib=~/^(additional.typing|other.typing|other.metadata|clinical|environmental|other.clinical|other.environmental|additional.metadata) *:/i){
				my ($attrib1,$attrib2) = $attrib=~/(.*) *: *(.*)/;
				$attrib1=~s/ +/_/g;
				$attrib=$attrib1;
				$value="$attrib2:$value";

				$attrib="other_typing" if $attrib=~/typing|additional.typing|other.typing/i;
				$attrib="other_clinical" if $attrib=~/^(clinical|other.clinical)/i;
				$attrib="other_environmental" if $attrib=~/^(environmental|other.environmental)/i;
				$attrib="additional_metadata" if $attrib=~/^(other.metadata|additional.metadata)/i;
			}elsif($attrib=~/amr/){
				$attrib="antimicrobial_resistance" if $attrib eq "amr";
				$attrib="antimicrobial_resistance_evidence" if $attrib=~/amr.evidence/;	
			}else{
				$attrib=~s/ *: */:/g;
				$attrib=~s/ +/_/g;
				$attrib="host_name" if $attrib eq "host";
			}


			if ($genome_id && $solr{$attrib}) {			
		
				my ($action, $value_json);

				$value = strftime "%Y-%m-%dT%H:%M:%SZ", localtime str2time($value) if $attrib=~/completion_date/;

				print TAB "$genome_id\t$genome_name\t$attrib\t$action\t$value\n";
			
				$action = $solr{$attrib} eq "multivalue"? "add" : "set";

				$value_json = "\"$value\"" if $action eq "set";
				$value_json = "[\"$value\"]" if $action eq "add";
				$value_json =~s/, */","/g if $action eq "add" && $attrib eq "antimicrobial_resistance";

				$metadata_record->{genome_id} = $genome_id;
				$metadata_record->{$attrib}->{$action} = $value;

			}elsif($genome_id && !$solr{$attrib}){ # Not a defined genome metadata attribute, must be AMR  

				my ($antibiotic, $mic, $mic_sign, $mic_value, $sir, $method);

				# AMR Panel metadata 
				if($attrib=~/source|laboratory.typing.method|laboratory.typing.method.version|laboratory.typing.platform|vendor|testing.standard|testing.standard.year/i){

					$amr_info->{laboratory_typing_method} = $value if $attrib=~/laboratory.typing.method\s*$/i;	
					$amr_info->{laboratory_typing_method_version} = $value if $attrib=~/laboratory.typing.method.version/i;	
					$amr_info->{laboratory_typing_platform} = $value if $attrib=~/laboratory.typing.platform/i;	
					$amr_info->{vendor} = $value if $attrib=~/vendor/i;	
					$amr_info->{testing_standard} = $value if $attrib=~/testing.standard\s*$/i;	
					$amr_info->{testing_standard_year} = $value if $attrib=~/testing.standard.year\s*$/i;	
					$amr_info->{source} = $value if $attrib=~/source\s*/i;	

					#print "$genome_id\t$attrib\t$value\n";

				}else{ # AMR Panel Data, mic/sir

					#print "$genome_name\t$attrib\t$value\n";

					# check if the attrib is the same as prev atrib, must be SIR interpretation for prev MIC value
					if ($attrib eq $prev_attrib){

						$sir=$value if $value=~/^(resistant|r|susceptible|s|sensitive|intermediate|i|non-susceptible|ns)$/i;
		
						$sir="Resistant" if $value=~/^(resistant|r)$/i;
						$sir="Susceptible" if $value=~/^(susceptible|s|sensitive)$/i;
						$sir="Intermediate" if $value=~/^(intermediate|i)$/i;
						$sir="Non-susceptible" if $value=~/^(non-susceptible|ns|nonsusceptible)$/i;
	
						push @amr, "Resistant" if $sir=~/Resistant/ && !(grep {$_ eq "Resistant"} @amr); 
						push @amr, "Susceptible" if $sir=~/Susceptible/ && !(grep {$_ eq "Susceptible"} @amr); 
						push @amr, "Intermediate" if $sir=~/Intermediate/ && !(grep {$_ eq "Intermediate"} @amr); 
						push @amr, "Non-susceptible" if $sir=~/Non-susceptible/ && !(grep {$_ eq "Non-susceptible"} @amr); 

						$amr_record->{resistant_phenotype} = $sir? $sir: "";

						push @amr_records, $amr_record if $amr_record->{measurement} || $amr_record->{resistant_phenotype};
						print AMR "$genome_id\t$genome_name\t$taxon_id\t$antibiotic\t$mic\t$mic_sign\t$mic_value\t$sir\t$method\n";
			
						$prev_attrib = $attrib;
						$amr_record = {};
						next;
					}else{ # next record, push the prev record first
						push @amr_records, $amr_record if $amr_record->{measurement} || $amr_record->{resistant_phenotype};
						print AMR "$genome_id\t$genome_name\t$taxon_id\t$antibiotic\t$mic\t$mic_sign\t$mic_value\t$sir\t$method\n";	
					}

					$amr_record = {};

					if ($attrib=~/([^:]*)\s*:\s*(.*)/){
						$amr_record->{antibiotic} = lc($1);
						$amr_record->{laboratory_typing_method} = $2;
						if ($amr_record->{laboratory_typing_method}=~/([^:]*)\s*:\s*(.*)/){
							$amr_record->{laboratory_typing_method}=$1;
							$amr_record->{laboratory_typing_method_version}=$2;
						}
						$amr_record->{laboratory_typing_method}=~s/_/ /g;
						$amr_record->{laboratory_typing_method_version}=~s/_/ /g;
					}else{
						$amr_record->{antibiotic} = $attrib;
						$antibiotic = $attrib;
					}		
				
					#print "###\t$genome_id\t$attrib\t$antibiotic\t$amr_record->{antibiotic}\t$value\n";

					if ($value=~/^(resistant|r|susceptible|s|sensitive|intermediate|i|non-susceptible|ns)$/i){
						$sir=$value;
					}else{
						($mic_sign, $mic_value, $sir) = $value=~/([>=<]*)\s*([\d\/\.]*)\s*[\(\[]*([SIRsir]*)[\)\]]*/;
      		}
		
					$mic = $mic_sign.$mic_value;
					$sir="Resistant" if $value=~/^(resistant|r)$/i;
					$sir="Susceptible" if $value=~/^(susceptible|s|sensitive)$/i;
					$sir="Intermediate" if $value=~/^(intermediate|i)$/i;
					$sir="Non-susceptible" if $value=~/^(non-susceptible|ns|nonsusceptible)$/i;
	
					push @amr, "Resistant" if $sir=~/Resistant/ && !(grep {$_ eq "Resistant"} @amr); 
					push @amr, "Susceptible" if $sir=~/Susceptible/ && !(grep {$_ eq "Susceptible"} @amr); 
					push @amr, "Intermediate" if $sir=~/Intermediate/ && !(grep {$_ eq "Intermediate"} @amr); 
					push @amr, "Non-susceptible" if $sir=~/Non-susceptible/ && !(grep {$_ eq "Non-susceptible"} @amr); 
			

					$amr_record->{genome_id} = $genome_id;
					$amr_record->{genome_name} = $genome_name;
					$amr_record->{taxon_id} = $taxon_id;
					$amr_record->{measurement} = $mic? $mic : "";
					$amr_record->{measurement_sign} = $mic_sign? $mic_sign: "";
					$amr_record->{measurement_value} = $mic_value? $mic_value: "";
					$amr_record->{measurement_unit} = "mg/L" if $mic_value && !($amr_record->{laboratory_typing_method}=~/(disk|disc).diffusion/i);
					$amr_record->{measurement_unit} = "mm" if $mic_value && $amr_record->{laboratory_typing_method}=~/(disk|disc).diffusion/i;
					$amr_record->{resistant_phenotype} = $sir? $sir: "";
					
					$amr_record->{laboratory_typing_method} = $amr_info->{laboratory_typing_method} unless $amr_record->{laboratory_typing_method};	
					$amr_record->{laboratory_typing_method_version} = $amr_info->{laboratory_typing_method_version} unless $amr_record->{laboratory_typing_method_version};	
					$amr_record->{laboratory_typing_platform} = $amr_info->{laboratory_typing_platform} if $amr_info->{laboratory_typing_platform};	
					$amr_record->{vendor} = $amr_info->{vendor};	
					$amr_record->{testing_standard} = $amr_info->{testing_standard};	
					$amr_record->{testing_standard_year} = $amr_info->{testing_standard_year};	
					$amr_record->{source} = $amr_info->{source}? $amr_info->{source}: $pmid;	

					$amr_record->{owner} = $owner;
					$amr_record->{public} = $public;

					$prev_attrib = $attrib;

				}		

			}else{

				print LOG "$genome_id\t$genome_name\t$attrib\t$value\n";

			}
			
		
		}

		# print last AMR record for a row
		#print "\tpush from outer loop\n";
		push @amr_records, $amr_record if $amr_record->{measurement} || $amr_record->{resistant_phenotype};
		
		if (scalar @amr > 0){
			$metadata_record->{antimicrobial_resistance}->{set} = \@amr;
			$metadata_record->{antimicrobial_resistance_evidence}->{set} = "AMR Panel";
		}
			
		push @metadata_records, $metadata_record if $metadata_record->{genome_id};
			
	}

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

