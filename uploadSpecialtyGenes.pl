#!/usr/bin/env perl 

use strict;
use warnings;
no warnings 'uninitialized';
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
    ["genome_list=s", "File containing list of annotation files"],
    ["spgene_ref=s", "File containing spgene_ref data"],
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->genome_list && $opt->spgene_ref;


#my $spgeneRef = $solrh->getSpGeneRef();

#source,source_id,property,locus_tag,organism,function,classification,pmid,assertion

my $spgeneRef;
 
open SPGREF, $opt->spgene_ref or die "Can't open spgene_ref file: $opt->spgene_ref!!\n" if $opt->spgene_ref;
while (my $line = <SPGREF>){
	chomp $line;
	my ($property,$source,$source_id,$gene_name,$locus_tag,$gene_id,
			$gi,$genus,$species,$organism,$product,$function,$classification,
			$antibiotics_class,$antibiotics,$pmid,$assertion) = split /\t/, $line;

	my $key = $source."_".$source_id;
	$spgeneRef->{$key} = "$property\t$source\t$source_id\t$locus_tag\t$organism\t".
		"$function\t$classification\t$antibiotics_class\t$antibiotics\t$pmid\t$assertion";
	
	#print "$key\t$spgeneRef->{$key}\n";
	
}
close SPGREF;


open LIST, $opt->genome_list or die "Can't open genome_list file: $opt->genome_list!!\n";

while (my $file_name = <LIST>) {
	chomp $file_name;

	my $genome_id = $file_name;
	$genome_id =~s/.*\///;
	$genome_id=~s/.PATRIC.faa.hit//;
	print "Processing $genome_id\n";

	# process each genome files
	if (-f $file_name){ 
		open GENOME, "$file_name" or next;
	}else {
		print "\t$file_name doesn't exist!!\n";
	}  

	# global arrays to hold spgene records
	my @spgenes;

	# new genome, get primary identifiers and existing features
	my $features = getFeatures($genome_id);

	while (my $line = <GENOME>){
		chomp $line;
		my ($patric_id, $spgene_match) = $line=~/([^\t]+)\t(.*)/;
		next unless $features->{$patric_id};

		my $spgene = prepareSPGene($features->{$patric_id},$spgene_match);
  	push @spgenes, $spgene if $spgene;
	}
	close GENOME;

	prepareJsonFile($genome_id, \@spgenes);
	postJsonFile($genome_id) if $opt->commit eq "true";

}

close LIST;


sub getFeatures {

	my ($genome_id) = @_;

	my $core = "/genome_feature";
	my $query = "/select?q=annotation:PATRIC AND feature_type:CDS AND genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,taxon_id,sequence_id,accession,annotation,feature_id,patric_id,alt_locus_tag,refseq_locus_tag,gene,product,owner,public";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrFormat="&wt=json";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my $result = `wget -q -O - "$solrQuery"`;
	
	my $resultObj = decode_json($result);
	my $features;
	
	foreach my $record(@{$resultObj->{'response'}->{'docs'}}){
		$features->{$record->{patric_id}} = $record;

	}

	return $features;

}

sub prepareSPGene {

	my ($feature, $spgene_match) = @_;

	my $spgene;
	my ($msource, $msource_id, $qcov, $scov, $identity, $evalue) = split /\t/, $spgene_match;

	$msource_id=~s/^\S*\|//;

	my $key = $msource.'_'.$msource_id;

	my ($property,$source,$source_id,$locus_tag,$organism,$function, $classification,$antibiotics_class,$antibiotics,$pmid, $assertion) 
		= split /\t/, $spgeneRef->{$key};

	return unless $source_id;

	my ($qgenus, $qspecies, $sgenus, $sspecies);
	($qgenus) = $feature->{genome_name}=~/^(\S+)/;
	($qspecies) = $feature->{genome_name}=~/^(\S+ +\S+)/;
	($sgenus) = $organism=~/^(\S+)/ if $organism;
	($sspecies) = $organism=~/^(\S+ +\S+)/ if $organism;

	my ($same_genus, $same_species, $same_genome, $evidence); 
	$same_genus = 1 if ($qgenus eq $sgenus && $sgenus ne "");
	$same_species = 1 if ($qspecies eq $sspecies && $sspecies ne ""); 
	$same_genome = 1 if ($feature->{genome} eq $organism && $organism ne "") ;

	$evidence = ($feature->{refseq_locus_tag} && $locus_tag && $feature->{refseq_locus_tag} eq $locus_tag)? 'Literature':'BLAT';

	$spgene->{owner} = $feature->{owner};
	$spgene->{public} = $feature->{public};
		
	$spgene->{genome_id} = $feature->{genome_id};	
	$spgene->{genome_name} = $feature->{genome_name};	
	$spgene->{taxon_id} = $feature->{taxon_id};	
		
	$spgene->{feature_id} = $feature->{feature_id};
	$spgene->{patric_id} = $feature->{patric_id};	
	$spgene->{alt_locus_tag} = $feature->{alt_locus_tag};
	$spgene->{refseq_locus_tag} = $feature->{refseq_locus_tag};
		
	$spgene->{gene} = $feature->{gene};
	$spgene->{product} = $feature->{product};

	$spgene->{property} = $property;
	$spgene->{source} = $source;
	$spgene->{property_source} = $property.': '.$source;
	
	$spgene->{source_id} = $source_id;
	$spgene->{organism} = $organism;
	$spgene->{function} = $function;
	$spgene->{classification} = [split /;|,/, $classification]; 
	
	$spgene->{antibiotics_class} = $antibiotics_class; 
	$spgene->{antibiotics} = [split /;/, $antibiotics]; 

	$spgene->{pmid} = [ split /,|;/, $pmid]; 
	$spgene->{assertion} = ""; # $assertion;

	$spgene->{query_coverage} = $qcov; 
	$spgene->{subject_coverage} =  $scov;
	$spgene->{identity} = $identity;
	$spgene->{e_value} = $evalue;

	$spgene->{same_genus} = $same_genus;
	$spgene->{same_species} = $same_species;
	$spgene->{same_genome} = $same_genome;
	$spgene->{evidence} = $evidence;	

	return $spgene;

}


sub prepareJsonFile {

	my($genome_id, $spgenes) = @_;

	print "\tPrepare $genome_id.sp_gene.json\n";
	my $spgene_json = $json->pretty->encode($spgenes);
	open SPG, ">$genome_id.sp_gene.json";
	print SPG "$spgene_json";
	close SPG;

}


sub postJsonFile {

	my($genome_id) = @_;

	print "\tPost $genome_id.sp_gene.json\n";
	`post.update.sh sp_gene $genome_id.sp_gene.json`;
	#`rm $genome_id.genome_feature.json`;
	
}

