#!/usr/bin/env perl

my $PATH = "/home/mshukla/patric3/patric3_data_infrastructure/data_processing/src";

`mkdir tmp`;
chdir("tmp");
`wget -q "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"`;
`tar xvfzp taxdump.tar.gz`;
chdir("../");

use Bio::Taxon;
use Bio::Tree::Tree;
use  Bio::DB::Taxonomy;

open OUT, ">taxonomy.json";

my $taxon_id_cellular = 131567;
my $taxon_id_bacteria = 2;
my $taxon_id_archaea = 2157;
my $taxon_id_virus = 10239;

my @taxon_id_eukaryota = qw(10090 10116 6239 7227 7955 9031 9544 9606 9669 9823);

my $dbh = Bio::DB::Taxonomy->new(	-source   => 'flatfile',
                                 	-directory=> '/tmp',
																	-nodesfile=> 'tmp/nodes.dmp',
																	-namesfile=> 'tmp/names.dmp');


my $tree = Bio::Tree::Tree->new();

my @taxa = ();
my $taxon;

$taxon = $dbh->get_taxon(-taxonid => $taxon_id_cellular);
push @taxa, $taxon;

$taxon = $dbh->get_taxon(-taxonid => $taxon_id_bacteria);
push @taxa, $taxon;
push @taxa, $taxon->db_handle->get_all_Descendents($taxon);

$taxon = $dbh->get_taxon(-taxonid => $taxon_id_archaea);
push @taxa, $taxon;
push @taxa, $taxon->db_handle->get_all_Descendents($taxon);

foreach my $taxon_id (@taxon_id_eukaryota){
	$taxon = $dbh->get_taxon(-taxonid => $taxon_id);
	push @taxa, $taxon;
	push @taxa, $tree->get_lineage_nodes($taxon);
}

%tax = ();
$taxon = $dbh->get_taxon(-taxonid => $taxon_id_virus);
push @taxa, $taxon;
$tax{$taxon->id} = 1;

print "Getting descendents\n";
foreach my $taxon1 ($taxon->db_handle->each_Descendent($taxon)){
	next if $tax{$taxon1->id};
	$tax{$taxon1->id} = 1;
	push @taxa, $taxon1;
	print $taxon1->id."\n";
	foreach my $taxon2 ($taxon->db_handle->each_Descendent($taxon1)){
		next if $tax{$taxon2->id};
		$tax{$taxon2->id} = 1;
		
		push @taxa, $taxon2;
		print $taxon2->id."\n";
		push @taxa, $taxon->db_handle->get_all_Descendents($taxon2);
		print "count:", scalar @taxa, "\n";
	}
}

print "count:", scalar @taxa, "\n";

print OUT "[\n";

my $count = 0;

foreach my $taxon (@taxa) {

	$count++;

	my $taxon_id = $taxon->id;
	my $taxon_name = $taxon->scientific_name;
	$taxon_name=~s/"/\"/g;
	my $taxon_rank = $taxon->rank;
	my $genetic_code = $taxon->genetic_code;
	my $parent_id = $taxon->parent_id;
	my $division = $taxon->division;

	my @taxon_common_names = ();

	foreach my $name ($taxon->common_names){
		$name=~s/"|^ *| *$//g;
		push @taxon_common_names, $name; 	
	}

	my $tree = Bio::Tree::Tree->new();
	my @lineage = ($tree->get_lineage_nodes($taxon), $taxon);

	my@taxon_lineage_ids = ();
	my @taxon_lineage_names = ();
	my @taxon_lineage_ranks = ();

	foreach my $node (@lineage) { 
		push @taxon_lineage_ids, $node->id; 
		push @taxon_lineage_names, $node->scientific_name; 
		push @taxon_lineage_ranks, $node->rank; 
	} 

	my $taxon_common_names = "[\"".join("\",\"", @taxon_common_names)."\"]" if scalar @taxon_common_names;
	my $taxon_lineage = join(",", @taxon_lineage_names); 
	my $taxon_lineage_ids = "[\"".join("\",\"", @taxon_lineage_ids)."\"]";
	my $taxon_lineage_names = "[\"".join("\",\"", @taxon_lineage_names)."\"]";
	my $taxon_lineage_ranks = "[\"".join("\",\"", @taxon_lineage_ranks)."\"]";

	#print "$taxon_id\t$taxon_name\t$taxon_rank\t$taxon_lineage_ids\t$taxon_lineage_names\t$taxon_lineage_ranks\n";

	## TO DO: 
	# Make "parent id":null , when not available for root taxon
	# Some taxon nmaes have embedded double quotes, they need to be escaped

	print OUT ",\n" unless $count == 1;
	print OUT "", 	
"	{	\"taxon_id\":\"$taxon_id\",
		\"taxon_name\":\"$taxon_name\",
		\"taxon_rank\":\"$taxon_rank\",
		\"genetic_code\":\"$genetic_code\",
		\"division\":\"$division\",";
		
	print OUT "\n",
"		\"parent_id\":$parent_id," if $parent_id;

	print OUT	"\n",
"		\"other_names\":$taxon_common_names," if $taxon_common_names;

	print OUT	"\n",
"		\"lineage\":\"$taxon_lineage\",
		\"lineage_ids\":$taxon_lineage_ids,
		\"lineage_names\":$taxon_lineage_names,
		\"lineage_ranks\":$taxon_lineage_ranks
	}";

}

print OUT "\n]";

close OUT;

`rm -rf tmp`;

#`rm *.dmp gc.parts taxdump.tar.gz readme.txt`;
