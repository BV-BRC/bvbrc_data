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
my $usage = $self. qq( --genome_list genome_list.tab --annotation PATRIC|RefSeq --public true|false --help);

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my ($help, $genome_list, $genome_id, $public);
my @genome_ids = ();
my @annotations = ();
my %uid_names = (		'genome' => 'genome_id', 
										'genome_sequence' => 'sequence_id',
										'genome_feature' => 'feature_id', 
										'pathway' => 'id',
										'subsystem' => 'id',
										'sp_gene' => 'id',
										'genome_amr' => 'id'
									);


my $opts = GetOptions(
    "help" => \$help,
    "genome_list=s" => \$genome_list,
    "genome_id=s" => \$genome_id,
    "annotation=s" => \@annotations,
		"public=s" => \$public
);

if (!$opts || $help || !($genome_list || $genome_id)){
  warn qq(\n   usage: $usage\n\n);
  exit(0);
}

@annotations = split(/,/, join(',', @annotations)) if scalar @annotations;

if ($genome_list){
	open LIST, $genome_list;
	@genome_ids = <LIST>;
	close LIST; 
}

if ($genome_id){
	push @genome_ids, $genome_id;
} 

$public = 1 if $public=~/true/i;
$public = 0 if $public=~/false/i;

foreach my $genome_id (@genome_ids){

	chomp $genome_id;
	next if $genome_id=~/genome_id/;

	foreach my $annotation (@annotations){

		my $core = "/genome_feature";
		my $query = "/select?q=genome_id:$genome_id AND annotation:$annotation";
		my $fields = "&fl=genome_name";
		my $rows = "&rows=1";
		my $sort = "";
		my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

		my $genome_name = `wget -q -O - "$solrQuery" | grep -v genome_name`;

		chomp $genome_name;

		#if ($genome_name_old && $genome_name ne $genome_name_old){

		print "$genome_id\t$genome_name\t$annotation\tPublic:$public\n";

		for my $core_name (keys %uid_names){

			my $uid_name = $uid_names{$core_name};

			next if ($annotation eq "RefSeq" && $core_name=~/sp_gene|genome_amr/);

			my @records;

			my $core = "/$core_name";
			my $query = "/select?q=genome_id:$genome_id";
			$query .= " AND annotation:$annotation" if ($core_name=~/genome_feature|pathway/);
			my $fields = "&fl=$uid_name";
			my $rows = "&rows=100000";
			my $sort = "";
			my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

			#print "$solrQuery\n";

			push @records, `wget -q -O - "$solrQuery"`;

			next unless scalar @records > 0;

			my $jsonFile = "$genome_id.$annotation.$core_name.json";
			open JSON, ">$jsonFile";

			my $count = 0;
			print JSON "[";

			foreach my $uid (@records) {
			
				chomp $uid;
				next if $uid=~/$uid_name/;

				$count++;
				print JSON ",\n" unless $count == 1;
				print JSON "{\"$uid_name\":\"$uid\",\"public\":{\"set\":\"$public\"}}";
		
			}

			print JSON "]";

			close JSON;

			`post.update.sh $core_name $jsonFile`;

			#`rm $jsonFile`;
	
		 }

		#}

	}

}
