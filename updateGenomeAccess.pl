#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long::Descriptive;
use JSON;

$0 =~ m/([^\/]+)$/;
my $self = $1;

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";
my $json = JSON->new->allow_nonref;

my @genome_ids = ();
my %uid_names = (		'genome' => 'genome_id', 
										'genome_sequence' => 'sequence_id',
										'genome_feature' => 'feature_id', 
										'pathway' => 'id',
										'sp_gene' => 'id',
										'genome_amr' => 'id'
									);


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["genome_list=s", "File containing list of genome ids"],
    ["genome_id=s", "Genome id"],
		["public=s", "public, true"],
		["owner=s", "PATRIC user id for new owner"],
		["add_user=s", "list of user ids, separated by comma"],
		["remove_user=s", "list of user ids, separated by comma"],
		["commit=s", "commit, true|false", {default => "false"}],
		["clean=s", "delete index files, true | false", {default => "false"}],
		[],
		["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless ($opt->genome_list || $opt->genome_id);
die($usage->text) unless ($opt->public || $opt->owner || $opt->add_user || $opt->remove_user);

if ($opt->genome_list){
	open LIST, $opt->genome_list or die "Can't open genome list: $opt->genome_list!";
	@genome_ids = <LIST>;
	close LIST; 
}elsif($opt->genome_id){
	push @genome_ids, $opt->genome_id;
}


foreach my $genome_id (@genome_ids){

	chomp $genome_id;
	next if $genome_id=~/genome_id/;

	my $core = "/genome_feature";
	my $query = "/select?q=genome_id:$genome_id";
	my $fields = "&fl=genome_name";
	my $rows = "&rows=1";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my $genome_name = `wget -q -O - "$solrQuery" | grep -v genome_name`;

	chomp $genome_name;

	print "$genome_id\t$genome_name\n";

	for my $core_name (keys %uid_names){

		my $uid_name = $uid_names{$core_name};

		my @records;
		my @updates;

		my $core = "/$core_name";
		my $query = "/select?q=genome_id:$genome_id";
		my $fields = "&fl=$uid_name";
		my $rows = "&rows=100000";
		my $sort = "";
		my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

		push @records, `wget -q -O - "$solrQuery"`;

		next unless scalar @records > 0;

		foreach my $uid (@records) {
			
			chomp $uid;
			next if $uid=~/$uid_name/;

			my $record;
			$record->{$uid_name} = $uid;
			$record->{public}->{set} = $opt->public=~/true|yes/i ? 1 : 0 if $opt->public;
			$record->{owner}->{set} = $opt->owner if $opt->owner;

			$record->{user_read}->{add} = [split /,/, $opt->add_user] if $opt->add_user; 
			$record->{user_write}->{add} = [split /,/, $opt->add_user] if $opt->add_user;
 
			$record->{user_read}->{remove} = [split /,/, $opt->remove_user] if $opt->remove_user; 
			$record->{user_write}->{remove} = [split /,/, $opt->remove_user] if $opt->remove_user; 
		
			push @updates, $record;

		}

		my $json = $json->pretty->encode(\@updates);
		my $jsonFile = "$genome_id.$core_name.json";
		open FH, ">$jsonFile";			
		print FH "$json";
		close FH;			

		`post.update.sh $core_name $jsonFile` if $opt->commit=~/true|yes/i;

		`rm $jsonFile` if $opt->clean=~/true|yes/;
	
	}

}
