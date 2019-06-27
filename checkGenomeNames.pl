#!/usr/bin/env perl

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;&rows=1000000";

my $core = "/genome";
my $query = "/select?q=genome_id:*.* AND public:true";
my $fields = "&fl=genome_id,genome_name,strain";
my $sort = "&sort=genome_name ASC";
my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

#print "$solrQuery\n";

my @genomes = `wget -q -O - "$solrQuery"`;

foreach my $genome (@genomes){

	my ($genome_id, $genome_name, $strain) = $genome=~/(.*)\t(.*)\s*\t\s*(.*)\s*\n/;

	#print "$genome_id\t$genome_name\n" if scalar (split / /, $genome_name) <=2;
	#next;

	next unless $strain;
	next if $strain=~/^ *(-|missing|na|n\/a|not|nd|unknown) */i;
	next if $genome_name=~/ (strain|st|str|str\.|sp\.) /;

	#$genome_name=~s/[\W_-]+//g;
	#$strain=~s/^.*[ :]+//;
	#$strain=~s/[\W_-]+//g;
	$strain=~s/^\s*(strain|st|str|str\.) *//i;
	$strain=~s/\?*//g;

	next if $genome_name=~/$strain/i;

	my $genome_name_new = $genome_name." strain $strain";

	#print "$genome_id\t$genome_name\t$strain\t$genome_name_new\n" if scalar (split / /, $genome_name) <=2;
	print "$genome_id\t$genome_name\t$strain\t$genome_name_new\n";

}

