#!/usr/bin/perl

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;&rows=100000";

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open IN, $infile;
open JSON, ">$outfile";

my $count = 0;

foreach my $genome (<IN>){

	chomp $genome;
	next if $genome=~/genome_id/;
	$count ++;
	
	my ($genome_length, $gc_content, $gc_count);
	my ($genome_id, $genome_name) = $genome=~/(.*)\t(.*)/;

	print "$genome_id\t$genome_name\n";

	my $core = "/genome_sequence";
	my $query = "/select?q=genome_id:$genome_id";
	my $fields = "&fl=sequence";
	my $solrQuery = $solrServer.$core.$query.$fields.$solrFormat;

	print "$solrQuery\n";

	my @sequences = `wget -q -O - "$solrQuery"`;

	foreach my $sequence (@sequences){
	
		chomp $sequence;
		next if $sequence=~/sequence|^\s*$/;

		$genome_length += length($sequence);
    $gc_count += $sequence=~tr/GCgc//;
	
	}
	
	$gc_content = sprintf("%.2f", ($gc_count*100/$genome_length)) if $genome_length > 0; 

	print JSON "[" if $count == 1;
	print JSON ",\n" unless $count == 1;
  print JSON "{\"genome_id\":\"$genome_id\",\"gc_content\":{\"set\":\"$gc_content\"}}";

}

print JSON "]";
close JSON;
close IN;

#`post.update.dayhoff.sh $core_name $jsonFile`;

