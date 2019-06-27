#!/usr/bin/perl 


my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my $famFile = $ARGV[0];
open FAM, "$famFile";

open TAB, ">patricfams";

my $taxon = getTaxon();
#my $genomes = getGenomes();

my $prev_genome_id = "";

while (my $row = <FAM>){

		chomp $row;

		my ($pgfam_id, $patric_id, $function, $genus, $plfam_idx, $plfam_id);

		my @cols = split(/\t/, $row);
		$pgfam_id = $cols[0];
		$patric_id = $cols[3];
		$fam_function = $cols[5];
		$plfam_idx = $cols[6];
		$genus = $cols[7];

		my ($genome_id) = $patric_id=~/fig\|(.*).peg/;

		if ($genome_id ne $prev_genome_id){
			print "$genome_id\n";

			print JSON "]" if $prev_genome_id;
			close JSON if $prev_genome_id;
			open JSON, ">$genome_id.json";
			print JSON "[";

			$FID = getFeatureIDs($genome_id);

			$count=0;
		
		}
		
		$feature_id = $FID->{$patric_id};
		$pgfam_id=~s/GF/PGF_/;
		$plfam_id="PLF_".$taxon->{$genus}."_".sprintf("%08d", $plfam_idx);

		print TAB "$genus\t$patric_id\t$pgfam_id\t$plfam_id\t$feature_id\n";
		
		print JSON ",\n" unless $count == 0;		
		print JSON "{\"feature_id\":\"$feature_id\",\"pgfam_id\":{\"set\":\"$pgfam_id\"},\"plfam_id\":{\"set\":\"$plfam_id\"}}";

		$prev_genome_id = $genome_id;
		$count++;

}

print JSON "]";
close JSON;
close TAB;

sub getTaxon {

	my $core = "taxonomy";
	my $query = "/select?q=taxon_rank%3Agenus";
	my $fields = "&fl=taxon_name%2Ctaxon_id";
	my $rows = "&rows=100000";
	my $sort = "&sort=taxon_name asc";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my @genera = `wget -q -O - "$solrQuery" | grep -v taxon_name`;

	my %taxon=();
	foreach my $genus (@genera){
		my ($taxon_name, $taxon_id) = $genus=~/(.*)\t(.*)/;
		$taxon{$taxon_name}=$taxon_id;
	}

	return \%taxon;

}


sub getGenomes {

	my $core = "genome";
	my $query = "/select?q=owner:PATRIC AND public:1";
	my $fields = "&fl=genome_id";
	my $rows = "&rows=1000000";
	my $sort = "&sort=genome_name asc";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my @genomes = `wget -q -O - "$solrQuery" | grep -v genome_id`;

	return \@genomes;

}


sub getFeatureIDs {

	my ($genome_id) = @_;

	my $core = "genome_feature";
	my $query = "/select?q=annotation:PATRIC AND feature_type:CDS AND genome_id:$genome_id";
	my $fields = "&fl=patric_id,feature_id";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my @features = `wget -q -O - "$solrQuery" | grep -v patric_id`;

	my %FID = ();
	foreach my $feature (@features){
		my ($patric_id, $feature_id) = $feature=~/(.*)\t(.*)/;
		$FID{$patric_id} = $feature_id;
	}

	return \%FID;

}
