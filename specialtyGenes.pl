#!/usr/bin/env perl

use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Bio::SearchIO; 

$0 =~ m/([^\/]+)$/;
my $self = $1;

my $usage = $self. qq( --in InFile --out OutFile [--program (blast|blat)] [--db (VFDB|Victors|PATRIC_VF|ARDB|CARD|DrugBank|TTD|Human|all)] [--filter (Y|N)] [--description (Y|N)]);

if (not @ARGV){
	warn qq(\n   usage: $usage\n\n);
 	exit(0);	
}	

my $blastdir="$Bin";
#my @BlastDBs=('VFDB', 'Victors', 'PATRIC_VF', 'CARD', 'NDARO', 'DrugBank', 'TTD', 'TCDB');
my @BlastDBs=('VFDB', 'Victors', 'PATRIC_VF', 'ARDB', 'CARD', 'NDARO', 'Human', 'DrugBank', 'TTD', 'TCDB');

my ($help, $in, $out, $program, $db, $filter, $description) = ('', '', '', 'blat', 'all', 'Y', 'N');

my $opts = GetOptions(
		'help' => \$help,
		'in=s' => \$in,
		'out=s' => \$out,
		'program=s' => \$program,
		'db=s' => \$db,
		'filter=s' => \$filter,
		'description=s' => \$description
	);

if (!$opts || $help || !$in){
	warn qq(\n   usage: $usage\n\n);
	exit(0);
}

$out = "$in.hit" unless $out;
$out =~s/^.*\///;

my $blastreport=$out;
$blastreport=~s/\.hit$|$/.br/;

`rm $out $blastreport`;

my @blastdbs=split(/\s|,/, $db);
@blastdbs = @BlastDBs if ($db eq "all");

print "Processing $in\n";

open(HIT,">$out");

my $header;
if ($description eq 'Y'){
	$header = "QID\tQAnnotation\tQOrganism\tDatabase\tSubID\tSubAnnotation\tSubOrganism\tQueryCoverage\tSubCoverage\tIdentity\tP-Value\n";
}else{
	$header = "QID\tDatabase\tSubID\tQueryCoverage\tSubCoverage\tIdentity\tP-Value\n" ;
}
print HIT $header;


foreach my $blastdb (sort @blastdbs){

	$blastdb = "$blastdir/$blastdb";

	my $blastcmd = "blastp -query $in -db $blastdb -num_alignments 1 -num_descriptions 1 -evalue 0.0001 -outfmt 0 -out $blastreport";
	my $blatcmd = "blat $blastdb.faa $in -prot -out=blast $blastreport";

	if ($program=~/blast/i){
		print "\t$blastcmd\n";
		`$blastcmd`;
	}elsif ($program=~/blat/i){
		print "\t$blatcmd\n";
		`$blatcmd`;
	}else{

	}

	my $in = new Bio::SearchIO(-format => 'blast', -file=> "$blastreport");

	while( my $result = $in->next_result ) {
		# $result is a Bio::Search::Result::ResultI compliant object

		while( my $hit = $result->next_hit ) {
			# $hit is a Bio::Search::Hit::HitI compliant object
		
			while( my $hsp = $hit->next_hsp ) {
				# $hsp is a Bio::Search::HSP::HSPI compliant object
=pod
				$result->query_name;
				$result->database_name;
				$hit->name;
				$hsp->query;
				$hsp->hit;
				$hsp->length('query');
				$hsp->length('hit');
				$hsp->length('total');
				$hsp->percent_identity;
				$hsp->start('query');
				$hsp->end('query');
				$hsp->start('hit');
				$hsp->start('hit');
				$hsp->pvalue;
				$hsp->significance;
				$hsp->score;
				$hsp->bits;
=cut

        my $qid = $result->query_name;
				$qid = $1 if ($qid=~/(fig\|[^|]*)|/); 
				my ($qAnnot, $qOrg) = ($result->query_description, "");
				($qAnnot, $qOrg) = $result->query_description=~/(.*)\s*\[(.*)\]/ if $result->query_description=~/(.*)\s*\[(.*)\]/;
												        
				my $database = $result->database_name;
				$database =~s/^.*\/|.faa$//g;
				
				my $sid = $hit->name;
				my ($sAnnot, $sOrg) = ($hit->description, "");
				($sAnnot, $sOrg) = $hit->description=~/(.*)\s*\[(.*)\]/ if $hit->description=~/(.*)\s*\[(.*)\]/;
			
				my $id=int(abs($hsp->percent_identity));	
				my $qcov=int(abs( $hsp->length('query') * 100 / $result->query_length ));
				my $scov=int(abs( $hsp->length('hit') * 100 / $hit->length ));
				my $pvalue = $hsp->significance;

				my $sim = "";
				if ($description eq 'Y'){
					$sim = "$qid\t$qAnnot\t$qOrg\t$database\t$sid\t$sAnnot\t$sOrg\t$qcov\t$scov\t$id\t$pvalue\n";
				}else{
					$sim = "$qid\t$database\t$sid\t$qcov\t$scov\t$id\t$pvalue\n";
				}

				if ($filter eq 'Y'){
					print HIT $sim if (( ($qcov>=80 || $scov>=80) && $id>=80 ) && !($database=~/Human/i) );
					print HIT $sim if ( ($qcov>=25 || $scov>=25) && $id>=40 && ($database=~/Human/i) );
				}else{
					print HIT $sim;
				}
				last; # get only first hit from every query  

			}
			last; # get only first hit from every query  
		}
	}
	`rm $blastreport`;
}

close HIT;
