#!/usr/bin/perl

###########################################################
=pod


=cut
###########################################################

use Connection;
use Utils;
use strict;

my $dbh = Connection->getDBHandleFromConfigFile();
my $params = {};

main();

sub main
{ 
	initialize();
 
 	Utils->usage($params, ['algorithm']);
  
	my @algorithms = split(/\s/, $params->{algorithm});
	
	my %format = ();
	my @formats = split(/\s/, $params->{format});
	foreach my $form(@formats){ $format{$form} = 1};

	my $genome = "\'". $params->{genome} ."\'";
	
	if ($params->{genome_list}){
		open(LIST, $params->{genome_list}) or die "Can't open genome_list file!!\n";
		$genome = "";
		while (my $entry = <LIST>){
			chomp $entry;
			$entry =~s/'/''/g;
			$genome .= "\'". $entry ."\', " ;
		}
		$genome=~s/, $//;
	}

  my @genomes = getGenomes($genome);

	foreach my $genome (@genomes){

		my($genome_name, $common_name)=$genome=~/(.*)\t(.*)/;
		#$genome_name=~s/'/''/g;

		print "Preparing files for $genome_name...\n";
		
		`mkdir $common_name`;

		foreach my $algorithm(@algorithms){
				
			prepareFSA($genome_name, $common_name, $algorithm);
			prepareTBL($genome_name, $common_name, $algorithm);
			
			prepareGB($genome_name, $common_name, $algorithm);
			#prepareGFF($genome_name, $common_name, $algorithm);

			`rm $common_name\/*.ecn $common_name\/*.f* $common_name\/*.sqn $common_name\/*.tbl`;

		}

	}	

}

sub initialize{

	$params->{algo_file} = { "PATRIC" => "PATRIC", "RefSeq" => "RefSeq", "BRC" => "BRC"};
	
	$params->{algo_name} = { "PATRIC" => "RAST", "RefSeq" => "RefSeq", "BRC" => "Curation"};

}

sub getGenomes {

  my ($genome) = @_;

	my @genomes;

	#$genome=~s/^'|'$//g; ####

	my $sql = "
	select distinct gi.genome_name, gi.common_name
	from cas.genomeinfo gi, cas.sequenceinfo si 
	where  gi.genome_info_id = si.genome_info_id ";
	
	$sql.="and genome_name in ($genome) " unless($genome eq "\'all\'");
	#$sql.="and genome_name like '%$genome%' ";
	$sql.="order by genome_name";

	#print "$sql\n";

	my $st = $dbh->prepare($sql);
   	$st->execute;

	while (my ($genome_name, $common_name) = $st->fetchrow_array){
		push(@genomes, "$genome_name\t$common_name");
	}

	return @genomes;
}


sub prepareFSA {

	my ($genome_name, $common_name, $algorithm) = @_;

	my $outFile = "$common_name/$common_name.".$params->{algo_file}->{$algorithm}.".fsa";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";
 
	print "\t$outFile\n";

  my $count = 0;
  my $sql = "
		select si.sequence_info_id, si.accession, ns.description, gi.genome_name, ns.sequence 
		from cas.genomeinfo gi, cas.sequenceinfo si, dots.nasequence ns 
		where gi.genome_info_id = si.genome_info_id
		and si.na_sequence_id = ns.na_sequence_id
		and ns.sequence is not null
		and genome_name = ?
		order by accession
  ";

  my $sth = $dbh->prepare($sql);
	$sth->execute($genome_name);
 
 	my $rows = 0;	
	while (my($sequence_info_id, $accession, $description, $genome_name, $sequence) = $sth->fetchrow_array){
		
		$sequence=createFastaSequence($sequence);

		print OUT ">$accession [organism=$genome_name] $description\n$sequence\n" if ($sequence ne "");
   
		$rows++;
	
	}
	
	close OUT;

	`rm $outFile` if ($rows == 0);

}


sub	prepareTBL {

	my ($genome_name, $common_name, $algorithm) = @_;

	my $outFile = "$common_name/$common_name.".$params->{algo_file}->{$algorithm}.".tbl";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";

	my @features = prepareFeatures($genome_name, $algorithm);
	my $go = prepareGO($genome_name, $algorithm);
	my $ec = prepareEC($genome_name, $algorithm);


	my $prevAccession = "";

	foreach my $feature (@features){
	
		my ($genome_name, $accession, $annotation, $feature_type, $na_feature_id, $locus_tag, $start, $end, $strand, $na_length, $gene, $product, $label, $is_pseudo, $bound_moiety, $anticodon, $protein_id, $pseed_id)
				= @$feature;
	
		if ($accession ne $prevAccession){
			print OUT ">Feature $accession\n";
			$prevAccession = $accession;
		}		

		print OUT "$start\t$end\t$feature_type\n" if ($strand eq '+');
		print OUT "$end\t$start\t$feature_type\n" if ($strand eq '-');

		print OUT "\t\t\tlocus_tag\t$locus_tag\n" if $locus_tag;		
		print OUT "\t\t\tproduct\t$product\n" if $product;		
		print OUT "\t\t\tgene\t$gene\n" if $gene;		
		print OUT "\t\t\tlabel\t$label\n" if $label;		
		print OUT "\t\t\tpseudo\t \n" if $is_pseudo;		
		print OUT "\t\t\tbound_moiety\t$bound_moiety\n" if $bound_moiety;		
		print OUT "\t\t\tanticodon\t$anticodon\n" if $anticodon;
		print OUT "\t\t\tprotein_id\t$protein_id\n" if $protein_id;		
		print OUT "\t\t\tdb_xref\tSEED:$pseed_id\n" if $pseed_id;		
		print OUT "\t\t\tdb_xref\tPATRIC:$na_feature_id\n" if $pseed_id;		


		foreach my $go_term (@{$go->{$na_feature_id}}){
			print OUT "\t\t\tgo_function\t$go_term\n";
		}
		
		foreach my $ec_no (@{$ec->{$na_feature_id}}){
			print OUT "\t\t\tEC_number\t$ec_no\n";
		}

	}

	close OUT;

}

sub prepareFeatures(){
	
	my ($genome_name, $algorithm) = @_;
	my @features;

	my $sql = "
		select genome_name, accession, 
			decode(df.algorithm,'Curation','BRC','RAST','PATRIC','RefSeq') as annotation, name as feature_type, 
			na_feature_id, source_id as locus_tag, 
			start_max, end_min, decode(is_reversed, 1, '-', '+') as strand,
			na_length, gene, product,
			label, is_pseudo, bound_moiety, anticodon, 
			protein_id, pseed_id
		from app.dnafeature df
		where df.genome_name = ?
		and df.algorithm = ?
		order by df.accession, start_max
	";	
  
	my $sth = $dbh->prepare($sql);
  $sth->execute($genome_name, $params->{algo_name}->{$algorithm});

	while (my @row=$sth->fetchrow_array){
		push @features, [@row];	
	}

	return @features;

}


sub	prepareGO(){

	my ($genome_name, $algorithm) = @_;
	my %go = ();

	my $sql = "
		select distinct df.na_feature_id, go_id, go_term
		from app.dnafeature df, app.gosummary gs
		where df.na_feature_id = gs.na_feature_id
		and df.name = 'CDS'
		and df.genome_name = ?
		and df.algorithm = ?
		order by df.na_feature_id
	";	
 
	my $sth = $dbh->prepare($sql);
  $sth->execute($genome_name, $params->{algo_name}->{$algorithm});

	my $prevFID;
	my @go;
	while (my @row = $sth->fetchrow_array){
		
		my ($na_feature_id, $go_id, $go_term) = @row;			

		if ($na_feature_id != $prevFID){
			$go{$prevFID} = [@go];
			@go=();
		}
		$prevFID = $na_feature_id;
		$go_id=~s/GO://;
		push @go, "$go_term|$go_id";

	}

	return \%go;

}


sub	prepareEC(){

	my ($genome_name, $algorithm) = @_;
	my %ec = ();

	my $sql = "
		select distinct df.na_feature_id, es.ec_number
		from app.dnafeature df, app.ecsummary es
		where df.na_feature_id = es.na_feature_id
		and df.name = 'CDS'
		and df.genome_name = ?
		and df.algorithm = ?
		order by df.na_feature_id
	";	
	
	my $sth = $dbh->prepare($sql);
  $sth->execute($genome_name, $params->{algo_name}->{$algorithm});

	my $prevFID;
	my @ec;
	while (my @row = $sth->fetchrow_array){
		
		my ($na_feature_id, $ec_no) = @row;

		if ($na_feature_id != $prevFID){
			$ec{$prevFID} = [@ec];
			@ec=();
		}
		$prevFID = $na_feature_id;
		push @ec, "$ec_no";

	}

	return \%ec;	

}


sub prepareGB {

	my ($genome_name, $common_name, $algorithm) = @_;
			
	`tbl2asn -p $common_name -V b -a s`;

	my $gbFile = "$common_name/$common_name.".$params->{algo_file}->{$algorithm}.".gbf";

	open GB, "$gbFile";
	open TMP, ">tmp.gbf";

	my $accession;
	while ( my $line = <GB>){

		($accession) = $line=~/^LOCUS\s+(\S+)\s+/ if $line=~/^LOCUS/;			

		print TMP $line unless($line=~/LOCUS|ACCESSION|VERSION|ORGANISM|SOURCE|\/organism/);
	
		if ($line=~/^LOCUS/){
			($accession) = $line=~/^LOCUS\s+(\S+)\s+/;
			$line =~s/ bp    DNA     circular     / bp    DNA     circular BCT /;
			$line =~s/ bp    DNA     linear       / bp    DNA     linear   BCT /;
			print TMP $line;
		}

		print TMP "ACCESSION   $accession\n" if ($line=~/^ACCESSION/);
		print TMP "VERSION     $accession\n" if ($line=~/^VERSION/);

		print TMP "SOURCE      $genome_name\n" if ($line=~/^SOURCE/);
		print TMP "  ORGANISM  $genome_name\n" if ($line=~/^  ORGANISM/);
		print TMP "                     /organism=\"$genome_name\"\n" if ($line=~/^                     .organism=.unknown./);


=pod
		if ($line=~/^SOURCE/){
			formline <<END, $genome_name;
SOURCE      ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END
		}
		print $^A;  $^A = "";

		if ($line=~/ORGANISM/){
		formline <<END, $genome_name;		  
  ORGANISM  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END
		}
		print $^A;  $^A = "";

		if ($line=~/^ *\/organism/){
		my $tmp = "/organism=\"$genome_name\"";
    formline <<END, $tmp;
		                     ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
END
		}

		print $^A;  $^A = "";
=cut

	}

	close GB,
	close TMP;

	`mv tmp.gbf $gbFile`;
	
}


sub prepareGFF {

	my ($genome_name, $common_name, $algorithm) = @_;

	my $gbFile = "$common_name.".$params->{algo_file}->{$algorithm}.".gbf";
	my $gffFile = "$common_name.".$params->{algo_file}->{$algorithm}.".gff";

	`mv $common_name/$gbFile ./`; 	

  `bp_genbank2gff3.pl $gbFile`;

	`mv $gbFile.gff $gffFile`;
	`perl -pi -e 's/\tGenBank\t/\tPATRIC\t/' $gffFile`;
	`perl -pi -e 's/.*\tPATRIC\texon\t.*\n//' $gffFile`;
	`perl -pi -e 's/.*\tPATRIC\tmRNA\t.*\n//' $gffFile`;
	`perl -pi -e 's/.t01;/;/' $gffFile`;

	`mv $common_name\.* $common_name\/.`;

}


sub createFastaSequence {
        my ($seq) = @_;
        my $fasta = "";
        my $seqLength = length($seq);
        while ($seqLength > 60){
                my $seqTemp = substr ($seq, 0, 60);
                $fasta .= $seqTemp."\n";
                $seq = substr ($seq, 60);
                $seqLength = length($seq);
        }
        $fasta .= $seq."\n";
        return $fasta;
}

