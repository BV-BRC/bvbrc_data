#!/usr/bin/env perl

###########################################################
=pod


=cut
###########################################################

use strict;
#use warnings;
use Getopt::Long::Descriptive;

my $solrServer=$ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;&rows=1000000";


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["genome_list=s", "File containing list of genome ids"],
    ["genome_id=s", "Genome id"],
		["format=s", "File formats, separated by comma", {default => "all"}],
		["annotation=s", "Annotation, PATRIC|RefSeq", {default => "PATRIC"}],
		["metadata!", "Prepare metadata files"],
		["blast!", "Build BLAST DBs"],
		["mash!", "Build Mash sketches"],
		["path=s", "Download directory path", {default => "/vol/patric3/downloads"}],
		[],
		["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless ($opt->genome_list || $opt->genome_id || $opt->metadata);


my $genomeDir = "$opt->{path}/genomes";
my $blastDir = "$opt->{path}/blastdb";
my $mashDir = "$opt->{path}/mash";

`mkdir $genomeDir` unless (-d $genomeDir);
`mkdir $blastDir` if ($opt->{blast} && ! -d $blastDir);
`mkdir $mashDir` if ($opt->{mash} && ! -d $mashDir);


main();

sub main
{ 
	#GetOptions(\%params, "genome_summary", "genome_lineage", "genome_metadata", "annotation=s", "format=s", "genome=s", "genomes=s");
  
	my @annotations = split(/\s|,/, $opt->{annotation});
	
	my %format = ();
	my @genomes = ();
	my @formats = split(/\s|,/, $opt->{format});
	foreach my $form(@formats){ $format{$form} = 1};

	@genomes = getGenomes();
  
	if ($opt->{metadata}){
		prepareGenomeSummary();
		prepareGenomeLineage();
		prepareGenomeMetadata();
		prepareAMRSummary();	
	}

	exit if $opt->{format}=~/none/ || scalar @genomes == 0;
	
	#exit unless $params{format}=~/fna|ffn|frn|faa|ftab|rtab|ptab|path|spgene|gff3|gbf/i;

	foreach my $genome (@genomes){

		my($genome_id, $genome_name)=$genome=~/(.*)\t(.*)/;
		next if $genome_name=~/genome_name/i;
		
		`mkdir $genomeDir/$genome_id`;

		print "Preparing files for $genome_id - $genome_name...\n";

		prepareSequenceFasta($genome_id, $genome_name);

		#next; # only fna files

		foreach my $annotation (@annotations){

      my $features = getFeatures($genome_id, $genome_name, $annotation);
			prepareFeatureTab($genome_id, $genome_name, $annotation, $features);
			prepareGFF3($genome_id, $genome_name, $annotation, $features);	
			prepareFeatureFasta($genome_id, $genome_name, $annotation, $features);

			next unless $annotation=~/PATRIC/;

			preparePathwayTab($genome_id, $genome_name, $annotation);
			prepareSubsystemTab($genome_id, $genome_name, $annotation);
			prepareSPGeneTab($genome_id, $genome_name, $annotation);
			
			buildBlastDB($genome_id, $genome_name, $annotation) if $opt->{blast};
			buildMashSketch($genome_id, $genome_name, $annotation) if $opt->{mash}
			
		}
		
	}	

}


sub getGenomes {

	my @genomes;

	if ($opt->{genome_list} && $opt->{genome_list} ne "all"){
		open(LIST, $opt->{genome_list}) or die "Can't open genomes file!!\n";
		
		while (my $genome_id = <LIST>){
			chomp $genome_id;
			print "$genome_id\n";
			
			my $core = "/genome";
			#my $query = "/select?q=genome_id:$genome_id AND public:1";
			my $query = "/select?q=genome_id:$genome_id";
			my $fields = "&fl=genome_id,genome_name";
			my $sort = "&sort=genome_name+asc";
			my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

			push @genomes, `wget -q -O - "$solrQuery"`;
		}
		close LIST;
	}

	if ($opt->{genome_list} && $opt->{genome_list} eq "all"){
			my $core = "/genome";
			my $query = "/select?q=genome_id:*.* AND public:1 AND taxon_lineage_ids:*";
			my $fields = "&fl=genome_id,genome_name";
			my $sort = "&sort=genome_name+asc";
			my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
	
			push @genomes, `wget -q -O - "$solrQuery"`;
	}
	
	if ($opt->{genome_id}){
			print "$opt->{genome_id}\n";
			my $core = "/genome";
			my $query = "/select?q=genome_id:$opt->{genome_id}";
			my $fields = "&fl=genome_id,genome_name";
			my $solrQuery = $solrServer.$core.$query.$fields.$solrFormat;
		
			push @genomes, `wget -q -O - "$solrQuery"`;
	}

	return @genomes;
}

sub prepareGenomeSummary(){

	my $outFile = "genome_summary";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";
	
	my $core = "/genome";
	my $query = "/select?q=genome_id:*.* AND public:1 AND taxon_lineage_ids:131567";
	my $fields = "&fl=genome_id,genome_name,genome_id,taxon_id,genome_length,genome_status,chromosomes,plasmids,contigs,patric_cds,refseq_cds";
	my $sort = "&sort=genome_name+asc";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
			
	`wget -q -O $outFile "$solrQuery"`;

	my $lines = `wc -l $outFile`;
	`rm $outFile` if $lines=~/^0 |^1 /; 

}


sub prepareGenomeLineage(){

	my $outFile = "genome_lineage";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";
	
	my $core = "/genome";
	my $query = "/select?q=genome_id:*.* AND public:1 AND taxon_lineage_ids:131567";
	my $fields = "&fl=genome_id,genome_name,taxon_id,kingdom,phylum,class,order,family,genus,species,taxon_lineage_ids,taxon_lineage_names";
	my $sort = "&sort=genome_name+asc";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
			
	`wget -q -O $outFile "$solrQuery"`;

	my $lines = `wc -l $outFile`;
	`rm $outFile` if $lines=~/^0 |^1 /; 

}

sub prepareGenomeMetadata(){

	my $outFile = "genome_metadata";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";	

	my $core = "/genome";
	my $query = "/select?q=genome_id:*.* AND public:1 AND taxon_lineage_ids:131567";
	my $fields = "&fl=".join(',', genomeAttributes());
	my $sort = "&sort=genome_name+asc";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
		
	#print "$solrQuery\n";
	
	`wget -q -O $outFile "$solrQuery"`;

	my $lines = `wc -l $outFile`;
	`rm $outFile` if $lines=~/^0 |^1 /; 

}


sub prepareAMRSummary(){

	my $outFile = "PATRIC_genomes_AMR.txt";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";
	
	my $core = "/genome_amr";
	my $query = "/select?q=public:1 AND NOT laboratory_typing_method:Computational";
	my $fields = "&fl=genome_id,genome_name,taxon_id,antibiotic,resistant_phenotype,".
		"measurement,measurement_sign,measurement_value,measurement_unit,".
		"laboratory_typing_method,laboratory_typing_method_version,laboratory_typing_platform,".
		"vendor,testing_standard,testing_standard_year,source";
	my $sort = "&sort=genome_name ASC,genome_id ASC,antibiotic ASC";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
	
	`wget -q -O $outFile "$solrQuery"`;

	my $lines = `wc -l $outFile`;
	#`rm $outFile` if $lines=~/^0 |^1 /; 

}


sub	prepareSequenceFasta(){

	my ($genome_id, $genome_name) = @_;

	my $outFile = "$genomeDir/$genome_id/$genome_id.fna";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";
	
	my $core = "/genome_sequence";
	my $query = "/select?q=genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,accession,gi,description,sequence";
	my $sort = "&sort=accession+asc";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

	my @result = `wget -q -O - "$solrQuery"`;

	return unless scalar @result;

	foreach my $record (@result) {
		next if $record=~/^genome_id/i;
		my ($genome_id, $genome_name, $accession, $gi, $description, $sequence) = split(/\t/, $record);
		next unless $sequence;
		$sequence=createFastaSequence($sequence);
		print OUT ">$accession   $description   [$genome_name | $genome_id]\n$sequence\n";
	}
	close OUT;

	my $lines = `wc -l $outFile`;
	`rm $outFile` if $lines=~/^0 |^1 /; 

}


sub	getFeatures(){

	my ($genome_id, $genome_name, $annotation) = @_;

	my $core = "/genome_feature";
	my $query = "/select?q=genome_id:$genome_id+AND+annotation:$annotation";
	my $fields = "&fl=genome_id,genome_name,accession,".
								"annotation,feature_type,".
								"patric_id,refseq_locus_tag,".
								"start,end,strand,na_length,".
								"gene,product,".
								"plfam_id,pgfam_id,".
								"na_sequence_md5,aa_sequence_md5";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

	my @result = `wget -q -O - "$solrQuery"`;
	return unless scalar @result;

	# get NA and AA sequences
	my %md5_hash = ();
	foreach my $record (@result){
		chomp $record;
		next if $record=~/^genome_id/i;
		my @fields = split(/\t/, $record);
		$md5_hash{$fields[15]} = 1 if $fields[15];
		$md5_hash{$fields[16]} = 1 if $fields[16];
	}

	getFeatureSequences(\%md5_hash); 	

	my %features = (); 
	foreach my $record (@result){
		if ($record=~/^genome_id/i){
			$features{"header"} = $record."\tna_sequence\taa_sequence";
			next;
		}	
		my @fields = split(/\t/, $record);
		my $key = $fields[2].".".sprintf("%08d", $fields[7]).".$fields[4]"; #accession.start.feature_type
		$features{$key} = $record."\t$md5_hash{$fields[15]}\t$md5_hash{$fields[16]}";
	}

	return \%features;

}

sub prepareFeatureTab(){

	my ($genome_id,$genome_name,$annotation,$features) = @_;

	my $outFile = "$genomeDir/$genome_id/$genome_id.$annotation.features.tab";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	my ($header) = $features->{"header"}=~/^(.*)\t\S*\t\S*\t\S*\t\S*$/;
	print OUT "$header\n";

	foreach my $feature (sort keys %$features){
		my ($info) = $features->{$feature}=~/^(.*)\t\S*\t\S*\t\S*\t\S*$/;
		next if $info=~/genome_id/i;
		print OUT "$info\n";
	}

	close OUT;

}


sub	prepareFeatureFasta(){

	my ($genome_id,$genome_name,$annotation,$features) = @_;

	print "\t$genomeDir/$genome_id/$genome_id.$annotation.ffn|frn|faa\n";
	
  my $ffnFile = "$genomeDir/$genome_id/$genome_id.$annotation.ffn";
	open FFN, ">$ffnFile" or die "Sorry, could not open output file $ffnFile !\n";
  
	my $frnFile = "$genomeDir/$genome_id/$genome_id.$annotation.frn";
	open FRN, ">$frnFile" or die "Sorry, could not open output file $frnFile !\n";
  
	my $faaFile = "$genomeDir/$genome_id/$genome_id.$annotation.faa";
	open FAA, ">$faaFile" or die "Sorry, could not open output file $faaFile !\n";
  
	foreach my $feature (sort keys %$features) {

		my (  $genome_id, $genome_name, $accession, $annotation, $feature_type,
			$patric_id, $refseq_locus_tag,
			$start, $end, $strand, $na_length, $gene, $product,
			$plfam_id, $pgfam_id, 
			$na_sequence_md5,$aa_sequence_md5,
			$na_sequence, $aa_sequence) = split(/\t/, $features->{$feature});	

		next if $genome_id=~/genome_id/;
	
		$na_sequence=createFastaSequence($na_sequence) if $na_sequence;
		$aa_sequence=createFastaSequence($aa_sequence) if $aa_sequence;
		
		my $fasta_id;	
		if ($annotation eq 'PATRIC'){
			$fasta_id = "$patric_id". ($refseq_locus_tag? "|$refseq_locus_tag":"");
		}elsif($annotation eq 'RefSeq'){
			$fasta_id = $refseq_locus_tag;
		}

		print FFN ">$fasta_id   $product   [$genome_name | $genome_id]\n$na_sequence" if $feature_type=~/CDS/ && $na_sequence;
		print FRN ">$fasta_id   $product   [$genome_name | $genome_id]\n$na_sequence" if $feature_type=~/RNA/ && $na_sequence;
		print FAA ">$fasta_id   $product   [$genome_name | $genome_id]\n$aa_sequence" if $feature_type=~/CDS/ && $aa_sequence;
	
	}

	close FFN;
	close FRN;
	close FAA;

}


sub	prepareGFF3 {

	my ($genome_id, $genome_name, $annotation, $features) = @_;

	my $outFile = "$genomeDir/$genome_id/$genome_id.$annotation.gff";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";

	print OUT "##gff-version 3\n";
	print OUT "#Genome: $genome_id|$genome_name\n";
	print OUT "#Date:".`date "+%m/%d/%Y"`."\n";

	foreach my $feature (sort keys %$features) {
		
		chomp $feature;

		my (	$genome_id, $genome_name, $accession, $annotation, $feature_type,
					$patric_id, $refseq_locus_tag,
					$start, $end, $strand, $na_length, $gene, $product,
					$plfam_id, $pgfam_id,
					$na_sequence_md5,$aa_sequence_md5,
					$na_sequence, $aa_sequence) = split(/\t/, $features->{$feature});
		
		next if $genome_id=~/genome_id/;
	
		$feature_type = 'region' if $feature_type eq 'source';
		$feature_type = 'transcript' if $feature_type eq 'misc_RNA';
		#$feature_type = 'pseudogene' if ($feature_type eq 'gene' && $is_pseudo == 1 );
		
		$product=~s/;/%3B/g;
		$product=~s/,/%2C/g;
		$product=~s/=/%3D/g;
		$product=~s/&/%26/g;

		#$go=~s/\|([^;]*);/,/g;
		#$go=~s/\|.*$//;
		#$ec=~s/\|([^;]*);/,/g;
		#$ec=~s/\|.*$//;

		if ($feature_type eq 'region'){
				print OUT "##sequence-region\taccn|$accession\t$start\t$end\n";
				next;
		}

		print OUT "accn|$accession\t$annotation\t$feature_type\t$start\t$end\t.\t$strand\t0\t";
		
		print OUT "ID=$patric_id" if $annotation eq 'PATRIC';
		print OUT "ID=$refseq_locus_tag" if $annotation eq 'RefSeq';
		print OUT ";locus_tag=$refseq_locus_tag" if $refseq_locus_tag;
		print OUT ";product=$product" if $product;
		print OUT ";gene=$gene" if $gene;
		#print OUT ";Ontology_term=$go" if $go;
		#print OUT ";ec_number=$ec" if $ec;

		print OUT "\n";
	
	}

	close OUT;
	
	my $lines = `wc -l $outFile`;
	`rm $outFile` if $lines=~/^0 |^1 /; 

}


sub	preparePathwayTab(){

	my ($genome_id, $genome_name, $annotation) = @_;
	
	my $core = "/pathway";
	my $query = "/select?q=genome_id:$genome_id+AND+annotation:$annotation";
	my $fields = "&fl=genome_id,genome_name,".
								"patric_id,refseq_locus_tag,gene,product,".
								"ec_number,ec_description,pathway_id,pathway_name";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

  my @result = `wget -q -O - "$solrQuery"`;
  return unless scalar @result;
	
	my $outFile = "$genomeDir/$genome_id/$genome_id.$annotation.pathway.tab";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";

  my %features = ();
  foreach my $record (@result){
    if ($record=~/^genome_id/i){
      print OUT $record; #header
      next;
    }
    my @fields = split(/\t/, $record);
		my $rank = $fields[2];
		$rank=~s/.*peg.//;
    my $key = sprintf("%06d",$rank).".$fields[6].$fields[8]"; #accession.start
		$features{$key} = $record;
  }

	foreach my $key (sort keys %features){
		print OUT "$features{$key}";
	}

	close OUT;	
	
}


sub	prepareSubsystemTab(){

	my ($genome_id, $genome_name, $annotation) = @_;
	
	my $core = "/subsystem";
	my $query = "/select?q=genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,".
								"patric_id,refseq_locus_tag,gene,product".
								"role_name,superclass,class,subclass,subsystem_name";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;

  my @result = `wget -q -O - "$solrQuery"`;
  return unless scalar @result;
	
	my $outFile = "$genomeDir/$genome_id/$genome_id.$annotation.subsystem.tab";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";

  my %features = ();
  foreach my $record (@result){
    if ($record=~/^genome_id/i){
      print OUT $record; #header
      next;
    }
    my @fields = split(/\t/, $record);
		my $rank = $fields[2];
		$rank=~s/.*peg.//;
    my $key = sprintf("%06d",$rank).".$fields[6].$fields[10]"; #accession.start
		$features{$key} = $record;
  }

	foreach my $key (sort keys %features){
		print OUT "$features{$key}";
	}

	close OUT;	
	
}


sub	prepareSPGeneTab(){

	my ($genome_id, $genome_name, $annotation) = @_;
	
	my $core = "/sp_gene";
	my $query = "/select?q=genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,".
								"patric_id,refseq_locus_tag,alt_locus_tag,gene,product,".
								"property,source,evidence,source_id,organism,function,classification,pmid,".
								"query_coverage,subject_coverage,identity,e_value";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
	
  my @result = `wget -q -O - "$solrQuery"`;
  return unless scalar @result;
	
	my $outFile = "$genomeDir/$genome_id/$genome_id.$annotation.spgene.tab";
  open OUT, ">$outFile" or die "Sorry, could not open output file $outFile !\n";

	print "\t$outFile\n";

  my %features = ();
  foreach my $record (@result){
    if ($record=~/^genome_id/i){
      print OUT $record; #header
      next;
    }
    my @fields = split(/\t/, $record);
		my $rank = $fields[2];
		$rank=~s/.*peg.//;
    my $key = sprintf("%06d", $rank).".$fields[7].$fields[8]"; #accession.start
    $features{$key} = $record;
  }

	foreach my $key (sort keys %features){
		print OUT "$features{$key}";
	}

	close OUT;	

}


sub getFeatureSequences {
	my ($md5_hash_ref) = @_;

	my $count = 0;
	my $md5_values = "";
	my $size = keys %{$md5_hash_ref};

	foreach my $md5(keys %{$md5_hash_ref}){	
		$count++;
		next if $md5 eq "";
		if ($md5_values eq ""){
			$md5_values="$md5";
		}else{
			$md5_values.=" OR $md5";
		}

		if ( $count%100==0 || $count == $size){
			my $core = "/feature_sequence";
			my $query = "/select?q=md5:($md5_values)";
			my $fields = "&fl=md5,sequence";
			my $sort = "";
			my $solrQuery = $solrServer.$core.$query.$fields.$sort.$solrFormat;
			my @sequences = `wget -q -O - "$solrQuery"`;

			my $i = 0;	
			foreach my $record (@sequences){
				chomp $record;
				$i++;
				my ($md5, $sequence) = split(/\t/, $record);
				$md5_hash_ref->{$md5} = $sequence;
			}
			$md5_values = "";
		}
	}

}


sub createFastaSequence {
        my ($seq) = @_;
				chomp $seq;
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


sub	buildBlastDB(){

	my ($genome_id,$genome_name,$annotation) = @_;

	print "\tBuilding BLAST DBs for $genome_id - $genome_name\n";

  `makeblastdb -dbtype nucl -in $genomeDir/$genome_id/$genome_id.fna -out $blastDir/$genome_id.fna`;
	`makeblastdb -dbtype prot -in $genomeDir/$genome_id/$genome_id.PATRIC.faa -out $blastDir/$genome_id.PATRIC.faa`;
  `makeblastdb -dbtype nucl -in $genomeDir/$genome_id/$genome_id.PATRIC.ffn -out $blastDir/$genome_id.PATRIC.ffn`;
  `makeblastdb -dbtype nucl -in $genomeDir/$genome_id/$genome_id.PATRIC.frn -out $blastDir/$genome_id.PATRIC.frn`;

}	


sub	buildMashSketch(){

	my ($genome_id,$genome_name,$annotation) = @_;

	print "\tBuilding Mash sketch for $genome_id - $genome_name\n";

	`mash sketch -s 1000 -o $mashDir/$genome_id.msh $genomeDir/$genome_id/$genome_id.fna`;

}


sub genomeAttributes {

	my @genomeAttributes = ("genome_id", "genome_name", "organism_name", "taxon_id", "genome_status",
				"strain", "serovar", "biovar", "pathovar", "mlst", "other_typing",
				"culture_collection", "type_strain",
				"completion_date", "publication",
				"bioproject_accession", "biosample_accession", "assembly_accession", "genbank_accessions",
				"refseq_accessions",
				"sequencing_centers", "sequencing_status", "sequencing_platform", "sequencing_depth", "assembly_method",
				"chromosomes", "plasmids", "contigs", "sequences", "genome_length", "gc_content",
				"patric_cds", "brc1_cds", "refseq_cds",
				"isolation_site", "isolation_source", "isolation_comments", "collection_date",
				"isolation_country", "geographic_location", "latitude", "longitude", "altitude", "depth", "other_environmental",
				"host_name", "host_gender", "host_age", "host_health", "body_sample_site", "body_sample_subsite", "other_clinical",
				"antimicrobial_resistance", "antimicrobial_resistance_evidence",
				"gram_stain", "cell_shape", "motility", "sporulation", "temperature_range", "optimal_temperature", "salinity", 
				"oxygen_requirement","habitat","disease", "comments", "additional_metadata");
		
	return @genomeAttributes;

}
