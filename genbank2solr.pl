#!/usr/bin/env perl

use FindBin qw($Bin);
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Date::Parse;
use POSIX;
use JSON;
#use Bio::DB::EUtilities;
use XML::Simple;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);

use lib "$Bin";
use Solr;

my $owner = 'p3_viral@patricbrc.org';
my $public = 0;
my $commit = $ARGV[1]=~/commit/? 1 : 0;;
my $annotation = "PATRIC";

my $solrh = Solr->new();
my $json = JSON->new->allow_nonref;

my $infile = $ARGV[0];
#my ($genome_id, $annotation) = $infile=~/([\d\.]+)\.(PATRIC|RefSeq).gbf/;
my ($genome_id) = $infile=~/(.*).gbff/;


open IN, $infile or die "Can't open input GenBank file: $infile\n";
close IN;

print "Processing $infile\n";

my $genome;
my @sequences = ();
my @features = ();
my @pathwaymap = ();
my @spgenemap = (); 
my $featureIndex;


my $usage = "genbank2solr.pl genbankFile\n";
my $infile = $infile or die "No input file";
my $outfile = $infile;
$outfile=~s/(.gb|.gbf|.gbff)$//;	

getGenomeInfo();
getGenomeSequences();
getGenomeFeatures();
writeJson();
postJson() if $commit;



sub writeJson {

	my $genome_json = $json->pretty->encode($genome);
	my $sequence_json = $json->pretty->encode(\@sequences);
	my $feature_json = $json->pretty->encode(\@features);
	my $pathwaymap_json = $json->pretty->encode(\@pathwaymap);
	my $spgenemap_json = $json->pretty->encode(\@spgenemap);

	open FH, ">$outfile.sequence.json"; 
	print FH $sequence_json;
	close FH;

	open FH, ">$outfile.feature.json"; 
	print FH $feature_json;
	close FH;

	open FH, ">$outfile.pathway.json"; 
	print FH $pathwaymap_json;
	close FH;

	open FH, ">$outfile.spgene.json"; 
	print FH $spgenemap_json;
	close FH;
	
	open FH, ">$outfile.genome.json"; 
	print FH "[".$genome_json."]";
	close FH;

}

sub postJson {
	
	`post.update.sh genome_sequence $outfile.sequence.json`;
	`post.update.sh genome_feature $outfile.feature.json`;
	`post.update.sh pathway $outfile.pathway.json`;
	`post.update.sh sp_gene $outfile.spgene.json`;
	`post.update.sh genome $outfile.genome.json`;
	
	`rm $outfile*.json`;

}

sub getGenomeInfo {

	my $genomeObj = Bio::SeqIO->new( -file   => "<$infile", -format => 'GenBank');

	my ($chromosomes, $plasmids, $contigs, $sequences, $cds, $genome_length, $gc_count);

	while (my $seqObj = $genomeObj->next_seq) {
		
		$sequences++;
		$genome_length += $seqObj->length;
		$gc_count += $seqObj->seq=~tr/GCgc//;
		
		$chromosomes++ if ($seqObj->desc=~/chromosome|complete genome/i); 
		$plasmids++ if ($seqObj->desc=~/plasmid/i); 
		$contigs++ if ($seqObj->desc=~/contig|scaffold/i || !($seqObj->desc=~/chromosome|complete genome|plasmid/i)); 

		next unless  $sequences == 1;
			
		$genome->{owner} = $owner;
		$genome->{public} = $public;

		$genome->{genome_id} = $genome_id;
		$genome->{assembly_accession} = $genome_id;
		$genome->{genome_name} = $seqObj->species->node_name;
		$genome->{common_name} = $seqObj->species->node_name;
		$genome->{common_name}=~s/\W+/_/g;
		$genome->{taxon_id}    =  $seqObj->species->ncbi_taxid;

		($genome->{taxon_lineage_ids}, $genome->{taxon_lineage_names}, $taxon_lineage_ranks)  = $solrh->getTaxonLineage($genome->{taxon_id});

		prepareTaxonomy($genome->{taxon_lineage_ids}) if $public;

		my $i=0;
		for my $rank (@{$taxon_lineage_ranks}){
			$genome->{kingdom} = $genome->{taxon_lineage_names}[$i] if $rank=~/kingdom/i;
			$genome->{phylum} = $genome->{taxon_lineage_names}[$i] if $rank=~/phylum/i;
			$genome->{class} = $genome->{taxon_lineage_names}[$i] if $rank=~/class/i;
			$genome->{order} = $genome->{taxon_lineage_names}[$i] if $rank=~/order/i;
			$genome->{family} = $genome->{taxon_lineage_names}[$i] if $rank=~/family/i;
			$genome->{genus} = $genome->{taxon_lineage_names}[$i] if $rank=~/genus/i;
			$genome->{species} = $genome->{taxon_lineage_names}[$i] if $rank=~/species/i;
			$i++;
		
		}

	}

	$genome->{chromosomes} = $chromosomes if $chromosomes;
	$genome->{plasmids} = $plasmids if $plasmids;
	$genome->{contigs} = $contigs if $contigs ;
	$genome->{sequences} = $sequences;
	$genome->{genome_length} = $genome_length;
	$genome->{gc_content} = sprintf("%.2f", ($gc_count*100/$genome_length));

	#print "metadata from genbank file\n";
	getMetadataFromGenBankFile();
	#print "metadata from prokaryotes.txt\n";
	getMetadataFromGenomeList();
	#print "metadata from bioproject\n";
	getMetadataFromBioProject($genome->{bioproject_accession}) if $genome->{bioproject_accession};
	#print "metadata from biosample\n";
	getMetadataFromBioSample($genome->{biosample_accession}) if $genome->{biosample_accession};


	$genomeObj->close();

}


sub getGenomeSequences {

	my $genomeObj = Bio::SeqIO->new( -file   => "<$infile", -format => 'GenBank');

	while (my $seqObj = $genomeObj->next_seq){

		my $sequence;
		
		$sequence->{owner} = $owner;
		$sequence->{public} = $public;
	
		$sequence->{genome_id} = $genome_id; # Get genome id from input directory

		$sequence->{genome_name} = $seqObj->species->node_name;
		$sequence->{taxon_id} = $seqObj->species->ncbi_taxid;

		#$sequence->{taxon_lineage_names} = join(';', reverse $seqObj->species->classification);
		#my @taxon_lineage = reverse $seqObj->species->classification;
		#$sequence->{taxon_lineage_names} = \@taxon_lineage;
	
		$sequence->{sequence_id} = $seqObj->accession_number;
		$sequence->{accession} = $seqObj->accession_number;
		$sequence->{gi} =	$seqObj->primary_id if $seqObj->primary_id=~/^\d*$/;

		my ($sequence_type) = $seqObj->desc=~/(chromosome|plasmid|contig|scaffold)/i;
		$sequence_type = "chromosome" if (!$sequence_type && $seqObj->desc=~/complete genome/i);
		$sequence_type = "contig" unless $sequence_type;
		$sequence->{sequence_type} = lc($sequence_type) if $sequence_type;
		$sequence->{topology} =	'circular' if $seqObj->is_circular;
		$sequence->{description} = $seqObj->desc;

		# Get the chromosole and plasmid number/name from corresponding source feature
		$sequence->{chromosome} = $1 if $sequence->{description}=~/chromosome (\S*)\s*,/i;
		$sequence->{plasmid} = $1 if $sequence->{description}=~/plasmid (\S*)\s*,/i;

		$sequence->{gc_content} = sprintf("%.2f", ($seqObj->seq=~tr/GCgc//)*100/$seqObj->length);
		$sequence->{length} = $seqObj->length;
		$sequence->{sequence} = $seqObj->seq if ((grep {$_=~/Bacteria|Archaea|Viruses/} @{$genome->{taxon_lineage_names}})); #subseq(1, 50);

		for my $date ($seqObj->get_dates){
			$sequence->{release_date} = strftime "%Y-%m-%dT%H:%M:%SZ", localtime str2time($date); 
		}
		$sequence->{version} =	$seqObj->seq_version if $seqObj->seq_version;


		push @sequences, $sequence;

	}

	$genomeObj->close();

}


sub getGenomeFeatures{
	
	my $genomeObj = Bio::SeqIO->new( -file   => "<$infile", -format => 'GenBank');

	while (my $seqObj = $genomeObj->next_seq){
			
		# my @taxon_lineage = reverse $seqObj->species->classification;
	
		for my $featObj ($seqObj->get_SeqFeatures){

			my $feature, $pathways, $ecpathways;
			my (@go, @ec_no, @ec, @pathway, @ecpathways, @spgenes, @uniprotkb_accns, @ids);

			$feature->{owner} = $owner;
			$feature->{public} = $public;
			
			$feature->{genome_id} = $genome_id;

			$feature->{genome_name} = $seqObj->species->node_name;
			$feature->{taxon_id} = $seqObj->species->ncbi_taxid;
			#$feature->{taxon_lineage_names} = join(';', reverse $seqObj->species->classification);

			$feature->{sequence_id} = $seqObj->accession_number;
			$feature->{accession} = $seqObj->accession_number;

			$feature->{annotation} = $annotation;
			$feature->{feature_type} = $featObj->primary_tag;
			$featureIndex->{$feature->{feature_type}}++;

			$feature->{start} = $featObj->start;
			$feature->{end} = $featObj->end;
			$feature->{strand} = ($featObj->strand==1)? '+':'-';
			$feature->{location} = $featObj->location->to_FTstring;

			my @segments;
			if ($featObj->location->isa('Bio::Location::SplitLocationI')){
				for my $location ($featObj->location->sub_Location){
					push @segments, $location->start."..". $location->end;
				}
			}else{
				push @segments, $featObj->start."..". $featObj->end;	
			}
			$feature->{segments} = \@segments;

			$feature->{pos_group} = "$feature->{sequence_id}:$feature->{end}:+" if $feature->{strand} eq '+';
			$feature->{pos_group} = "$feature->{sequence_id}:$feature->{start}:-" if $feature->{strand} eq '-';

			#$feature->{na_length} = length($featObj->spliced_seq->seq) if $featObj->spliced_seq->seq;
			#$feature->{na_sequence} = $featObj->spliced_seq->seq unless $featObj->primary_tag eq "source";


			for my $tag ($featObj->get_all_tags){

				for my $value ($featObj->get_tag_values($tag)){

					$feature->{feature_type} 	= 'pseudogene' if ($tag eq 'pseudo' && $feature->{feature_type} eq 'gene');

					$feature->{patric_id} = $1 if ($tag eq 'db_xref' && $value=~/SEED:(fig.*)/);
					
					#$feature->{refseq_locus_tag} 	= $value if ($tag eq 'locus_tag' && $annotation eq "RefSeq");
					$feature->{refseq_locus_tag} 	= $value if ($tag eq 'locus_tag');
					$feature->{refseq_locus_tag} 	= $1 if ($tag eq 'db_xref' && $value=~/Refseq_locus_tag:(.*)/i);
					
					$feature->{protein_id} 	= $value if ($tag eq 'protein_id');
					$feature->{protein_id} 	= $1 if ($tag eq 'db_xref' && $value=~/protein_id:(.*)/);
			
					$feature->{gene_id} = $1 if ($tag eq 'db_xref' && $value=~/^GeneID:(\d+)/);
					$feature->{gi} = $1 if ($tag eq 'db_xref' && $value=~/^GI:(\d+)/);

					$feature->{aa_sequence} = $value if ($tag eq 'translation');
					$feature->{aa_length} 	= length($value) if ($tag eq 'translation');
					$feature->{aa_sequence_md5} = md5_hex($value) if ($tag eq 'translation');
					
					$feature->{gene} = $value if ($tag eq 'gene');
					$feature->{product} = $value if ($tag eq 'product');
					#$feature->{product}=~s/\"/''/g;

					$feature->{figfam_id} 	= $value if ($tag eq 'FIGfam');
					
					push @ec_no, $value if ($tag eq 'EC_number');
					push @ec, $solrh->getEC($value) if ($tag eq 'EC_number');
					#push @uniprotkb_accns, $solrh->getUniprotkbAccns('GI', $1) if ($tag eq 'db_xref' && $value=~/^GI:(\d+)/ && $annotation eq 'PATRIC');

					if (($tag=~/note/i && $value=~/^go_/i) || $tag=~/^go_/i){
						$value=~s/^go_\W+: *//;
						foreach (split(/;/, $value)){
							my ($go_id) = $_=~/ *(GO:\d+)/;
							push @go, $solrh->getGO($go_id) if $go_id;
						}
					}elsif($tag=~/note/i){
						push @{$feature->{notes}}, $value if $value; 
					}else{

					}
		
					if ($feature->{feature_type}=~/source/i && $tag=~/note/i){
						push @{$genome->{comments}}, $value if $value;
					}
					
					push @spgenes, $value if ($tag=~/note/i && $value=~/^Similar to specialty gene/); 

					push @ids, $value if ($tag eq 'db_xref');
			
				}

			}

			$feature->{ec} = \@ec if scalar @ec;
			$feature->{go} = \@go if scalar @go;

			my ($pathways, $ecpathways) = $solrh->getPathways(@ec_no) if scalar @ec_no;
		 	$feature->{pathway} = $pathways if scalar @$pathways;
			#$feature->{ecpathway} = $ecpathways if scalar @$ecpathways; ## for debug, remove later
			push @pathwaymap, preparePathways($feature, $ecpathways);
	
			push @spgenemap, prepareSpGene($feature, $_) foreach(@spgenes);

			$feature->{uniprotkb_accession} = \@uniprotkb_accns if scalar @uniprotkb_accns;
					
			#push @ids, $solrh->getIDs(@uniprotkb_accns) if scalar @uniprotkb_accns;
			#$feature->{ids} = \@ids if scalar @ids;
			
			my $strand = ($feature->{strand} eq '+')? 'fwd':'rev';
			$feature->{feature_id}		=	"$annotation.$genome_id.$feature->{accession}.".
																	"$feature->{feature_type}.$feature->{start}.$feature->{end}.$strand";

	
			#$feature->{patric_id} = 'fig|'.$feature->{genome_id}.'.'.$feature->{feature_type}.'.'.$featureIndex->{$feature->{feature_type}} if ($annotation eq 'PATRIC' && !$feature->{patric_id});
			$feature->{patric_id} = 'fig|'.$feature->{genome_id}.'.'.$feature->{feature_type}.'.'.$featureIndex->{$feature->{feature_type}} unless $feature->{patric_id};
		
			#print "$feature->{feature_id}\t$feature->{patric_id}\n";
	
			push @features, $feature  unless ($feature->{feature_type} eq 'gene' && (grep {$_=~/Bacteria|Archaea/} @{$genome->{taxon_lineage_names}}));

		}

	}

	$genomeObj->close();

}


sub prepareSpGene {

		my ($feature, $note) = @_;
		my $spgene;

		my ($source, $source_id, $qcov, $scov, $identity, $evalue)
			= $note=~/Similar to specialty gene, DB:(.*), ID:\S*\|(.*), Query Coverage:(\d+), Subject Coverage:(\d+), Identity:(\d+), E-value:(.*)/;	

		my ($property, $locus_tag, $organism, $function, $classification, $pmid, $assertion) 
			= split /\t/, $solrh->getSpGeneInfo($source, $source_id) if ($source && $source_id);

		my ($qgenus) = $feature->{genome_name}=~/^(\S+)/;
		my ($qspecies) = $feature->{genome_name}=~/^(\S+ +\S+)/;
		my ($sgenus) = $organism=~/^(\S+)/;
		my ($sspecies) = $organism=~/^(\S+ +\S+)/;

		my ($same_genus, $same_species, $same_genome, $evidence); 

		$same_genus = 1 if ($qgenus eq $sgenus && $sgenus ne "");
		$same_species = 1 if ($qspecies eq $sspecies && $sspecies ne ""); 
		$same_genome = 1 if ($feature->{genome} eq $organism && $organism ne "") ;
	 	$evidence = $feature->{refseq_locus_tag} eq $locus_tag? 'Literature':'BLASTP';	

		$spgene->{owner} = $owner;
		$spgene->{public} = $public;
		
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
		$spgene->{classification} = $classification; 

		$spgene->{pmid} = $pmid; 
		$spgene->{assertion} = $assertion;

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


sub preparePathways {

	my ($feature, $ecpathways) = @_;
	my @pathways = ();

	foreach my $ecpathway (@$ecpathways){
		my $pathway;
		my ($ec_number, $ec_description, $pathway_id, $pathway_name, $pathway_class) = split /\t/, $ecpathway;
		
		$pathway->{owner} = $owner;
		$pathway->{public} = $public;

		$pathway->{genome_id} = $feature->{genome_id};
		$pathway->{genome_name} = $feature->{genome_name};
		$pathway->{taxon_id} = $feature->{taxon_id};

		$pathway->{sequence_id} = $feature->{sequence_id};
		$pathway->{accession} = $feature->{accession};
		
		$pathway->{annotation} = $feature->{annotation};
		
		$pathway->{feature_id} = $feature->{feature_id};
		$pathway->{patric_id} = $feature->{patric_id};
		$pathway->{alt_locus_tag} = $feature->{alt_locus_tag};
		$pathway->{refseq_locus_tag} = $feature->{refseq_locus_tag};
		
		$pathway->{gene} = $feature->{gene};
		$pathway->{product} = $feature->{product};
		
		$pathway->{ec_number} = $ec_number;
		$pathway->{ec_description} = $ec_description;
		
		$pathway->{pathway_id} = $pathway_id;
		$pathway->{pathway_name} = $pathway_name;
		$pathway->{pathway_class} = $pathway_class;
		
		$pathway->{genome_ec} = $feature->{genome_id}.'_'.$ec_number;
		$pathway->{pathway_ec} = $pathway_id.'_'.$ec_number;

		push @pathways, $pathway;

	}

	return @pathways;

}


sub getMetadataFromGenBankFile {

	print "Getting genome metadata from genbank file: $infile ...\n";

	open GB, "<$infile" || return "Can't open genbank file: $infile\n";
	my @gb = <GB>;
	my $gb = join "", @gb;
	close GB;

	my $strain = $2 if $gb=~/\/(strain|isolate)="([^"]*)"/;
	$strain =~s/\n */ /g;
	$genome->{strain} = $strain unless $strain=~/^ *(-|missing|na|n\/a|not available|not provided|not determined|nd|unknown) *$/i;

	$genome->{genome_name} .= " strain $genome->{strain}" if ($genome->{strain} && (not $genome->{genome_name}=~/$genome->{strain}/i));

	$genome->{geographic_location} = $1 if $gb=~/\/country="([^"]*)"/;
	$genome->{geographic_location} =~s/\n */ /g;
	$genome->{isolation_country} = $1 if $genome->{geographic_location}=~/^([^:]*):.*/;
	
	$genome->{host_name} = $1 if $gb=~/\/host="([^"]*)"/;
	$genome->{host_name} =~s/\n */ /g;
  
	$genome->{isolation_source} = $1 if $gb=~/\/isolation_source="([^"]*)"/;
	$genome->{isolation_source} =~s/\n */ /g;
	
	$genome->{collection_date} = $1 if $gb=~/\/collection_date="([^"]*)"/;
	$genome->{collection_date} =~s/\n */ /g;
	$genome->{collection_year} = $1 if $genome->{collection_date}=~/(\d\d\d\d)/;

	$genome->{culture_collection} = $1 if $gb=~/\/culture_collection="([^"]*)"/;
	$genome->{culture_collection} =~s/\n */ /g;
	
	$genome->{assembly_method} = $1 if $gb=~/Assembly Method\s*:: (.*)/;
	$genome->{assembly_method} =~s/\n */ /g;
  
	$genome->{sequencing_depth} = $1 if $gb=~/Genome Coverage\s*:: (.*)/;
	$genome->{sequencing_depth} =~s/\n */ /g;
  
	$genome->{sequencing_platform} = $1 if $gb=~/Sequencing Technology\s*:: (.*)/;
	$genome->{sequencing_platform} =~s/\n */ /g;

}


sub getMetadataFromGenomeList {

	my $assembly_accession = $genome->{genome_id};

	#print "$assembly_accession\n";

	print "Getting metadata from genome list\n";

	`wget -q "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"` unless (-e 'assembly_summary_genbank.txt');

	open PROK, "assembly_summary_genbank.txt" or "Can't open genome list file\n";

	while (my $line = <PROK>){

		chomp $line;			

		next unless ($line=~/$assembly_accession/i);

		my @a = split(/\t/, $line);

		#print join ("\n", @a);

		$genome->{bioproject_accession} = $a[1];
		$genome->{biosample_accession} = $a[2];

=pod
		#$genome->{domain} = ($a[4]=~/archae/i)? 'Archaea' : 'Bacteria';
		$genome->{genome_length} = $a[6] * 1000000 unless ($a[6] eq '-' || $genome->{genome_length});
		$genome->{gc_content} = $a[7] if (!$genome->{gc_content} && $a[7] ne '-');

		$genome->{refseq_accessions} = ($a[10] eq '-')? $a[8] : $a[8].','.$a[10];
		$genome->{genbank_accessions} = ($a[11] eq '-')? $a[9] : $a[9].','.$a[11];
		$genome->{genbank_accessions} = $a[12] if ($a[12] ne '-' && $genome->{genbank_accessions} eq '-');

		$genome->{chromosomes} = split(',', $a[9]) unless ($a[9] eq '-' || $genome->{chromosomes});
		$genome->{plasmids} = split(',', $a[11]) unless ($a[11] eq '-' || $genome->{plasmids});
		$genome->{contigs} = $a[13] unless ($a[13] eq '-' || $genome->{contigs});

		$genome->{completion_date} = $a[16];
		$genome->{completion_date}=~s/\//-/g;
		$genome->{completion_date}=~s/$/T00:00:00Z/g;

		$genome->{genome_status} = ( ($a[18]=~/complete|chromosome|gapless/i)? 'complete': ($a[18]=~/contig|scaffold/i)? 'WGS' : $a[18]);
		$genome->{sequencing_centers} = $a[19];
			
		$genome->{biosample_accession} = $a[20];
		$genome->{assembly_accession} = $a[21];
		$genome->{reference_genome} = $a[22];
		$genome->{publication} = $a[24];

		last if ($line=~/$genome_name.*$bioproject_accession.*$assembly_accession/);		
=cut
		last;

	}

	close PROK;

}



sub getMetadataFromBioProject {
	
	my($bioproject_accn) = @_;

	print "Getting genome metadata from BioProject: $bioproject_accn...\n";

  my $xml = `wget -q -O - "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=$bioproject_accn"`;
  $xml=~s/\n//;
  my ($bioproject_id) = $xml=~/<Id>(\d+)<\/Id>/;

	`wget -q -O "$outfile.bioproject.xml" "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=bioproject&retmode=xml&id=$bioproject_id"`;

	return unless -f "$outfile.bioproject.xml";
	
	my ($projID, $projDescription, $subgroup, $organism, $description, $var, $serovar, $biovar, $pathovar, $strain, $cultureCollection, $typeStrain);
	my ($isolateComment, $source, $month, $year, $country, $method, $person, $epidemic, $location, $altitude, $depth);
	my ($hostName, $hostGender, $hostAge, $hostHealth);
	my ($publication, $taxonID, $epidemiology);


	my $xml = XMLin("$outfile.bioproject.xml", ForceArray => ["Organization"]);
	my $root = $xml->{DocumentSummary};

	#print Dumper $xml;
	
	$organism = $root->{Project}->{ProjectType}->{ProjectTypeSubmission}->{Target}->{Organism};

	my $organization = $root->{Submission}->{Description}->{Organization};	
	if (ref($organization->[0]->{Name}) eq "HASH"){
		$genome->{sequencing_centers} = $organization->[0]->{Name}->{content};
	}else{
		$genome->{sequencing_centers} = $organization->[0]->{Name};
	}

	$genome->{strain} = $organism->{Strain} if $organism->{Strain} && not $genome->{strain};
	$genome->{genome_name} .= " strain $genome->{strain}" if ($genome->{strain} && (not $genome->{genome_name}=~/$genome->{strain}/) );
	
	$genome->{disease} = $organism->{BiologicalProperties}->{Phenotype}->{Disease};
	$genome->{temperature_range} = $1 if $organism->{BiologicalProperties}->{Environment}->{TemperatureRange}=~/^e*(.*)/;
	$genome->{optimal_temperature} = $organism->{BiologicalProperties}->{Environment}->{OptimumTemperature};
	$genome->{oxygen_requirement} = $1 if $organism->{BiologicalProperties}->{Environment}->{OxygenReq}=~/^e*(.*)/;
	$genome->{habitat} = $1 if $organism->{BiologicalProperties}->{Environment}->{Habitat}=~/^e*(.*)/;	
	$genome->{cell_shape} = $1 if $organism->{BiologicalProperties}->{Morphology}->{Shape}=~/^e*(.*)/;
	$genome->{motility} = $1 if $organism->{BiologicalProperties}->{Morphology}->{Motility}=~/^e*(.*)/;
	$genome->{gram_stain} = $1 if $organism->{BiologicalProperties}->{Morphology}->{Gram}=~/^e*(.*)/;

	if (ref($root->{Project}->{ProjectDescr}->{Publication}->{id}) eq "HASH"){
		foreach my $pmid (keys %{$root->{Project}->{ProjectDescr}->{Publication}}){
			$genome->{publication} .= ",$pmid" unless $genome->{publication}=~/$pmid/; 
		}
	}else{
		my $pmid = $root->{Project}->{ProjectDescr}->{Publication}->{id};
	}
	$genome->{publication}=~s/^,|,$//g;

	$description = $root->{Project}->{ProjectDescr}->{Description};

	$description=~s/&lt;\/*.&gt;|&#x0D;|<\/*.>//g;
	$description=~s/&lt;.*?&gt;|<.*?>//g;
	$description=~s/^ *| *$//g;
	$description=~s/\t+/ /g;
	$description=~s/Dr\.\s*/Dr /g;

	$typeStrain="Yes" if($description=~/type str/i);
		
	my($var1, $var2) = "$organism $description"=~/(serovar|sv\.|sv|serotype|biovar|bv\.|bv|biotype|pathovar|pv\.|pv|pathotype)\s*([\w-]*)/i;
	$var = "$var1 $var2" if($var2 && !($var2=~/^of|and$/i));
	$var =~s/serovar|sv\.|sv|serotype/serovar/i;
	$var =~s/biovar|bv\.|bv|biotype/biovar/i;
	$var =~s/pathovar|pv\.|pv|pathotype/pathovar/i;
	$var =~s/^\s*$//;
	$serovar = $var if($var=~/serovar/);
	$biovar = $var if($var=~/biovar/);
	$pathovar = $var if($var=~/pathovar/);

	my($cc1, $cc2) = $description=~/(ATCC|NCTC|CCUG|DSM|LMG|CIP|NCIB|BCCM|NRRL)\s*([\w-]*)/;
	$cultureCollection = "$cc1 $cc2" if ($cc1 && $cc2);
	
	#($isolateComment) = $description=~/(isolated\s*.*?\S\S\.)/i;
	my ($com1, $com2) = $description=~/(isolated|isolate from|isolate obtained from|derived|is from|was from|came from|recovered from)\s+(.*?\S\S\.)/i;
	$isolateComment = "$com1 $com2";
	$isolateComment =~s/\s*$|\s*\.$//;
	$isolateComment =~s/&lt;|&gt;//g;

	$source = $1 if $isolateComment=~/from\s*(.*)/;
	$source =~s/^(a|an)\s+//;
	$source =~s/( in | by | and will | and is | and is |\s*\.).*//g;
	
	($month, $year) = $isolateComment=~/in\s*([a-zA-Z]*)\s*,*\s+(\d\d\d\d)/;

	($method) = $isolateComment=~/isolated by ([a-z].*)/;
	$method =~s/ (from|in) .*//;

	($person) = $isolateComment=~/isolated by ([A-Z].*)/;
	$person =~s/ (at|in|during) .*//;

	($epidemiology) = $isolateComment=~/(outbreak|epidemic|pandemic)/i;

	$isolateComment=~/(^|\s)(in|in the|at|at the|near|near the)\s([A-Z]+.*)/;	
	$location = $3;
	$location =~s/\..*//;
	
	$location =~s/\W+/ /g;
	$location =~s/^\s*|\s*$|\s\s//g;
	$location =~s/ [a-z0-9]+.*//g;

	my ($num, $unit);
	if ($isolateComment=~/depth of/i){
		($num, $unit) = $isolateComment=~/depth of (\d+\s*[\.\-to]*\s*\d*)\s*(\w+)/; 
		$depth = "$num $unit";
	}elsif($isolateComment=~/below the surface/i){
		($num, $unit) = $isolateComment=~/(\d+\s*[\.\-to]*\s*\d*)\s*(\w+) below the surface/;
		$depth = "$num $unit";
		$depth =~s/\s$//;
	}

	$hostGender = $1 if $isolateComment=~/(male|female)/i;	# man|woman	
	$hostAge = "$1 $2" if $isolateComment=~/(\d+).(year|years|month|months|day|days).old/; # month|months|day|days	

	$hostHealth = $2 if $isolateComment=~/(patient with |suffering from )(.*)/i; # suffered|diagnosed with 
	$hostHealth =~s/^(a|an)\s//;
	$hostHealth =~s/\s+(and is|and will be|and has|and was|by|in|and contains)\s+.*//;

	$hostName = $1 if $isolateComment=~/ (pig|sheep|goat|dog|cat |cattle|chicken|cow|mouse|rat|buffalo|tick|mosquito)/i;
	$hostName ="Human, Homo sapiens" if ($isolateComment=~/ (human|man|woman|infant|child|patient|homo sapiens)/i);


	# Organism Info

	$genome->{serovar} = $serovar if $serovar;
	$genome->{biovar} = $biovar if $biovar;
	$genome->{pathovar} = $pathovar if $pathovar;
	$genome->{type_strain} = $typeStrain if $typeStrain;
	$genome->{culture_collection} = $cultureCollection if $cultureCollection;
	push @{$genome->{comments}}, $description if $description;
	$genome->{publication} = $publication if $publication;


	# Isolate / Environmental Metadata
	
	#$genome->{isolation_site} = ""; #$source;
	$genome->{isolation_source} = $source if $source;
	$genome->{isolation_comments} = $isolateComment if $isolateComment;
	$genome->{collection_date} = $year if $year;
	#$genome->{isolation_country} = ""; #$location;
	$genome->{geographic_location} = $location if $location;
	#$genome->{latitude} = "";
	#$genome->{longitude} = "";
	#$genome->{altitude} = "";
	#$genome->{depth} = $depth;

	# Host Metadata
	
	$genome->{host_name} = $hostName if $hostName;
	$genome->{host_gender} = $hostGender if $hostGender;
	$genome->{host_age} = $hostAge if $hostAge ;
	$genome->{host_health} = $hostHealth if $hostHealth;
	#$genome->{body_sample_site} = "";
	#$genome->{body_sample_substitute} = "";

}


sub getMetadataFromBioSample {

	my($biosample_accn) = @_;

	print "Getting genome metadata from BioSample: $biosample_accn ...\n";
  
	my $xml = `wget -q -O - "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term=$biosample_accn"`;
  $xml=~s/\n//;
  my ($biosample_id) = $xml=~/<Id>(\d+)<\/Id>/;

	`wget -q -O "$outfile.biosample.xml" "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=$biosample_id"`;

	return unless -f "$outfile.biosample.xml";
	
	my $xml = XMLin("$outfile.biosample.xml", ForceArray => ["Row"]);

	return unless ref $xml->{BioSample}->{Attributes}->{Attribute} eq 'ARRAY';

	foreach my $attrib (@{$xml->{BioSample}->{Attributes}->{Attribute}}){
	
		my $attrib_name = $attrib->{harmonized_name};
		my $attrib_value = $attrib->{content};

		if ($attrib_name=~/lat_lon/i){ 
			my ($latitude, $longitude) = $attrib=~/(\d+\.*\d* [NS])\s+(\d+\.*\d* [EW])/;
			$genome->{latitude} = $latitude;
			$genome->{longitude} = $longitude;	
		}else{
			my $patric_attrib_name = biosample2patricAttrib($attrib_name);
			next unless ($patric_attrib_name && $attrib_value && not $attrib_value=~/^ *(-|missing|na|n\/a|not available|not provided|not determined|nd|unknown) *$/i);
		
			if ($patric_attrib_name=~/other|additional/i){
				my ($attrib1, $attrib2) = $patric_attrib_name=~/^([^:]*):(.*)/;
				push @{$genome->{$attrib1}}, "$attrib2:$attrib_value";
			}elsif ($patric_attrib_name=~/comments/i){
				push @{$genome->{$patric_attrib_name}}, $attrib_value;
			}else{
				$genome->{$patric_attrib_name} = $attrib_value;
			} 
		}

	}

	$genome->{genome_name} .= " strain $genome->{strain}" if ($genome->{strain} && (not $genome->{genome_name}=~/$genome->{strain}/) ) ;
	
	# parse AMR metadata

	return unless $xml->{BioSample}->{Description}->{Comment}->{Table}->{class}=~/Antibiogram/i;

	foreach my $row (@{$xml->{BioSample}->{Description}->{Comment}->{Table}->{Body}->{Row}}){
	
		my @amr1 = @{$row->{Cell}};

		my $amr;

		$amr->{owner} = $genome->{owner};
		$amr->{public} = $public;
		$amr->{genome_id} = $genome->{genome_id};
		$amr->{genome_name} = $genome->{genome_name};
		$amr->{taxon_id} = $genome->{taxon_id};

		$amr->{antibiotic} = lc $amr1[0]; 	
		$amr->{resistant_phenotype} = ucfirst $amr1[1];	
		$amr->{measurement_sign} = $amr1[2] unless ref $amr1[2] eq ref {};	
		$amr->{measurement_value} = $amr1[3] unless ref $amr1[3] eq ref {};
		$amr->{measurement} = $amr->{measurement_sign}.$amr->{measurement_value}; 	
		$amr->{measurement_unit} = $amr1[4] unless ref $amr1[4] eq ref {}; 	
		$amr->{laboratory_typing_method} = $amr1[5] unless ref $amr1[5] eq ref {}; 	
		$amr->{laboratory_typing_platform} = $amr1[6] unless ref $amr1[6] eq ref {}; 	
		$amr->{vendor} = $amr1[7] unless ref $amr1[7] eq ref {}; 	
		$amr->{laboratory_typing_method_version} = $amr1[8] unless ref $amr1[8] eq ref {}; 	
		$amr->{testing_standard} = $amr1[9] unless ref $amr1[9] eq ref {}; 	

		push @{$genome->{antimicrobial_resistance}}, ucfirst $amr->{resistant_phenotype} unless (grep {$_ eq ucfirst $amr->{resistant_phenotype}} @{$genome->{antimicrobial_resistance}}); 
		$genome->{antimicrobial_resistance_evidence} = "AMR Panel";
	
		push @genome_amr, $amr;

	}

}

sub biosample2patricAttrib{

	my ($attribute) = @_;

	my %biosample2patric = (

		"altitude" => "altitude",
		"biomaterial_provider" => "additional_metadata:biomaterial_provider",
		"collected_by" => "additional_metadata:collected_by",
		"collection_date" => "collection_date",			
		"culture_collection" => "culture_collection", 
		"depth" => "depth", 
		"description" => "comments", 
		"env_biome" => "other_environmental:env_biome", 
		"genotype" => "other_typing:genotype", 
		"geo_loc_name" => "geographic_location", 
		"host" => "host_name", 
		"host_age" => "host_age", 
		"host_description" => "other_clinical:host_description", 
		"host_disease" => "host_health",
		"host_disease_outcome" => "other_clinical:host_disease_outcome",
		"host_disease_stage" => "other_clinical:host_disease_stage",
		"host_health_state" => "other_clinical:host_health_state",
		"host_sex" => "host_gender",
		"host_subject_id" => "other_clinical:host_subject_id", 
		"host_tissue_sampled" => "isolation_source",
		"identified_by" => "additional_metadata:identified_by",
		"isolate" => "additional_metadata:isolate", 
		"isolation_source" => "isolation_source", 
		"lab_host" => "additional_metadata:lab_host", 
		"lat_lon" => "other_environmental:lat_lon", 
		"mating_type" => "additional_metadata:mating_type", 
		"organism" => "",
		"passage_history" => "additional_metadata:passage_history",
		"pathotype" => "pathovar",
		"sample_name" => "",
		"sample_title" => "",
		"sample_type" => "additional_metadata:sample_type",
		"samp_size" => "",
		"serotype" => "serovar",
		"serovar" => "serovar",
		"specimen_voucher" => "additional_metadata:specimen_voucher", 
		"strain" => "strain",
		"isolate" => "strain",
		"subgroup" => "",
		"subtype" => "",
		"temp" => "other_environmental:temperature"
	);

	return $biosample2patric{$attribute} if $attribute;	

}

sub prepareTaxonomy {
	
	"Preparing taxonomy update ...\n";

	my ($taxon_lineage_ids) = @_;

	foreach my $taxon_id (@{$taxon_lineage_ids}){
	
		my $taxon;	
		$taxon->{taxon_id} = $taxon_id;
		$taxon->{genomes}->{inc} = 1;		

		push @taxonomy, $taxon

	}

}


