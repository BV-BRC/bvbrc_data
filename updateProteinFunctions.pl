#!/usr/bin/perl 

use FindBin qw($Bin);
use JSON;
use Data::Dumper;

use lib "$Bin";
use SolrAPI;

my $solrServer = $ENV{PATRIC_SOLR};
my $data_api_url = $ENV{PATRIC_DATA_API};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my $solrh = SolrAPI->new($data_api_url);
my $json = JSON->new->allow_nonref;

# Get reference data
my $ecRef = $solrh->getECRef();
my $pathwayRef = $solrh->getPathwayRef();

my $infile = $ARGV[0];

# process functions and prepare update files
processFunctions($infile);



sub processFunctions {

	my ($infile) = @_;

	open IN, "$infile";

	my $prev_genome_id = "";

  while (my $row = <IN>){

		chomp $row;
			print "$genome_id\n";

		my ($patric_id, $pgfam_id, $product) = $row=~/(.*)\t(.*)\t(.*)/;

		next unless ($pgfam_id || $product);

		my ($genome_id) = $patric_id=~/fig\|(.*).peg/;

		if ($genome_id ne $prev_genome_id || !$prev_genome_id){

			# post updates for the last genome
			prepareUpdateFiles($prev_genome_id, \@update_features, \@update_spgenes, \@delete_pathways, \@add_pathways) if $prev_genome_id || last;

			print "$genome_id\n";

			# global arrays to record all the updates
			my @update_features;
			my @update_spgenes;
			my @delete_pathways;
			my @add_pathways;

			# new genome, get primary identifiers and existing features
			$fids = getFeatureIDs($genome_id);
			$spgids = getSPGeneIDs($genome_id);
			$pathway_ids = getPathwayIDs($genome_id);
			$feature_records = getFeatures($genome_id);

		}

		# process entry, compare functions and prepare updates 

		my ($feature, @ec_no, @ec, @go, @pathways, @ecpathways);	
		my $feature_id = $fids->{$patric_id};

		next unless $feature_id; 
	
		$product=~s/^\s*|\s*$//g;
		@ec_no = $product=~/\( *EC[: ]([\d-\.]+) *\)/g;		
    
		foreach my $ec_number (@ec_no){

      my $ec_description = $ecRef->{$ec_number}->{ec_description};
      push @ec, $ec_number.'|'.$ec_description unless (grep {$_ eq $ec_number.'|'.$ec_description} @ec);

      foreach my $go_term (@{$ecRef->{$ec_number}->{go}}){
        push @go, $go_term unless (grep {$_ eq $go_term} @go);
      }

      foreach my $pathway (@{$pathwayRef->{$ec_number}->{pathway}}){
        my ($pathway_id, $pathway_name, $pathway_class) = split(/\t/, $pathway);
        push @pathways, $pathway_id.'|'.$pathway_name unless (grep {$_ eq $pathway_id.'|'.$pathway_name} @pathways);
        $ecpathway = "$ec_number\t$ec_description\t$pathway_id\t$pathway_name\t$pathway_class";
        push @ecpathways, $ecpathway unless (grep {$_ eq $ecpathway} @ecpathways);
      }

    }

		# update genome_feature
		$feature->{feature_id} = $feature_id;
		$feature->{pgfam_id}->{set} = $pgfam_id;

		# if product:
		# 	- update product, ec, go and pathways
		# 	- add new pathways
			
		if($product){

			# New function is assigned by function recall

			# update genome_feature
			$feature->{product}->{set} = $product;
    	$feature->{ec}->{set} = \@ec; # if scalar @ec;
    	$feature->{go}->{set} = \@go; # if scalar @go;
    	$feature->{pathway}->{set} = \@pathways; # if scalar @pathways;
		
			push @update_features, $feature;

			# update sp_genes
			foreach my $id (@{$spgids->{$patric_id}}){
				my $spg;
				$spg->{id}=$id;
				$spg->{product}->{set}=$product;
				push @update_spgenes, $spg;
			}

			# delete existing pathways for a given gene
			push @delete_pathways, "curl \"$solrServer/pathway/update?stream.body=<delete><query>patric_id:$patric_id</query></delete>&commit=true\"\n" if $pathway_ids->{$patric_id};	

			# add new pathways
			$feature_records->{$patric_id}->{product}=$product;
			push @add_pathways, preparePathways($feature_records->{$patric_id}, \@ecpathways) if scalar @ecpathways;

		}else{

			# new function is not called by kmers
			# 	- leave the existing function, ec, go, and pathways
			# 	- just need to update the family assignment

			push @update_features, $feature;	
		
		}

			$prev_genome_id = $genome_id;

  }

	prepareUpdateFiles($prev_genome_id, \@update_features, \@update_spgenes, \@delete_pathways, \@add_pathways);

  close IN;

}

sub prepareUpdateFiles{

	my($genome_id, $update_features, $update_spgenes, $delete_pathways, $add_pathways) = @_;

	my $feature_json = $json->pretty->encode($update_features);
	open GF, ">$genome_id.genome_feature.json";
	print GF "$feature_json";
	close GF;

	my $spgenes_json = $json->pretty->encode($update_spgenes);
	open SPG, ">$genome_id.sp_gene.json";
	print SPG "$spgenes_json";
	close SPG;

	open DEL_PATH, ">$genome_id.delete_pathway.sh";
	print DEL_PATH join("", @$delete_pathways);
	close DEL_PATH;

	my $pathway_json = $json->pretty->encode($add_pathways);
	open PATH, ">$genome_id.pathway.json";
	print PATH "$pathway_json";
	close PATH;

}


sub getFeatureIDs {

	my ($genome_id) = @_;

	my $core = "/genome_feature";
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


sub getFeatures {

	my ($genome_id) = @_;

	my $core = "/genome_feature";
	my $query = "/select?q=annotation:PATRIC AND feature_type:CDS AND genome_id:$genome_id";
	my $fields = "&fl=genome_id,genome_name,taxon_id,sequence_id,accession,annotation,feature_id,patric_id,alt_locus_tag,refseq_locus_tag,gene,product,owner,public";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrFormat="&wt=json";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my $result = `wget -q -O - "$solrQuery"`;

	my $resultObj = decode_json($result);
	my $features;

	foreach my $record(@{$resultObj->{'response'}->{'docs'}}){
		$features->{$record->{patric_id}} = $record;

	}

	return $features;

}


sub getPathwayIDs {

	my ($genome_id) = @_;

	my $core = "/pathway";
	my $query = "/select?q=annotation:PATRIC AND genome_id:$genome_id";
	my $fields = "&fl=patric_id,id";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my @result = `wget -q -O - "$solrQuery" | grep -v patric_id`;

	my %pathway_ids = ();
	foreach my $record (@result){
		my ($patric_id, $id) = $record=~/(.*)\t(.*)/;
		push @{$pathway_ids{$patric_id}}, $id;
	}

	return \%pathway_ids;

}



sub getSPGeneIDs {

	my ($genome_id) = @_;

	my $core = "/sp_gene";
	my $query = "/select?q=genome_id:$genome_id";
	my $fields = "&fl=patric_id,id";
	my $rows = "&rows=1000000";
	my $sort = "";
	my $solrQuery = $solrServer.$core.$query.$fields.$rows.$sort.$solrFormat;

	my @result = `wget -q -O - "$solrQuery" | grep -v patric_id`;

	my %spgids = ();
	foreach my $record (@result){
		my ($patric_id, $id) = $record=~/(.*)\t(.*)/;
		push @{$spgids{$patric_id}}, $id;
	}

	return \%spgids;

}


sub preparePathways {

  my ($feature, $ecpathways) = @_;
  my @pathways = ();

  foreach my $ecpathway (@$ecpathways){
    my $pathway;
    my ($ec_number, $ec_description, $pathway_id, $pathway_name, $pathway_class) = split /\t/, $ecpathway;

    $pathway->{owner} = $feature->{owner};
    $pathway->{public} = $feature->{public};

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
