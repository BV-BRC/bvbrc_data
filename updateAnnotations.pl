#!/usr/bin/env perl 

use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long::Descriptive;
use JSON;
use Data::Dumper;

use lib "$Bin";
use SolrAPI;

my $solrServer = $ENV{PATRIC_SOLR};
my $solrFormat="&wt=csv&csv.separator=%09&csv.mv.separator=;";

my $solrh = SolrAPI->new($ENV{PATRIC_DATA_API}, $ENV{PATRIC_REFERENCE_DATA});
my $json = JSON->new->allow_nonref;


my ($opt, $usage) = describe_options(
		"%c %o",
		[],
    ["genome_list=s", "File containing list of annotation files"], 
    ["commit=s", "Commit updates to Solr, true|false", { default => "false"}],
		[],
    ["help|h", "Print usage message and exit"]
);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless $opt->genome_list;


# Get reference data
my $ecRef = $solrh->getECRef();
my $pathwayRef = $solrh->getPathwayRef();
my $spgeneRef = $solrh->getSpGeneRef();


# process functions and prepare update files
processAnnotationFiles($opt->genome_list);



sub processAnnotationFiles {

 my ($genome_list) = @_;
 open LIST, "$genome_list" or die "Can't open genome_list file: $genome_list!!\n";

 while (my $file_name = <LIST>) {

	chomp $file_name;
	next unless $file_name;

 	my $genome_id = $file_name;
	$genome_id =~s/.*\///;		
	print "Processing $genome_id\n";

	# global arrays to record all the updates
	my ($fids, $spgids, $pathway_ids, $feature_records);
	my @update_features;
	my @update_spgenes;
	my @delete_pathways;
	my @add_pathways;

	# new genome, get primary identifiers and existing features
	$fids = getFeatureIDs($genome_id);
	$spgids = getSPGeneIDs($genome_id);
	$pathway_ids = getPathwayIDs($genome_id);
	$feature_records = getFeatures($genome_id);

	# process each genome files
	if (-f $file_name){ 
		open GENOME, "$file_name" or next;
	}else {
		print "\t$file_name doesn't exist!!\n";
	}  

	while (my $row = <GENOME>){

		chomp $row;

		my ($patric_id, $update_code, $old_product, $old_pgfam_id, $old_plfam_id, 
				$product, $pgfam_id, $plfam_id) = split /\t/, $row;

		if ($plfam_id){	
			my ($plfam_prefix, $plfam_idx) = $plfam_id=~/(PLF_\d+)_(\d+)/;
			$plfam_id = $plfam_prefix."_".sprintf("%08d", $plfam_idx);
		}

		# no changes in function or families
		next if $update_code eq "000";
		next if ($old_product eq $product && $old_pgfam_id eq $pgfam_id && $old_plfam_id eq $plfam_id);

		# process entry, compare functions and prepare updates 

		my ($feature, @ec_no, @ec, @go, @pathways, @ecpathways);	
		my $feature_id = $fids->{$patric_id};

		next unless $feature_id; 
	
		# update genome_feature
		$feature->{feature_id} = $feature_id;
		$feature->{pgfam_id}->{set} = $pgfam_id; #unless $old_pgfam_id eq $pgfam_id;
		$feature->{plfam_id}->{set} = $plfam_id; # unless $old_plfam_id eq $plfam_id;

		# if product:
		# 	- update product, ec, go and pathways
		# 	- add new pathways
			
		if($product && $product ne $old_product){

			# New function is assigned by the recall
			$product=~s/^\s*|\s*$//g;
			@ec_no = $product=~/\( *EC[: ]([\d\-\.]+) *\)/g;		
    
			# if product has EC number, get corresponding EC,GO and pathways
			foreach my $ec_number (@ec_no){

      	my $ec_description = $ecRef->{$ec_number}->{ec_description};
				next unless $ec_description;
      	push @ec, $ec_number.'|'.$ec_description unless (grep {$_ eq $ec_number.'|'.$ec_description} @ec);

      	foreach my $go_term (@{$ecRef->{$ec_number}->{go}}){
       	 push @go, $go_term unless (grep {$_ eq $go_term} @go);
      	}

      	foreach my $pathway (@{$pathwayRef->{$ec_number}->{pathway}}){
       	 my ($pathway_id, $pathway_name, $pathway_class) = split(/\t/, $pathway);
        	push @pathways, $pathway_id.'|'.$pathway_name unless (grep {$_ eq $pathway_id.'|'.$pathway_name} @pathways);
					my $ecpathway = "$ec_number\t$ec_description\t$pathway_id\t$pathway_name\t$pathway_class";
					push @ecpathways, $ecpathway unless (grep {$_ eq $ecpathway} @ecpathways);
      	}

			}

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
			push @delete_pathways, $patric_id if $pathway_ids->{$patric_id};	

			# add new pathways
			$feature_records->{$patric_id}->{product}=$product;
			push @add_pathways, preparePathways($feature_records->{$patric_id}, \@ecpathways) if scalar @ecpathways;

		}else{

			# No change in function from recall
			# 	- leave the existing function, ec, go, and pathways
			# 	- just need to update the family assignment

			push @update_features, $feature;	
		
		}

  }

	close GENOME;
	prepareUpdateFiles($genome_id, \@update_features, \@update_spgenes, \@delete_pathways, \@add_pathways);
	postUpdateFiles($genome_id) if $opt->commit eq "true";

 }

  close LIST;

}



sub prepareUpdateFiles {

	my($genome_id, $update_features, $update_spgenes, $delete_pathways, $add_pathways) = @_;

	print "\tPrepare $genome_id.genome_feature.json\n";
	my $feature_json = $json->pretty->encode($update_features);
	open GF, ">$genome_id.genome_feature.json";
	print GF "$feature_json";
	close GF;

	print "\tPrepare $genome_id.sp_gene.json\n";
	my $spgenes_json = $json->pretty->encode($update_spgenes);
	open SPG, ">$genome_id.sp_gene.json";
	print SPG "$spgenes_json";
	close SPG;

	print "\tPrepare $genome_id.pathway.json\n";
	my $pathway_json = $json->pretty->encode($add_pathways);
	open PATH, ">$genome_id.pathway.json";
	print PATH "$pathway_json";
	close PATH;
	
	print "\tPrepare $genome_id.delete_pathway.sh\n";
	open DEL_PATH, ">$genome_id.delete_pathway.sh";
	my $patric_ids = scalar @$delete_pathways? join(" OR ", @$delete_pathways) : "12345";
	my $query_str = "<delete><query>patric_id:($patric_ids)</query></delete>"; 
	print DEL_PATH "curl -X POST \"$solrServer/pathway/update?commit=true\" -d \"$query_str\"";
	close DEL_PATH;

}

sub postUpdateFiles {

	my($genome_id) = @_;

	
	print "\tPost $genome_id.delete_pathway.sh\n";
	`chmod 755 $genome_id.delete_pathway.sh`;
	`./$genome_id.delete_pathway.sh`;
	#`rm $genome_id.delete_pathway.sh`;
	
	print "\tPost $genome_id.genome_feature.json\n";
	`post.update.sh genome_feature $genome_id.genome_feature.json`;
	#`rm $genome_id.genome_feature.json`;
	
	print "\tPost $genome_id.pathway.json\n";
	`post.update.sh pathway $genome_id.pathway.json`;
	#`rm $genome_id.pathway.json`;

	print "\tPost $genome_id.sp_gene.json\n";
	`post.update.sh sp_gene $genome_id.sp_gene.json`;
	#`rm $genome_id.sp_gene.json`;

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
