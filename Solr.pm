package Solr;
use strict;
use JSON;
use Data::Dumper;


my %params = ();
my $solr;

my $json = JSON->new->allow_nonref;

sub new{

	my ($class, $solr_db) = @_;
	my $self = {};
	bless $self, $class;
	$solr = $ENV{PATRIC_SOLR};
	return $self;

}

sub getTaxonLineage {

	my ($self, $taxon_id) = @_;
	my ($lineage_ids, $lineage_names, $lineage_ranks);

	my $core = "/taxonomy";
	my $query = "/select?q=taxon_id:$taxon_id";
	my $fields = "&fl=lineage_ids,lineage_names,lineage_ranks";
	my $format = "&wt=json&rows=10000";


	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;
	
	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		$lineage_ids = $record->{lineage_ids};
		$lineage_names = $record->{lineage_names}; 
		$lineage_ranks = $record->{lineage_ranks}; 
	}

	return ($lineage_ids, $lineage_names, $lineage_ranks);

}


sub getEC {

	my ($self, $ec_no) = @_;
	my $ec;

	my $core = "/enzyme_class_ref";
	my $query = "/select?q=ec_number:$ec_no";
	my $fields = "&fl=ec_number,ec_description";
	my $format = "&wt=json";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;
	
	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		$ec = $record->{ec_number}.'|'.$record->{ec_description};
	}

	return $ec;

}


sub getGO {

	my ($self, $go_id) = @_;
	my $go;

	my $core = "/gene_ontology_ref";
	my $query = "/select?q=go_id:\\\"$go_id\\\"";
	my $fields = "&fl=go_id,go_name";
	my $format = "&wt=json";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;

	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		$go = $record->{go_id}.'|'.$record->{go_name};
	}

	return $go;

}

sub getECGO {

	my ($self, @ec_no) = @_;
	my @ec = ();
	my @go = ();

	my $core = "/enzyme_class_ref";
	my $query = "/select?q=ec_number:(".join(" OR ", @ec_no).")";
	my $fields = "&fl=ec_number,ec_description,go";
	my $format = "&wt=json&rows=10000";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;

	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		push @ec, $record->{ec_number}.'|'.$record->{ec_description};
		foreach my $go (@{$record->{go}}){
			push @go, $go unless (grep {$_ eq $go} @go);
		}
	}

	return (\@ec, \@go);

}


sub getPathways {

	my ($self, @ec_no) = @_;
	my @pathways = ();
	my @ecpathways = ();

	my $core = "/pathway_ref";
	my $query = "/select?q=ec_number:(".join(" OR ", @ec_no).")";
	my $fields = "&fl=ec_number,ec_description,pathway_id,pathway_name,pathway_class";
	my $format = "&wt=json&rows=10000";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;

	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		my $pathway = $record->{pathway_id}.'|'.$record->{pathway_name};
		my $ecpathway = "$record->{ec_number}\t$record->{ec_description}\t$record->{pathway_id}\t$record->{pathway_name}\t$record->{pathway_class}";
		push @pathways, $pathway unless (grep {$_ eq $pathway} @pathways);
		push @ecpathways, $ecpathway unless (grep {$_ eq $ecpathway} @ecpathways);
	}

	return (\@pathways, \@ecpathways);

}


sub getSpGeneInfo {

	my ($self, $source, $source_id) = @_;
	my $spgeneinfo;

	my $core = "/sp_gene_ref";
	my $query = "/select?q=source:$source AND source_id:$source_id";
	my $fields = "&fl=property,locus_tag,organism,function,classification,pmid,assertion";
	my $format = "&wt=json&rows=10000";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;

	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		$spgeneinfo = "$record->{property}\t$record->{locus_tag}\t$record->{organism}"
									."\t$record->{function}\t";
		$spgeneinfo .= join(',', @{$record->{classification}}) if $record->{classification};
		$spgeneinfo .= "\t$record->{pmid}\t$record->{assertion}";
	}

	return $spgeneinfo;

}


sub getUniprotkbAccns {

	my ($self, $id_type, $id) = @_;
	my @accns = ();

	my $core = "/id_ref";
	my $query = "/select?q=id_type:$id_type+AND+id_value:$id";
	my $fields = "&fl=uniprotkb_accession";
	my $format = "&wt=json&rows=10000";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;

	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		push @accns, $record->{uniprotkb_accession};
	}

	return @accns;

}


sub getIDs {

	my ($self, @accns) = @_;
	my @ids = ();

	my $core = "/id_ref";
	my $query = "/select?q=uniprotkb_accession:(". join(" OR ", @accns). ")";
	my $fields = "&fl=id_type,id_value";
	my $format = "&wt=json&rows=10000";

	my $solrQ = $solr.$core.$query.$fields.$format; 
	my $result = `wget -q -O - "$solrQ"`;

	my $resultObj = decode_json($result);

	foreach my $record (@{$resultObj->{'response'}->{'docs'}}){
		my $id_str = $record->{id_type}.'|'.$record->{id_value};
		push @ids, $id_str
			unless ($record->{id_type}=~/GI|GeneID|EnsemblGenome|EnsemblGenome_PRO|EnsemblGenome_TRS|PATRIC|EMBL|EMBL-CDS|KEGG|BioCyc|NCBI_TaxID|RefSeq_NT/ || (grep {$_ eq $id_str} @ids) );
	}

	return @ids;

}
