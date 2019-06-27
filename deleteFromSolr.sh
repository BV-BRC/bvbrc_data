echo "$PATRIC_SOLR/$1/update?stream.body=<delete><query>$2:$3</query></delete>&commit=true"
curl "$PATRIC_SOLR/$1/update?stream.body=<delete><query>$2:$3</query></delete>&commit=true"
