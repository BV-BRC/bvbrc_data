# post.jar is included in solr package
POSTJAR=/homes/mshukla/patric/git/patric3_data_infrastructure/data_processing/src/post.jar
URL=$PATRIC_SOLR/$1/update

java -Durl=$URL -Dauto=yes -Dfiletypes=json -Dcommit=yes -Dout=yes -Drecursive=yes -jar $POSTJAR $2
