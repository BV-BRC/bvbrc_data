# post.jar is included in solr package
POSTJAR=/homes/mshukla/bin/post.jar
URL=$PATRIC_SOLR_DEV/$1/update

java -Durl=$URL -Dauto=yes -Dfiletypes=json -Dcommit=yes -Dout=yes -Drecursive=yes -jar $POSTJAR $2
