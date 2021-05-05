cd `dirname $0`

cp ../network/full/weighted_nodes.tsv nodes.txt
cp ../network/full/raw_edges.tsv edgesr.txt
cp ../network/full/weighted_edges.tsv edgesw.txt

R --slave --args nodes.txt edgesw.txt resultsw.txt < random_walk.r
R --slave --args nodes.txt edgesr.txt resultsr.txt < random_walk.r

