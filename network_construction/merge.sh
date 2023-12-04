# build network respectively
mkdir ./tmp
Rscript ../bin/bray-curtis.R data/train_x.csv 100 ./tmp
Rscript ../bin/coat_p.R data/train_x.csv 100 1e-12 ./tmp
Rscript ../bin/huge.R data/train_x.csv 100 ./tmp
Rscript ../bin/MI.R data/train_x.csv 100 ./tmp
# merge 
Rscript ../bin/Modify_adj.R ./tmp 0.05 ./ case
Rscript ../bin/Merge.R ./tmp 0.05 case
