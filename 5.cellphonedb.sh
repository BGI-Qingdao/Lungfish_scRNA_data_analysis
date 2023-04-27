## 1.data_pre
Rscript data_pre.r lung.celltype.rds lungfish_eggnog_one2one.trans.txt AE_lung

## 2.run cellphonedb
cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt --counts-data=gene_name --threads 8
cellphonedb plot dot_plot --means-path out/means.txt --pvalues-path out/pvalues.txt --output-path ./
cellphonedb plot heatmap_plot cellphonedb_meta.txt --pvalues-path out/pvalues.txt --output-path ./

## 3.visualization
Rscript bubble_plot.r out/pvalues.txt out/means.txt all all out
Rscript net_plot.r count_network.txt 1 1 0.3 AE_lung
