# import libraries
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
os.chdir("/home/xavi/Documents/Lab/histonome-analysis/results_histones_phylo/histone_classification/")

# input files
pal_fn = ["euk.Histone.domains.diamond_to_human.csv","euk.Histone.domains.diamond_to_histdb.csv","vir.Histone.domains.diamond_to_histdb.csv","arc.Histone.domains.diamond_to_histdb.csv"]


for pal_fi in pal_fn:

	out_fi = pal_fi.replace(".csv", "")
	cla_fi = pal_fi.replace(".csv", ".clas.csv")

	# load
	pal = pd.read_csv(pal_fi, sep="\t", names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	cla = pd.read_csv(cla_fi, sep="\t")
	cla_seq = cla["qseqid"].values
	cla_cla = cla["classification"].values
	cla_clu, cla_cli = np.unique(cla_cla, return_inverse=True)
	unknown_label = cla_clu.shape[0] -1
	cla_cli[ cla_cli == unknown_label ] = -1

	cmap = plt.get_cmap('gist_rainbow')
	colors = cmap( np.linspace(0, 1, len(np.unique(cla_cla)) + 0 ) )
	cla_col_d = dict()
	for n,i in enumerate(np.unique(cla_cla)):
		cla_col_d[i] = colors[n,]

	cla_col = [ cla_col_d[i] for i in cla_cla ]

	# bitscore wide table
	pal_b = pal[["qseqid","sseqid","bitscore"]]
	pal_bw = pd.pivot(pal_b, index="qseqid", columns="sseqid", values="bitscore")
	pal_bw = pal_bw.fillna(0)

	# order data according to classification
	pal_bw = pal_bw.reindex(cla_seq)

	# supervised learning UMAP
	import umap

	# train UMAP embedding with test data
	print("# UMAP")
	pal_mapper = umap.UMAP(n_neighbors=50, n_components=2, metric="hamming").fit(X=pal_bw)
	pal_embedding = pal_mapper.embedding_

	# save pdf
	plt.figure(figsize=(5,5))
	plt.scatter(pal_embedding[:,0], pal_embedding[:,1], s=2, color=cla_col)
	plt.title("UMAP (n var = %i)" % pal_bw.shape[1])
	plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in cla_col_d.values() ] , labels = [ col for col in cla_col_d.keys() ], bbox_to_anchor=(1,1), loc="upper left")
	plt.savefig("%s.umap_projection.pdf" % out_fi, bbox_inches='tight')
	plt.close()


