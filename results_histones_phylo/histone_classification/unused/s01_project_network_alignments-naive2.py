# import libraries
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf
import os
import networkx as nx
import sys
# os.chdir("/home/xavi/Documents/Lab/histonome-analysis/results_histones_phylo/histone_classification/")

if len(sys.argv) < 3:
	print("Two arguments needed:")
	print("1. diamond self-alignments")
	print("2. alignments to a reference (HistDB, human...)")
	exit()

# input files
pal_fn = sys.argv[1]
cla_fn = sys.argv[2]
out_fn = cla_fn.replace(".csv.gz", "")

# load
pal = pd.read_csv(pal_fn, sep="\t", names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"], compression="gzip")
cla = pd.read_csv(cla_fn, sep="\t", names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"], compression="gzip")
cla["classification"] = cla["sseqid"].str.replace(".*\\|","")

# color dictionary
# list of histone types (includes main non-canonical histones)
his_list = np.unique(cla["classification"])
cmap = plt.get_cmap('viridis')
colors = cmap( np.linspace(0, 1, len(np.unique(his_list)) + 0 ) )
color_dict_cla = dict()
for n,his in enumerate(his_list):
	color_dict_cla[his] = colors[n]

# species dictionary for selected species for which we have PTMs
sps_list = ["Hsap","Dmel","Scil","Cowc","Cfra","Clim","Scer","Spom","Ncra","Spun","Ttra","Acas","Ddis","Atha","Ppat","Crei","Gthe","Ehux","Tpse","Phatri","Esil","Pinf","Tthe","Bnat","Ngru","arc","vir"]
# sps_list = ["arc","vir"]
cmap = plt.get_cmap('plasma')
colors = cmap( np.linspace(0, 1, len(np.unique(sps_list)) + 0 ) )
color_dict_sps = dict()
for n,sps in enumerate(sps_list):
	color_dict_sps[sps] = colors[n]

# subset alignments to top quality
print("# Network size = %i alignments" % len(pal))
pal = pal[pal["evalue"] < 1e-20]
pal = pal[pal["bitscore"] > 20]
print("# Pruned network size = %i alignments" % len(pal))

# subset classification to top quality
print("# Classification = %i alignments" % len(cla))
# cla = cla[cla["evalue"] < 1e-20]
cla = cla[cla["bitscore"] > 20]
print("# Pruned classification = %i alignments" % len(cla))

# create network
pal_n = nx.convert_matrix.from_pandas_edgelist(pal, source="qseqid", target="sseqid", edge_attr="bitscore")


#### Connected components ####
pdfc = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-cla.pdf" % (out_fn))
pdfs = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-sps.pdf" % (out_fn))

# identify connected components
pal_cc = nx.connected_components(pal_n)
print(type(pal_cc))
pal_cd = { c for c in pal_cc }
print(pal_cd)

n,members in enumerate(pal_cc)

# # ubset network to the largest components, and plot them separately
# for n,members in enumerate(pal_cc):

# 	if len(members) >= 10:

# 		print("# Network component %i (n=%i)" % (n, len(members)))

# 		pal_ns = pal_n.subgraph(members)

# 		pal_pos = nx.spring_layout(pal_ns, weight="bitscore", iterations=100, seed=111)
# 		# get list of edge weights
# 		pal_weight = list(nx.get_edge_attributes(pal_n,'bitscore').values())

# 		# vector of colors per histone class
# 		node_class = [ cla[cla["qseqid"] == node]["classification"].values.tolist() for node in pal_ns.nodes() ]
# 		# node_class = [ "" if cla == [] else np.unique(cla)[0] for cla in node_class ] # retain most common alignment
# 		node_class = [ "" if cla == [] else cla[0] for cla in node_class ] # retain FIRST alignment (best)
# 		clas_color = [ color_dict_cla[i] if i in set(color_dict_cla.keys()) else "gray" for i in node_class ]

# 		# plot network, colored by histone class
# 		plt.figure(figsize=(10,10))
# 		plt.title("Network component %i (n=%i)" % (n, len(members)))
# 		#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
# 		nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
# 		nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=clas_color, node_size=50)
# 		#nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
# 		plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_cla.values() ] , labels = [ col for col in color_dict_cla.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
# 		pdfc.savefig()
# 		plt.close()

# 		# vector of colors per species
# 		node_sps =   [ node.split("_")[0] for node in pal_ns.nodes() ]
# 		spss_color = [ color_dict_sps[i] if i in set(sps_list) else "slategray" for i in node_sps ]

# 		# positions of nodes to color
# 		pal_pos_brief = dict()
# 		for pi in pal_ns.nodes():
# 			pisp = pi.split("_")[0]
# 			if pisp in set(sps_list):
# 				pal_pos_brief[pi] = pal_pos[pi]
# 		pal_ns_brief = pal_n.subgraph(set(pal_pos_brief.keys()))

# 		# plot network, colored by species
# 		plt.figure(figsize=(10,10))
# 		plt.title("Network component %i (n=%i)" % (n, len(members)))
# 		#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
# 		nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
# 		nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=spss_color, node_size=50)
# 		nx.draw_networkx_labels(G=pal_ns_brief, pos=pal_pos_brief, font_size=6, alpha=0.5)
# 		plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_sps.values() ] , labels = [ col for col in color_dict_sps.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
# 		pdfs.savefig()
# 		plt.close()

# 		# TODO
# 		# report table with classification!

# pdfc.close()
# pdfs.close()



#### All components ####
print("# Full network (n=%i)" % (len(pal_n)))

pdfc = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-cla-allp.pdf" % (out_fn))
pdfs = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-sps-allp.pdf" % (out_fn))

pal_ns = pal_n

pal_pos = nx.spring_layout(pal_ns, pos=pal_cc, weight="bitscore", iterations=10000, seed=111)
# get list of edge weights
pal_weight = list(nx.get_edge_attributes(pal_n,'bitscore').values())

# vector of colors per histone class
node_class = [ cla[cla["qseqid"] == node]["classification"].values.tolist() for node in pal_ns.nodes() ]
# node_class = [ "" if cla == [] else np.unique(cla)[0] for cla in node_class ] # retain most common alignment
node_class = [ "" if cla == [] else cla[0] for cla in node_class ] # retain FIRST alignment (best)
clas_color = [ color_dict_cla[i] if i in set(color_dict_cla.keys()) else "gray" for i in node_class ]

# plot network, colored by histone class
plt.figure(figsize=(20,20))
plt.title("Full network (n=%i)" % (len(pal_n)))
#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=clas_color, node_size=50)
#nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_cla.values() ] , labels = [ col for col in color_dict_cla.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
pdfc.savefig()
plt.close()

# vector of colors per species
node_sps =   [ node.split("_")[0] for node in pal_ns.nodes() ]
spss_color = [ color_dict_sps[i] if i in set(sps_list) else "slategray" for i in node_sps ]

# positions of nodes to color
pal_pos_brief = dict()
for pi in pal_ns.nodes():
	pisp = pi.split("_")[0]
	if pisp in set(sps_list):
		pal_pos_brief[pi] = pal_pos[pi]
pal_ns_brief = pal_n.subgraph(set(pal_pos_brief.keys()))

# plot network, colored by species
plt.figure(figsize=(20,20))
plt.title("Full network (n=%i)" % (len(pal_n)))
#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=spss_color, node_size=50)
nx.draw_networkx_labels(G=pal_ns_brief, pos=pal_pos_brief, font_size=6, alpha=0.5)
plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_sps.values() ] , labels = [ col for col in color_dict_sps.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
pdfs.savefig()
plt.close()


pdfc.close()
pdfs.close()
