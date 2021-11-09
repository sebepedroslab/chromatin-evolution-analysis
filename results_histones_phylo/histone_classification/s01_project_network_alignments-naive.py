# import libraries
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf
import os
import networkx as nx
import sys
# os.chdir("/home/xavi/Documents/histonome-analysis/results_histones_phylo/histone_classification/")

if len(sys.argv) < 4:
	print("Two arguments needed:")
	print("1. diamond self-alignments")
	print("2. alignments to a reference (HistDB, human...)")
	print("3. folder for FASTA outputs")
	exit()

# input files
pal_fn = sys.argv[1]
cla_fn = sys.argv[2]
out_fo = sys.argv[3]
out_fn = cla_fn.replace(".csv.gz", "")

if not os.path.exists("split_networks/"):
	os.makedirs("split_networks")


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
pal = pal[pal["evalue"] < 1e-10]
pal = pal[pal["bitscore"] > 20]
pal = pal[pal["qseqid"] != pal["sseqid"]]
print("# Pruned network size = %i alignments" % len(pal))

# subset classification to top quality
print("# Classification = %i alignments" % len(cla))
# cla = cla[cla["evalue"] < 1e-10]
cla = cla[cla["bitscore"] > 20]
print("# Pruned classification = %i alignments" % len(cla))

# create network
pal_n = nx.convert_matrix.from_pandas_edgelist(pal, source="qseqid", target="sseqid", edge_attr="bitscore")
pal_d = pd.DataFrame()

#### Connected components ####
pdfc = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-cla.pdf" % (out_fn))
pdfs = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-sps.pdf" % (out_fn))

# ubset network to the largest components, and plot them separately
for n,members in enumerate(nx.connected_components(pal_n)):

	# store members of connected component

	pal_d_i = pd.DataFrame()
	pal_d_i["members"] = list(members)
	pal_d_i["component"] = n

	if len(members) >= 10:

		print("# Network component %i (n=%i)" % (n, len(members)))
		pal_ns = pal_n.subgraph(members)

		pal_pos = nx.spring_layout(pal_ns, weight="bitscore", iterations=100, seed=111)
		# get list of edge weights
		pal_weight = list(nx.get_edge_attributes(pal_ns,'bitscore').values())

		# pal_com = list(nx.algorithms.community.asyn_lpa_communities(pal_ns, weight="bitscore"))
		# pal_ns_edges = np.array(list(pal_ns.edges()))

		# pal_ns_edges[0][1]

		# for com in pal_com:
		# 	np.sum(np.isin(pal_ns_edges[:,0], members))
		# 	# np.sum(np.logical_or(np.isin(pal_ns_edges[:,0], com) , np.isin(pal_ns_edges[:,1], com)))


		# pal_coms = [ pal_com[com] for com in pal_com ]

		# vector of colors per histone class
		node_class = [ cla[cla["qseqid"] == node]["classification"].values.tolist() for node in pal_ns.nodes() ]
		# node_class = [ "" if cla == [] else np.unique(cla)[0] for cla in node_class ] # retain most common alignment
		node_class = [ "" if cla == [] else cla[0] for cla in node_class ] # retain FIRST alignment (best)
		clas_color = [ color_dict_cla[i] if i in set(color_dict_cla.keys()) else "gray" for i in node_class ]

		# plot network, colored by histone class
		plt.figure(figsize=(10,10))
		plt.title("Network component %i (n=%i)" % (n, len(members)))
		# nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
		nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=clas_color, node_size=25, edgecolors=None, linewidths=0)
		# #nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
		plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_cla.values() ] , labels = [ col for col in color_dict_cla.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
		plt.savefig("split_networks/%s.Spring-Naive-cla.%i.n.svg" % (out_fn,n))
		nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
		pdfc.savefig()
		plt.savefig("split_networks/%s.Spring-Naive-cla.%i.ne.svg" % (out_fn,n))
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
		plt.figure(figsize=(10,10))
		plt.title("Network component %i (n=%i)" % (n, len(members)))
		#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
		nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=spss_color, node_size=15, edgecolors=None, linewidths=0)
		nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", width=1)
		#nx.draw_networkx_labels(G=pal_ns_brief, pos=pal_pos_brief, font_size=6, alpha=0.5)
		plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_sps.values() ] , labels = [ col for col in color_dict_sps.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
		pdfs.savefig()
		plt.close()

		# classification table
		pal_d_i["classification"] = node_class
		pal_d_i["posx"] = [ pal_pos[node][0] for node in pal_ns.nodes() ]
		pal_d_i["posy"] = [ pal_pos[node][1] for node in pal_ns.nodes() ]

		# list of domains in this component
		if not os.path.exists("%s/%s/" % (out_fo,out_fn)):
			os.makedirs("%s/%s/" % (out_fo,out_fn))
		node_class_most_common_list = np.unique(node_class, return_counts=True)
		node_class_most_common = node_class_most_common_list[0][0]
		pal_d_i.to_csv("%s/%s/CC%i_%s.txt" % (out_fo,out_fn,n,node_class_most_common), sep="\t", index=False, columns=["members"], header=False)

	pal_d = pd.concat([pal_d, pal_d_i])

pdfc.close()
pdfs.close()

# output classification table
pal_d.to_csv("%s.Spring-Naive.csv" % (out_fn), sep="\t", index=False)


# #### All components ####
# print("# Full network (n=%i)" % (len(pal_n)))

# pdfc = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-cla-all.pdf" % (out_fn))
# pdfs = matplotlib.backends.backend_pdf.PdfPages("%s.Spring-Naive-sps-all.pdf" % (out_fn))

# pal_ns = pal_n

# pal_pos = nx.spring_layout(pal_ns, weight="bitscore", iterations=1000, seed=111)
# # get list of edge weights
# pal_weight = list(nx.get_edge_attributes(pal_n,'bitscore').values())

# # vector of colors per histone class
# node_class = [ cla[cla["qseqid"] == node]["classification"].values.tolist() for node in pal_ns.nodes() ]
# # node_class = [ "" if cla == [] else np.unique(cla)[0] for cla in node_class ] # retain most common alignment
# node_class = [ "" if cla == [] else cla[0] for cla in node_class ] # retain FIRST alignment (best)
# clas_color = [ color_dict_cla[i] if i in set(color_dict_cla.keys()) else "gray" for i in node_class ]

# # plot network, colored by histone class
# plt.figure(figsize=(20,20))
# plt.title("Full network (n=%i)" % (len(pal_n)))
# #nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
# nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
# nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=clas_color, node_size=50)
# #nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
# plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_cla.values() ] , labels = [ col for col in color_dict_cla.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
# pdfc.savefig()
# plt.close()

# # vector of colors per species
# node_sps =   [ node.split("_")[0] for node in pal_ns.nodes() ]
# spss_color = [ color_dict_sps[i] if i in set(sps_list) else "slategray" for i in node_sps ]

# # positions of nodes to color
# pal_pos_brief = dict()
# for pi in pal_ns.nodes():
# 	pisp = pi.split("_")[0]
# 	if pisp in set(sps_list):
# 		pal_pos_brief[pi] = pal_pos[pi]
# pal_ns_brief = pal_n.subgraph(set(pal_pos_brief.keys()))

# # plot network, colored by species
# plt.figure(figsize=(20,20))
# plt.title("Full network (n=%i)" % (len(pal_n)))
# #nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
# nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
# nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=spss_color, node_size=50)
# nx.draw_networkx_labels(G=pal_ns_brief, pos=pal_pos_brief, font_size=6, alpha=0.5)
# plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_sps.values() ] , labels = [ col for col in color_dict_sps.keys() ], bbox_to_anchor=(1,1), loc="upper left", fontsize=6)
# pdfs.savefig()
# plt.close()


# pdfc.close()
# pdfs.close()
