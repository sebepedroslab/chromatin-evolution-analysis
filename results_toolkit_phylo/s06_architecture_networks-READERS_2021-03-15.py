### Input ###
# libraries
import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends import backend_pdf
import string
import itertools
import logging

# logging
logging.basicConfig(
	level=logging.INFO,
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s"
	)

# files
arq_fn = "architectures.csv"
fam_fn = "../data/gene_families_hmm.csv"
ogs_fn = "orthogroups_euk.csv"
out_fn = "results_domains/"
skip_print = False
do_remove_central = True

### Functions ###
def duets_from_string(string, split_char=" "):

	array = np.unique(np.sort(string.split(split_char)))

	if array.shape[0] == 1:
		array = np.sort(np.append(array,array))

	duets = np.array(list(itertools.combinations(array, 2)))

	return duets

def legend_edge_widths(legend_widths, color="lightgray", alpha=1, factor=0.5):

	legend_lines = []
	legend_widths = sorted(legend_widths)
	for width in legend_widths:
		legend_lines.append(plt.Line2D([],[], linewidth=np.sqrt(width)*factor, color=color, alpha=alpha))

	return legend_lines, legend_widths


def legend_widths(edge_widths=None, node_widths=None, color_edge="gray", color_node="purple", alpha=1):

	legend_objects = []
	legend_titles  = []
	if edge_widths is not None:
		edge_widths = sorted(edge_widths)
		for width in edge_widths:
			legend_objects.append(plt.Line2D([],[], linewidth=np.sqrt(width), color=color_edge, alpha=alpha))
			legend_titles.append(width)
	if node_widths is not None:
		node_widths = sorted(node_widths)
		for width in node_widths:
			legend_objects.append(plt.Line2D([],[], linewidth=np.sqrt(width), color=color_node, alpha=alpha))
			# legend_objects.append(plt.Line2D([],[], markersize=np.sqrt(width),linewidth=np.sqrt(width), color=color_node, alpha=alpha))
			legend_titles.append(width)

	return legend_objects, legend_titles


def natural_sort(l): 
	import re
	convert = lambda text: int(text) if text.isdigit() else text.lower()
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
	return sorted(l, key = alphanum_key)



### Main ###

legend_objects, legend_titles = legend_widths(edge_widths=[1,5,25,10,20,50,100,200], node_widths=[1,5,25,10,20,50,100,200], alpha=0.3)

# load list of gene families:
fam = pd.read_csv(fam_fn, sep="\t", names=["class","type","family","domains","search_param","inflation","minsize"])

# load per-gene architectures
arq = pd.read_csv(arq_fn, sep="\t", names=["gene","architecture"])
arq = arq.dropna()

# load gene classification (per family)
ogs = pd.read_csv(ogs_fn, sep="\t")
ogs["family"] = [i[0] for i in np.char.split(ogs["orthogroup"].values.astype(str), sep=".")]

# retrieve orthogroup, and domains in orthogroup
# ogg = pd.read_csv("%s" % (ogs_fn), sep="\t")
families_in_readers = fam[ fam["type"] == "Readers" ] [ "family" ].values
genes_in_family = ogs[ np.isin( ogs["family"], families_in_readers ) ] [ "gene" ].values
genes_in_family = np.unique(genes_in_family)

if len(genes_in_family) > 0:

	# find archs
	archs_in_ogi = arq [ np.isin(arq["gene"], genes_in_family) ] ["architecture"].values

	# obtain all domain duets
	ogg_duets = np.empty(shape=(0,2) , dtype=str)
	for m,arch in enumerate(archs_in_ogi): 
		ogi_duets = duets_from_string(string=arch)
		ogg_duets = np.concatenate((ogg_duets, ogi_duets))

	# obtain unique domain duets
	ogg_duets_u = np.unique(ogg_duets, return_counts=True, axis=0)

	# subset: remove rare duets
	ix_rare_duet = np.where(ogg_duets_u[1] > 5)
	ogg_duets_u_counts = ogg_duets_u[1][ix_rare_duet]
	ogg_duets_u_ijs = ogg_duets_u[0][ix_rare_duet]
	ogg_duets_d = pd.DataFrame(data = { 
		"i": ogg_duets_u_ijs[:,0], 
		"j" : ogg_duets_u_ijs[:,1], 
		"weight" : ogg_duets_u_counts
		})

	# create network
	ogg_duets_n = nx.convert_matrix.from_pandas_edgelist(ogg_duets_d, source="i", target="j", edge_attr="weight")

	# log
	logging.info("Domain duets | %s | %i genes | %i nodes" % ("readers", len(genes_in_family), len(ogg_duets_n)))

	# # drop reader domains
	# for core_domain in families_in_readers:
	# 	# ix_to_drop = np.where(np.isin(np.array(ogg_duets_n.nodes()), core_domain))[0][0]
	# 	# node_size_scaled.pop(ix_to_drop)
	# 	ogg_duets_n.remove_node(core_domain)

	# find node degrees
	duet_degr = ogg_duets_n.degree()
	duet_degr = dict(duet_degr)

	# node sizes
	node_size = [ duet_degr[i] for i in ogg_duets_n.nodes() ]
	print(node_size, len(node_size))
	# add bounds to scale:
	node_size_scaled = node_size
	node_size_scaled = [ 5 if i <= 5 and i < 0.5 else i for i in node_size_scaled ]
	node_size_scaled = [ 200 if i >= 200 else i for i in node_size_scaled ]
	node_cent = "None (included in plot)"

	# for core_domain in families_in_readers:
	# 	ix_to_change = np.where(np.isin(np.array(ogg_duets_n.nodes()), core_domain))[0][0]
	# 	node_size_scaled [ ix_to_change ] = 50

	# now find communities
	duet_comm = nx.algorithms.community.label_propagation.label_propagation_communities(G=ogg_duets_n)
	#duet_comm = nx.algorithms.community.modularity_max.greedy_modularity_communities(G=ogg_duets_n)
	duet_coml = list(duet_comm)
	duet_coml = sorted(duet_coml, key=len, reverse=True)

	# map communities to nodes, in network node order
	duet_coml_dict = dict()
	for n,i in enumerate(duet_coml):
		for j in i:
			duet_coml_dict[j] = n
	duet_coml_list = [ duet_coml_dict[i] + 1 for i in ogg_duets_n ]
	
	# which elements are singletons? color them a single color
	singleton_count = np.unique(duet_coml_list, return_counts=True)
	singleton_commu = singleton_count[0][singleton_count[1] == 1]
	duet_coml_list2 = [ 0 if i in set(singleton_commu) else i for i in duet_coml_list  ]


	# create layout
	# ogg_duets_pos = nx.spring_layout(ogg_duets_n, weight=None, iterations=30)
	ogg_duets_pos = nx.nx_agraph.graphviz_layout(ogg_duets_n, prog="neato")
	# get list of edge weights
	ogg_duets_weight = list(nx.get_edge_attributes(ogg_duets_n,'weight').values())

	if not skip_print:

		pdf = matplotlib.backends.backend_pdf.PdfPages("%s/reader_network.%s.network.pdf" % (out_fn,"readers"))

		# plot network
		plt.figure(figsize=(3.6,3.6))
		_=plt.title("%s\nn = %i\nCentral domain is %s" % ("readers", archs_in_ogi.shape[0], node_cent))
		nx.draw_networkx_edges(G=ogg_duets_n, pos=ogg_duets_pos, edge_color="darkgray", alpha=0.3, width=1)
		nx.draw_networkx_nodes(G=ogg_duets_n, pos=ogg_duets_pos, node_color="gray", node_size=node_size_scaled, linewidths=0.1)
		nx.draw_networkx_labels(G=ogg_duets_n, pos=ogg_duets_pos, font_size=5, alpha=0.3, font_color="olive")
		# nx.draw_networkx_edge_labels(G=ogg_duets_n, pos=ogg_duets_pos, font_size=5, alpha=0.3)
		_=plt.legend(legend_objects, legend_titles, frameon=False, loc='upper left', bbox_to_anchor=(1.04,1))
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# plot network
		plt.figure(figsize=(9,9))
		_=plt.title("%s\nn = %i\nCentral domain is %s" % ("readers", archs_in_ogi.shape[0], node_cent))
		nx.draw_networkx_edges(G=ogg_duets_n, pos=ogg_duets_pos, edge_color="darkgray", alpha=0.3, width=1)
		nx.draw_networkx_nodes(G=ogg_duets_n, pos=ogg_duets_pos, node_color=duet_coml_list2, node_size=node_size_scaled, linewidths=0.1)
		nx.draw_networkx_labels(G=ogg_duets_n, pos=ogg_duets_pos, font_size=5, alpha=0.3, font_color="olive")
		# nx.draw_networkx_edge_labels(G=ogg_duets_n, pos=ogg_duets_pos, font_size=5, alpha=0.3)
		_=plt.legend(legend_objects, legend_titles, frameon=False, loc='upper left', bbox_to_anchor=(1.04,1))
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# plot degrees
		duet_degr_sorkey = sorted(duet_degr, key=duet_degr.get, reverse=True)
		duet_degr_sorval = [ duet_degr[i] for i in duet_degr_sorkey ]

		plt.figure(figsize=(6,6))
		top_degree = min(20, len(duet_degr_sorval))
		_=plt.title("%s node degree (top %i)" % ("readers",top_degree))
		plt.bar(range(top_degree), duet_degr_sorval[0:top_degree], align='center', color="blue")
		plt.xticks(range(top_degree), duet_degr_sorkey[0:top_degree], size=6, rotation='vertical')
		pdf.savefig()
		plt.close()

		pdf.close()

	# output networks
	ogg_duets_d.to_csv("%s/reader_network.%s.network_duets.csv" % (out_fn,"readers"), sep="\t", index=None, mode="w")

	# output architectures of individual genes within networks
	oga = ogs [ np.isin(ogs["gene"], genes_in_family) ]
	oga = oga.merge(arq, how="left", left_on="gene", right_on="gene")
	oga.to_csv("%s/reader_network.%s.network_genes.csv" % (out_fn,"readers"), sep="\t", index=None, mode="w")

else:

	logging.info("Domain duets | %s | No genes in this domain! Nothing to do." % ("readers general"))


