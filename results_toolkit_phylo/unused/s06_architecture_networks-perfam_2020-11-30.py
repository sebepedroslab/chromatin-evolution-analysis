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
import argparse

# logging
logging.basicConfig(
	level=logging.INFO,
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s"
	)

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-i", "--in", required=True, help="Table with genes and orthogroups, as in gene1 <tab> OG1 (one gene per line).", type=str)
arp.add_argument("-a", "--arch", required=True, help="Table with domain architectures per gene, as in gene1 <tab> dom1 dom2... (one gene per line).", type=str)
arp.add_argument("-o", "--out", required=False, default=None, help="Output folder. Defaults to folder where input file (-i/--in) is found.", type=str)
arp.add_argument("-skipprint", "--skipprint", required=False, action="store_true", help="Print network?")
arp.add_argument("-removecentral", "--removecentral", required=False, action="store_true", help="Remove central node and map node size to # of connections to this node? (works nicely with highly centralised networks)")
arl = vars(arp.parse_args())


# variables
ogs_fn = arl["in"]

if arl["out"] is not None:
	out_fn = arl["out"] + "/" +  os.path.basename(ogs_fn)
else:
	out_fn = ogs_fn

arq_fn = arl["arch"]
skip_print = arl["skipprint"]
do_remove_central = arl["removecentral"]

#################
### FUNCTIONS ###
#################


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
			# important: markersize needs to be divided by two, because node drawing function relies on pyplot's scatter function, which
			# multiplies node size by two for some reason.
			legend_titles.append(width)

	return legend_objects, legend_titles


def natural_sort(l): 
	import re
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)


########################
######### MAIN #########
########################


# create legend objects (network edge widths)
# legend_edges, legend_edges_width = legend_edge_widths(legend_widths=[1,5,10,20,50], color="gray", alpha=0.3, factor=width_factor_edge)
# legend_nodes, legend_nodes_width = legend_edge_widths(legend_widths=[1,5,10,20,50], color="purple", alpha=0.3, factor=width_factor_node)

legend_objects, legend_titles = legend_widths(edge_widths=[1,5,25,10,20,50,100], node_widths=[1,5,25,10,20,50,100], alpha=0.3)


# load per-gene architectures
arq = pd.read_csv("%s" % (arq_fn), sep="\t", names=["gene","architecture"])
# remove genes with empty architectures
arq = arq.dropna()

# retrieve orthogroup, and domains in orthogroup
ogg = pd.read_csv("%s" % (ogs_fn), sep="\t")
logging.info("Domain duets | %s" % (ogs_fn))

# find archs
archs_in_ogi = arq [ np.isin(arq["gene"], ogg["gene"].values) ] ["architecture"].values

# obtain all domain duets
ogg_duets = np.empty(shape=(0,2) , dtype=str)
for m,arch in enumerate(archs_in_ogi): 
	ogi_duets = duets_from_string(string=arch)
	ogg_duets = np.concatenate((ogg_duets, ogi_duets))

# obtain unique domain duets
ogg_duets_u = np.unique(ogg_duets, return_counts=True, axis=0)
ogg_duets_d = pd.DataFrame(data = { 
	"i": ogg_duets_u[0][:,0], 
	"j" : ogg_duets_u[0][:,1], 
	"weight" : ogg_duets_u[1]
	})

# create network
ogg_duets_n = nx.convert_matrix.from_pandas_edgelist(ogg_duets_d, source="i", target="j", edge_attr="weight")

if do_remove_central:

	# find central node and remove it:
	duet_degr = ogg_duets_n.degree()
	duet_degr = dict(duet_degr)
	node_cent = max(duet_degr, key=duet_degr.get)

	# find number of connections to central node (which will be used as node sizes)
	node_size = [ ogg_duets_n.get_edge_data(node_cent, i)["weight"] if ogg_duets_n.get_edge_data(node_cent, i) is not None else 0 for i in ogg_duets_n if i != node_cent ]
	# remove central node
	logging.info("Domain duets | %s | Node size according to central node %s, degree = %i / %i" % (ogs_fn,node_cent, duet_degr[node_cent], len(ogg_duets_n)))
	logging.info("Domain duets | %s | Drop central node %s" % (ogs_fn, node_cent))
	ogg_duets_n.remove_node(node_cent)

	# add bounds to scale:
	node_size_scaled = node_size
	node_size_scaled = [ 5 if i <= 5 and i != 0 else i for i in node_size_scaled ]
	node_size_scaled = [ 200 if i >= 200 else i for i in node_size_scaled ]

	pd.DataFrame(data = { 
		"domain": ogg_duets_n.nodes(),
		"size" : node_size,
		"size_scaled" : node_size_scaled
		}).to_csv("%s.network_nodes.csv" % (out_fn), sep="\t", index=None, mode="w")

else:

	node_size = 20
	node_size_scaled = node_size

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
np.unique(duet_coml_list, return_index=True)
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

	pdf = matplotlib.backends.backend_pdf.PdfPages("%s.network.pdf" % (out_fn))

	# plot network
	plt.figure(figsize=(8,8))
	_=plt.title("%s\nn = %i" % (ogs_fn, archs_in_ogi.shape[0]))
	nx.draw_networkx_edges(G=ogg_duets_n, pos=ogg_duets_pos, edge_color="gray", alpha=0.3, width=np.sqrt(ogg_duets_weight))
	nx.draw_networkx_nodes(G=ogg_duets_n, pos=ogg_duets_pos, node_color=duet_coml_list2, node_size=node_size_scaled, cmap=plt.cm.hsv, linewidths=0)
	nx.draw_networkx_labels(G=ogg_duets_n, pos=ogg_duets_pos, font_size=5, alpha=0.3)
	# nx.draw_networkx_edge_labels(G=ogg_duets_n, pos=ogg_duets_pos, font_size=5, alpha=0.3)
	_=plt.legend(legend_objects, legend_titles, frameon=False, loc='upper left', bbox_to_anchor=(1.04,1))
	plt.tight_layout()
	pdf.savefig()
	plt.close()

	pdf.close()

# output networks
ogg_duets_d.to_csv("%s.network_duets.csv" % (out_fn), sep="\t", index=None, mode="w")

# output architectures of individual genes within networks
oga = arq [ np.isin(arq["gene"], ogg["gene"].values) ]
oga = ogg.merge(arq, how="left", left_on="gene", right_on="gene")
oga.to_csv("%s.network_genes.csv" % (out_fn), sep="\t", index=None, mode="w")

