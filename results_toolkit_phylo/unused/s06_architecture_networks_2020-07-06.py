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
arp.add_argument("-cluname", "--cluname", required=False, default="cluster_name", help="Cluster column name (defaults to \"cluster_name\").", type=str)
arp.add_argument("-skipprint", "--skipprint", required=False, default=False, help="Print network?", type=bool)
arl = vars(arp.parse_args())


# variables
ogs_fn = arl["in"]

if arl["out"] is not None:
	out_fn = arl["out"] + "/" +  os.path.basename(ogs_fn)
else:
	out_fn = ogs_fn

cluname = arl["cluname"]
arq_fn = arl["arch"]
skip_print = arl["skipprint"]



#################
### FUNCTIONS ###
#################


# def duets_from_string(string, constant_element=None, split_char=" "):

# 	array = np.unique(np.sort(string.split(split_char)))

# 	if constant_element is not None:
# 		array = np.sort(np.append(array,constant_element))

# 	duets = np.array(list(itertools.combinations(array, 2)))
# 	# duets = pd.DataFrame( data = { "i" : duets[:,0], "j" : duets[:,1] } )

# 	return duets


def duets_from_string(string, split_char=" "):

	array = np.unique(np.sort(string.split(split_char)))

	if array.shape[0] == 1:
		array = np.sort(np.append(array,array))

	duets = np.array(list(itertools.combinations(array, 2)))

	return duets

def domain_counts_from_string(string, split_char=" "):

	array = np.unique(np.sort(string.split(split_char)))

	if array.shape[0] == 1:
		array = np.sort(np.append(array,array))

	duets = np.array(list(itertools.combinations(array, 2)))

	return duets


def legend_edge_widths(legend_widths, color="lightgray", alpha=1):

	legend_lines = []
	legend_widths = sorted(legend_widths)
	for width in legend_widths:
		legend_lines.append(plt.Line2D([],[], linewidth=np.sqrt(width), color=color, alpha=alpha))

	return legend_lines, legend_widths


def natural_sort(l): 
	import re
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)


########################
######### MAIN #########
########################


# create legend objects (network edge widths)
legend_lines, legend_widths = legend_edge_widths(legend_widths=[1,5,10,20,50,100], color="gray", alpha=0.3)


# load per-gene architectures
arq = pd.read_csv("%s" % (arq_fn), sep="\t", names=["gene","architecture"])
# remove genes with empty architectures
arq = arq.dropna()

# load orthologous groups
ogs = pd.read_csv("%s" % (ogs_fn), sep="\t")
ogs_list = natural_sort(np.unique(ogs[cluname]))
ogs_list.append("ALL")

# open pdf
if not skip_print:
	pdf = matplotlib.backends.backend_pdf.PdfPages("%s.network.pdf" % (out_fn))

# loop through orthogroup
archs_in_ogs_counts = pd.DataFrame()
ogg_duets_concat = pd.DataFrame()
for ogi in ogs_list:
	
	# retrieve orthogroup, and domains in orthogroup
	logging.info("Domain duets | %s | %s" % (ogs_fn,ogi))

	# select nodes to analyse
	if ogi is not "ALL":
		ogg = ogs[ogs[cluname] == ogi]
	else:
		ogg = ogs
	
	# find archs
	archs_in_ogi = arq [ np.isin(arq["gene"], ogg["gene"].values) ] ["architecture"].values

	# counts of common architectures (ignore singletons)
	archs_in_ogi_unique = np.unique(archs_in_ogi, return_counts=True)
	archs_in_ogi_counts = pd.DataFrame(data = { 
		"group" : ogi, 
		"architecture": archs_in_ogi_unique[0] [ archs_in_ogi_unique[1] > 0 ],
		"frequency":    archs_in_ogi_unique[1] [ archs_in_ogi_unique[1] > 0 ],
		})
	# concatenate
	archs_in_ogi_counts = archs_in_ogi_counts.sort_values(by="frequency", ascending=False)
	archs_in_ogs_counts = pd.concat((archs_in_ogs_counts, archs_in_ogi_counts))

	# obtain all domain duets
	ogg_duets = np.empty(shape=(0,2) , dtype=str)
	for m,arch in enumerate(archs_in_ogi): 
		ogi_duets = duets_from_string(string=arch)
		ogg_duets = np.concatenate((ogg_duets, ogi_duets))

	if len(ogg_duets) > 0 and not skip_print:

		# obtain unique domain duets
		ogg_duets_u = np.unique(ogg_duets, return_counts=True, axis=0)
		ogg_duets_d = pd.DataFrame(data = { 
			"i": ogg_duets_u[0][:,0], 
			"j" : ogg_duets_u[0][:,1], 
			"weight" : ogg_duets_u[1] 
			})
		# store them...
		if ogi is not "ALL":
			ogg_duets_concat = pd.concat((ogg_duets_concat, ogg_duets_d), sort=False)
			ogg_duets_concat["orthogroup"] = ogi

		# create network
		ogg_duets_n = nx.convert_matrix.from_pandas_edgelist(ogg_duets_d, source="i", target="j", edge_attr="weight")
		ogg_duets_pos = nx.spring_layout(ogg_duets_n)
		# get list of edge weights
		ogg_duets_weight = list(nx.get_edge_attributes(ogg_duets_n,'weight').values())

		# plot network
		_=plt.title("%s\nn = %i" % (ogi, archs_in_ogi.shape[0]))

		nx.draw_networkx_edges(G=ogg_duets_n, pos=ogg_duets_pos, edge_color="gray", alpha=0.3, width=np.sqrt(ogg_duets_weight))
		nx.draw_networkx_nodes(G=ogg_duets_n, pos=ogg_duets_pos, node_color="darkorange", node_size=200)
		nx.draw_networkx_labels(G=ogg_duets_n, pos=ogg_duets_pos)
		_=plt.legend(legend_lines, legend_widths, frameon=False, loc='upper left', bbox_to_anchor=(1.04,1))
		plt.tight_layout()
		pdf.savefig()
		plt.close()

# output networks
archs_in_ogs_counts.to_csv("%s.network.csv" % (out_fn), sep="\t", index=None, mode="w")
ogg_duets_concat.to_csv("%s.duets.csv" % (out_fn), sep="\t", index=None, mode="w")

# output architectures of individual genes within networks
oga = arq [ np.isin(arq["gene"], ogs["gene"].values) ]
oga = ogs.merge(arq, how="left", left_on="gene", right_on="gene")
oga.to_csv("%s.genes.csv" % (out_fn), sep="\t", index=None, mode="w")

if not skip_print:
	pdf.close()

