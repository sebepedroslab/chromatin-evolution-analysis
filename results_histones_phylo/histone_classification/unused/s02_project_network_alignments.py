# import libraries
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf
import os
import networkx as nx
os.chdir("/home/xavi/Documents/Lab/histonome-analysis/results_histones_phylo/histone_classification/")

# input files
pal_fn = "euk.Histone.domains.diamond_to_self.csv"
cla_fn = "euk.Histone.domains.diamond_to_histdb.clas.csv"
out_fn = pal_fn.replace(".csv", "")

# load
pal = pd.read_csv(pal_fn, sep="\t", names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
pal.shape
pal = pal[pal["evalue"] < 1e-5]
#pal = pal[pal["bitscore"] > 20]
pal.shape
cla = pd.read_csv(cla_fn, sep="\t")
cla["classification"] = cla["classification"].str.replace("bitscore.","")


# color dictionary
# list of histone types (includes main non-canonical histones)
his_list = np.unique(cla["classification"])
cmap = plt.get_cmap('viridis')
colors = cmap( np.linspace(0, 1, len(np.unique(his_list)) + 0 ) )
color_dict_cla = dict()
for n,his in enumerate(his_list):
	color_dict_cla[his] = colors[n]

# species dictionary for selected species for which we have PTMs
sps_list = ["Hsap","Dmel","Scil","Cowc","Cfra","Clim","Scer","Spom","Ncra","Spun","Ttra","Acas","Ddis","Atha","Ppat","Crei","Gthe","Ehux","Tpse","Phatri","Esil","Pinf","Tthe","Bnat","Ngru"]
cmap = plt.get_cmap('plasma')
colors = cmap( np.linspace(0, 1, len(np.unique(sps_list)) + 0 ) )
color_dict_sps = dict()
for n,sps in enumerate(sps_list):
	color_dict_sps[sps] = colors[n]

# define edge width legend
def legend_edge_widths(legend_widths, color="lightgray", alpha=1):

	legend_lines = []
	legend_widths = sorted(legend_widths)
	for width in legend_widths:
		legend_lines.append(plt.Line2D([],[], linewidth=np.sqrt(width), color=color, alpha=alpha))

	return legend_lines, legend_widths

legend_lines, legend_widths = legend_edge_widths(legend_widths=[1,5,10,20,50,100], color="gray", alpha=0.3)
_=plt.legend(legend_lines, legend_widths, frameon=False, loc='upper left', bbox_to_anchor=(1.04,1))



#### SPRING LAYOUT
pdf = matplotlib.backends.backend_pdf.PdfPages("%s.Spring.pdf" % (out_fn))
for his in ["H2A","H2B","H3","H4","H2AZ","macroH2A"]:

	# filter
	cla_i = cla[cla["classification"] == his]
	cla_i_list = cla_i["qseqid"].values
	pal_i = pal[ np.logical_or( np.isin(pal["qseqid"], cla_i_list) , np.isin(pal["sseqid"], cla_i_list) ) ]

	# create network
	pal_n = nx.convert_matrix.from_pandas_edgelist(pal_i, source="qseqid", target="sseqid", edge_attr="bitscore")
	
	# subset network to the largest components, and plot them separately
	for n,members in enumerate(nx.connected_components(pal_n)):
		

		if len(members) > 20:
	
			print("Network %s - component %i" % (his, n))

			pal_ns = pal_n.subgraph(members)

			pal_pos = nx.spring_layout(pal_ns, weight="bitscore")
			# get list of edge weights
			pal_weight = list(nx.get_edge_attributes(pal_n,'bitscore').values())

			# vector of colors per histone class
			node_class = [ cla[cla["qseqid"] == node]["classification"].values.tolist() for node in pal_ns.nodes() ]
			node_class = [ "unknown" if cla == [] else cla[0] for cla in node_class ]
			clas_color = [ color_dict_cla[cla] for cla in node_class ]

			# plot network, colored by histone class
			plt.figure(figsize=(8,8))
			plt.title("network based on %s, component %i" % (his,n))
			#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
			nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
			nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=clas_color, node_size=50)
			#nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
			plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_cla.values() ] , labels = [ col for col in color_dict_cla.keys() ], bbox_to_anchor=(1,1), loc="upper left")
			pdf.savefig()
			plt.close()

			# vector of colors per species
			node_sps =   [ node.split("_")[0] for node in pal_ns.nodes() ]
			spss_color = [ color_dict_sps[sps] if sps in set(sps_list) else "slategray" for sps in node_sps ]


			# plot network, colored by histone class
			plt.figure(figsize=(8,8))
			plt.title("network based on %s, component %i" % (his,n))
			#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
			nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
			nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=spss_color, node_size=50)
			#nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
			plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_sps.values() ] , labels = [ col for col in color_dict_sps.keys() ], bbox_to_anchor=(1,1), loc="upper left")
			pdf.savefig()
			plt.close()

pdf.close()


exit()

#### Kamada-Kawai LAYOUT
pdf = matplotlib.backends.backend_pdf.PdfPages("%s.KamadaKawai.pdf" % (out_fn))
for his in ["H2A","H2B","H3","H4","H2AZ","macroH2A"]:

	# filter
	cla_i = cla[cla["classification"] == his]
	cla_i_list = cla_i["qseqid"].values
	pal_i = pal[ np.logical_or( np.isin(pal["qseqid"], cla_i_list) , np.isin(pal["sseqid"], cla_i_list) ) ]

	# create network
	pal_n = nx.convert_matrix.from_pandas_edgelist(pal_i, source="qseqid", target="sseqid", edge_attr="bitscore")
	
	# subset network to the largest components, and plot them separately
	for n,members in enumerate(nx.connected_components(pal_n)):
		
		if len(members) > 20:
			pal_ns = pal_n.subgraph(members)

			pal_pos = nx.kamada_kawai_layout(pal_ns, weight="bitscore")
			# get list of edge weights
			pal_weight = list(nx.get_edge_attributes(pal_n,'bitscore').values())

			# vector of colors per histone class
			node_class = [ cla[cla["qseqid"] == node]["classification"].values.tolist() for node in pal_ns.nodes() ]
			node_class = [ "unknown" if cla == [] else cla[0] for cla in node_class ]
			clas_color = [ color_dict_cla[cla] for cla in node_class ]

			# plot network, colored by histone class
			plt.figure(figsize=(8,8))
			plt.title("network based on %s, component %i" % (his,n))
			#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
			nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
			nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=clas_color, node_size=50)
			#nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
			plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_cla.values() ] , labels = [ col for col in color_dict_cla.keys() ], bbox_to_anchor=(1,1), loc="upper left")
			pdf.savefig()
			plt.close()

			# vector of colors per species
			node_sps =   [ node.split("_")[0] for node in pal_ns.nodes() ]
			spss_color = [ color_dict_sps[sps] if sps in set(sps_list) else "slategray" for sps in node_sps ]


			# plot network, colored by histone class
			plt.figure(figsize=(8,8))
			plt.title("network based on %s, component %i" % (his,n))
			#nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="gray", alpha=0.3, width=np.sqrt(pal_weight)/np.max(np.sqrt(pal_weight))*10)
			nx.draw_networkx_edges(G=pal_ns, pos=pal_pos, edge_color="silver", alpha=0.5, width=1)
			nx.draw_networkx_nodes(G=pal_ns, pos=pal_pos, node_color=spss_color, node_size=50)
			#nx.draw_networkx_labels(G=pal_ns, pos=pal_pos)
			plt.legend( handles = [ plt.Line2D([0,0],[0,0],color=col, marker='o', linestyle='') for col in color_dict_sps.values() ] , labels = [ col for col in color_dict_sps.keys() ], bbox_to_anchor=(1,1), loc="upper left")
			pdf.savefig()
			plt.close()

pdf.close()

