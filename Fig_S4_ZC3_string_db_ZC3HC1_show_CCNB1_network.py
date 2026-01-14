import requests
import webbrowser
import wget
from time import sleep
import networkx as nx
import pandas as pd
from pandas import json_normalize
import matplotlib.pyplot as plt
import json
from sklearn import cluster
import numpy as np
from matplotlib import cm
from matplotlib import patches
from matplotlib.colors import to_rgba
import datetime
import seaborn as sns

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 80)
pd.set_option('display.width', 1200)
pd.set_option('max_colwidth', 80)

server = "https://string-db.org/api/"

ext = "image/network?identifiers=LRRK2%0dSYT10"

today = str(datetime.date.today())
today = today.replace('-', '')

# webbrowser.open(server + ext)
# wget.download(server + ext, out="String_DB_out")


def get_network_interactions() -> 'Getting the STRING network interactions':
    """
    The network API method also allows you to retrieve your STRING interaction network
    for one or multiple proteins in various text formats. It will tell you the combined score
    and all the channel specific scores for the set of proteins. You can also extend the network
    neighborhood by setting "add_nodes", which will add, to your network, new interaction partners
    in order of their confidence.

    Format 	Description
    tsv 	tab separated values, with a header line
    tsv-no-header 	tab separated values, without header line
    json 	JSON format
    xml 	XML format
    psi-mi 	PSI-MI XML format
    psi-mi-tab 	PSI-MITAB format

    Parameter 	Description
    identifiers 	required parameter for multiple items, e.g. DRD1_HUMAN%0dDRD2_HUMAN
    species 	NCBI taxon identifiers (e.g. Human is 9606, see: STRING organisms).
    required_score 	threshold of significance to include a interaction,
                    a number between 0 and 1000 (default depends on the network)
    network_type 	network type: functional (default), physical
    add_nodes 	adds a number of proteins with to the network based on their confidence score
    show_query_node_labels 	when available use submitted names in the preferredName column when (0 or 1) (default:0)
    caller_identity 	your identifier for us.

    score 	combined score
    nscore 	gene neighborhood score
    fscore 	gene fusion score
    pscore 	phylogenetic profile score
    ascore 	coexpression score
    escore 	experimental score
    dscore 	database score
    tscore 	textmining score
    """

    output_format = "json"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])
    my_gene_temp = my_genes
    my_gene_temp.append("CCNB1")

    params = {
        "identifiers":  "%0d".join(my_gene_temp),  # your proteins
        "species": 9606,  # species NCBI identifier = Homo sapiens
        "add_white_nodes": 20,  # add 15 white nodes to my protein (40)
        "network_type": "functional",  # network type: functional (default), physical
        "network_flavor": "confidence",  # the style of edges in the network: evidence, confidence (default), actions
        "caller_identity": "Tobias_Reinberger_Uni_Luebeck"  # your app name
    }

    response = requests.post(request_url, data=params)
    decoded = response.json()
    df = json_normalize(decoded)
    print(df)
    # df = df[df.score > 0.4]
    # df = df[df["escore"] > 0]
    df.to_csv("String_DB_out/%s_network_2.csv" % "_".join(my_genes[:5]))
    G = nx.from_pandas_edgelist(df, "preferredName_A", "preferredName_B",
                                [
                                 # "score",
                                 "escore",
                                 "ascore",
                                 "dscore",
                                 # "tscore"
                                 ])

    #  Get CCBN1 subnetwork
    CCNB1_network_temp = [n for n in G.neighbors("CCNB1")]
    CCNB1_network_temp2 = [n for n in G.neighbors("CCNB1")]
    CCNB1_network = CCNB1_network_temp.copy()
    print(CCNB1_network_temp)
    print(CCNB1_network_temp2)
    for gene in CCNB1_network_temp:
        if gene != "CCNB1":
            print(gene)
            extra_genes = [n for n in G.neighbors(gene)]
            # print(extra_genes)
            CCNB1_network.extend(extra_genes)
    CCNB1_network = sorted(np.unique(CCNB1_network))
    print(CCNB1_network)

    params = {
        "identifiers":  "%0d".join(list(CCNB1_network)),  # your proteins
        "species": 9606,  # species NCBI identifier = Homo sapiens
        "add_white_nodes": 8,  # add 15 white nodes to my protein (40)
        "network_type": "functional",  # network type: functional (default), physical
        "network_flavor": "confidence",  # the style of edges in the network: evidence, confidence (default), actions
        "caller_identity": "Tobias_Reinberger_Uni_Luebeck"  # your app name
    }
    response = requests.post(request_url, data=params)
    decoded = response.json()
    df = json_normalize(decoded)
    print(df)
    # df = df[df.score > 0.4]
    # df = df[df["escore"] > 0]
    df.to_csv("String_DB_out/%s_network_2.csv" % "_".join(my_genes[:5]))
    G = nx.from_pandas_edgelist(df, "preferredName_A", "preferredName_B")

    G1 = G
    to_be_removed = []
    for x in G1.nodes():
        is_network = False
        for n in G1.neighbors(x):
            if G1.degree(n) > 1:
                is_network = True
                break
        if not is_network:
            to_be_removed.append(x)
    # to_be_removed = [x for x in G1.neighbors(x) if G1.degree(x) <= 1]
    print(to_be_removed)
    for x in to_be_removed:
        G1.remove_node(x)
    G = G1
    CCNB1_network_temp = [n for n in G.neighbors("CCNB1")]
    # G = nx.ego_graph(G, "CCNB1")
    clusters_sorted, cluster_dict = get_cluster(G)
    print(clusters_sorted)
    print(cluster_dict)
    enriched_clusters = []
    for idx, clusters in enumerate(cluster_dict.values()):
        enrichment_out = get_enrichment(clusters)
        enrichment_out.sort_values("p_value")
        enrichment_out = enrichment_out.head(2)
        enrichment_out["Cluster"] = idx + 1
        enrichment_out = enrichment_out[["Cluster", "number_of_genes", "inputGenes", "p_value", "fdr", "description", "category"]]
        enrichment_out = enrichment_out.rename(columns={"number_of_genes": "#Genes",
                                                        "fdr": "FDR",
                                                        "description": "GO process"})
        enrichment_out["inputGenes"] = enrichment_out["inputGenes"].apply(lambda x: ", ".join(x))
        enriched_clusters.append(enrichment_out)

    enrichment_out_sum = pd.concat(enriched_clusters)

    cluster_colors_temp = ["gray", "green", "slateblue", "firebrick", "orange", "darkgray",
                           "darkcyan", "lightsteelblue", "peru", "firebrick", "lightgreen", "c", "orchid",
                           "royalblue"]

    # cluster_colors_temp = sns.color_palette("Set2")
    enrichment_out_sum["cluster_color"] = enrichment_out_sum["Cluster"].apply(lambda x: cluster_colors_temp[x - 1])
    print(enrichment_out_sum)
    enrichment_out_sum.to_excel("Enrichment/%s_%s.xlsx" % (today, "_".join(my_genes[:4])))

    edge_colors, edge_width, node_colors = [], [], []
    node_sizes = []
    color_values = []
    relabels = []

    for gene in list(G.nodes):

        if gene == "ZC3HC1":
            foldchange = file.loc[file.gene == gene, "log2FoldChange"].values[0]
            relabels.append(gene)
            node_sizes.append(1400)
            color_values.append(foldchange)
        elif gene == "CCNB1":
            relabels.append(gene)
            node_sizes.append(1400)
            color_values.append(0)
        elif gene in my_genes:
            foldchange = file.loc[file.gene == gene, "log2FoldChange"].values[0]
            color_values.append(foldchange)
            relabels.append(gene)
            node_sizes.append(400)
        elif gene in CCNB1_network_temp:
            relabels.append(gene)
            node_sizes.append(400)
            color_values.append(0)
        else:
            color_values.append(0)
            relabels.append("")
            node_sizes.append(20)

    node_colors, max_beta = get_colormap_for_values(color_values)

    temp = []
    print(list(cluster_dict.values()))
    for elem in list(cluster_dict.values())[1:]:
        temp.extend(elem)

    for gene_edges in list(G.edges):
        if "CCNB1" in gene_edges:
            edge_colors.append("gray")
            edge_width.append(2.4)
        elif gene_edges[0] in temp and gene_edges[1] in temp:
            for k, v in cluster_dict.items():
                if gene_edges[0] in v or gene_edges[1] in v:
                    edge_colors.append(cluster_colors_temp[k])
                    edge_width.append(1)
                    break
        else:
            edge_colors.append("silver")
            edge_width.append(0.3)

    print(len(edge_colors))
    print(len(list(G.edges)))
    edge_colors = edge_colors[:len(list(G.edges))]

    fig = plt.figure(1, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, right=0.8, top=1, bottom=0)
    ax = fig.add_subplot(111, aspect='auto')
    positions = nx.kamada_kawai_layout(G, dim=2)
    cluster_pos = get_cluster_center(positions, np.asarray(clusters_sorted))

    nx.draw(G,
            with_labels=True,
            labels=dict(zip(list(G.nodes), relabels)),
            pos=positions,  # default: kamada_kawai_layout
            width=edge_width,
            edge_color=edge_colors,
            node_color=node_colors,
            node_size=node_sizes,
            alpha=1,
            # fontdict={"size": 1, "weight": "bold"},
            # font_weight="bold",
            font_size=8)
    """
    cluster_name_ypos = 0.8
    for i in range(2, num_cluster + 1):
        cluster_name = "\n".join(enrichment_out_sum.loc[enrichment_out_sum.Cluster == i, "GO process"].values)
        xcenter = cluster_pos[i][0]
        ycenter = cluster_pos[i][1] + 0.02
        width = cluster_pos[i][2] + 0.13
        height = cluster_pos[i][3] + 0.16
        angle = 0
        e1 = patches.Ellipse((xcenter, ycenter), width, height,
                             angle=angle, linewidth=0,
                             fill=True, zorder=2,
                             facecolor=cluster_colors_temp[i - 1], alpha=0.00,
                             edgecolor=None)
        e2 = patches.Ellipse((xcenter, ycenter), width, height,
                             angle=angle, linewidth=3,
                             fill=False, zorder=2,
                             facecolor=None, alpha=0.8,
                             edgecolor=cluster_colors_temp[i - 1])
        plt.text(x=0.85, y=cluster_name_ypos, s=cluster_name, color=cluster_colors_temp[i - 1],
                 fontdict={"size": 8, "weight": "bold"}, ha="left")
        cluster_name_ypos -= 0.08
        ax.add_patch(e1)
        ax.add_patch(e2)
    """
    """ plt legend """
    legend_pos = "bottom right"
    max_beta = round(max_beta, 1)
    xticklabels_in = list(np.linspace(-max_beta, max_beta, num=5))
    xticks_in = list(np.linspace(0, 255, num=5))  # needs to be 255 not 256
    xticklabels_in[:] = [round(x, 1) for x in xticklabels_in]
    xticklabels_in[2] = 0
    width, height = 0.15, 0.01
    if legend_pos == "bottom left":
        left, bottom = 0.1, 0.1
    elif legend_pos == "bottom right":
        left, bottom = 0.75, 0.1
    elif legend_pos == "top left":
        left, bottom = 0.1, 0.9
    else:
        left, bottom = 0.4, 0.92
    ax = plt.axes([left, bottom, width, height],
                  yticks=[],
                  xticks=xticks_in,
                  xticklabels=xticklabels_in,
                  title="log$_{2}$(fold change)",
                  aspect='auto', in_layout=True)
    fontsize = 10
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(int(fontsize))
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(int(fontsize * 0.8))
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name="coolwarm"))

    fig2, ax = plt.subplots()
    ax.tick_params(axis="y", direction="in", pad=-10
                   )
    ax.tick_params(axis="x", direction="in"
                  )
    plt.yticks(ha="left")
    print(enrichment_out_sum)
    enrichment_out_sum = enrichment_out_sum[enrichment_out_sum.Cluster > 1]
    go_terms = enrichment_out_sum["GO process"].values
    go_terms[:] = [x.capitalize() for x in go_terms]
    plt.barh(y=enrichment_out_sum["GO process"],
             height=0.8,
             color=enrichment_out_sum["cluster_color"].values,
             width=enrichment_out_sum["p_value"].apply(lambda x: -1 * np.log10(x)),
             alpha=0.6
             )
    plt.ylim(len(enrichment_out_sum["GO process"]), -0.5)
    plt.ylabel("GO process")
    plt.xlabel("-log(P-value) (gene enrichment)")

    cluster_foldchange = get_cluster_foldchange_values(color_values, np.asarray(clusters_sorted))

    fig4, ax4 = plt.subplots()
    plt.subplots_adjust(right=0.4)
    ax4.tick_params(axis="y", direction="in")
    ax4.tick_params(axis="x", direction="in")

    my_palette = dict(zip(list(np.arange(0, num_cluster)), cluster_colors_temp[1: num_cluster + 1]))
    sns.violinplot(data=cluster_foldchange[1:], palette=my_palette, orient="horizontal",
                   alpha=0.6, width=1.1, inner=None, linewidth=0.0
                   )
    for idx, elem in enumerate(cluster_foldchange[1:]):
        plt.text(x=np.median(elem), y=idx, s=str(len(elem)), ha="center", va="center")
    plt.vlines(x=0, ymin=-0.5, ymax=num_cluster - 1.5, colors="gray", linestyles="dotted", linewidth=1.5)

    plt.xlabel("log$_{2}$(fold change)")
    plt.ylabel("Cluster")
    plt.yticks(list(np.arange(0, num_cluster - 1)), labels=list(np.arange(1, num_cluster)))
    plt.show()
    return


def get_cluster_foldchange_values(values, clusters):
    values = np.asarray(values)
    cluster_foldchange_ = []
    for i in range(0, num_cluster):
        cluster_foldchange_.append(list(values[clusters == i]))
    return cluster_foldchange_


def get_cluster_center(pos, clusters):
    pos = np.asarray(list(pos.values()))
    # print(pos)
    cluster_pos = {}
    for i in range(0, num_cluster):
        pos_temp = pos[clusters == i]
        cluster_pos[i + 1] = [np.mean(pos_temp[:, 0]),
                              np.mean(pos_temp[:, 1]),
                              abs(np.max(pos_temp[:, 0]) - np.min(pos_temp[:, 0])),
                              abs(np.max(pos_temp[:, 1]) - np.min(pos_temp[:, 1]))]
    # print(cluster_pos)
    return cluster_pos


def graph_to_edge_matrix(G, genes_in):
    """Convert a networkx graph into an edge matrix.
    See https://www.wikiwand.com/en/Incidence_matrix for a good explanation on edge matrices

    Parameters
    ----------
    G : networkx graph
    """
    # Initialize edge matrix with zeros
    edge_mat = np.zeros((len(G), len(G)), dtype=int)

    # Loop to set 0 or 1 (diagonal elements are set to 1)
    for node in G:
        for neighbor in G.neighbors(node):
            edge_mat[genes_in.index(node)][genes_in.index(neighbor)] = 1
        edge_mat[genes_in.index(node)][genes_in.index(node)] = 1

    return edge_mat


def get_cluster(graph):
    print("\n")
    genes_in = list(graph.nodes())
    """ https://www.learndatasci.com/tutorials/k-means-clustering-algorithms-python-intro/ """
    algorithms = {}
    algorithms['kmeans'] = cluster.KMeans(n_clusters=num_cluster, n_init=200)
    algorithms['agglom'] = cluster.AgglomerativeClustering(n_clusters=num_cluster, linkage="ward")
    algorithms['spectral'] = cluster.SpectralClustering(n_clusters=num_cluster, affinity="precomputed", n_init=50)
    algorithms['affinity'] = cluster.AffinityPropagation(damping=0.6)
    clustering = algorithms['agglom'].fit(graph_to_edge_matrix(graph, genes_in))
    clusters = (list(clustering.labels_))
    cluster_dict = {}
    for i in range(0, num_cluster):
        cluster_dict[i] = []
    for i in range(0, len(clusters)):
        if not genes_in[i] in cluster_dict[clusters[i]]:
            cluster_dict[clusters[i]].append(genes_in[i])
    sorted_clusters = sorted(cluster_dict, key=lambda k: len(cluster_dict[k]), reverse=True)
    cluster_dict_sorted = {}
    for i in range(0, num_cluster):
        cluster_dict_sorted[i] = cluster_dict[sorted_clusters[i]]
    clusters_sorted = []
    for i in range(0, len(clusters)):
        clusters_sorted.append(sorted_clusters.index(clusters[i]))
    return clusters_sorted, cluster_dict_sorted


def get_enrichment(gene_set):
    ##############################################################
    ## The following script retrieves and prints out
    ## significantly enriched (FDR < 1%) GO Processes
    ## for the given set of proteins.

    output_format = "json"
    method = "enrichment"
    request_url = "/".join([string_api_url, output_format, method])
    my_genes = gene_set

    params = {
        "identifiers":  "%0d".join(my_genes),  # your proteins
        "species": 9606,  # species NCBI identifier = Homo sapiens
        "caller_identity": "Tobias_Reinberger_Uni_Luebeck"  # your app name
    }

    response = requests.post(request_url, data=params)

    # Read and parse the results
    data = json.loads(response.text)

    df_out = json_normalize(data)
    df_out = df_out
    print(df_out["category"].unique())
    # print(df_out)

    df_out = df_out[(df_out["category"] == "Process")
                    | (df_out["category"] == "Function")
                    | (df_out["category"] == "SMART")
                    | (df_out["category"] == "InterPro")
                    | (df_out["category"] == "Pfam")
                    | (df_out["category"] == "RCTM")
                    | (df_out["category"] == "NetworkNeighborAL")
                    ]

    # df_out.to_excel("enrichment_NIPA_STRING_DB_2.xlsx", index=None)
    """
    for row in data:
        term = row["term"]
        preferred_names = ",".join(row["preferredNames"])
        fdr = float(row["fdr"])
        description = row["description"]
        category = row["category"]

        if category == "Process" and fdr < 0.01:
            #  print significant GO Process annotations
            print("\t".join([term, preferred_names, str(fdr), description]))
    """
    #df_out[df_out["category"] == "Process"]
    return df_out


def get_colormap_for_values(values):
    """
    :param values:
    :return colormap:
    """
    beta_colors = []
    rgb = list(cm.get_cmap(name='coolwarm')(np.arange(0, 256)))
    max_beta = max(values)
    min_beta = min(values)
    if abs(min_beta) > abs(max_beta):
        max_beta = abs(min_beta)
    for beta in values:
        temp_beta = beta / max_beta
        index = 127 + int(128 * temp_beta)
        if index == -1:
            index = 0
        beta_colors.append(rgb[index])
    return beta_colors, max_beta


if __name__ == '__main__':
    string_api_url = "https://string-db.org/api"

    file = pd.read_csv(r"Z:\Projects\Neointima_ZC3HC1_Cells\Redouane_NIPA_RNA_Seq\Treatment_Batch_Knockdown_results.csv",
                       )
    file2 = pd.read_csv(r"Z:\Projects\Neointima_ZC3HC1_Cells\Redouane_NIPA_RNA_Seq\NOTreatment_Batch_Knockdown_results.csv",
                        )

    # file = pd.concat([file, file2])
    # file = file2
    file.rename(columns={"Unnamed: 0": "gene"}, inplace=True)
    file.dropna(inplace=True)
    file.drop_duplicates(inplace=True)
    print(file)
    file = file[file["padj"] < 0.05]
    file = file[abs(file["log2FoldChange"]) > 0.4]
    file_under = file[file["log2FoldChange"] < -0.5]
    file_over = file[file["log2FoldChange"] > 0.5]

    my_genes = file["gene"].to_list()
    genes_underexpr = file_under["gene"].values
    genes_overexpr = file_over["gene"].values
    # a = get_enrichment(genes_overexpr)
    # b = get_enrichment(genes_underexpr)
    # print(a)
    # print(b)

    print(my_genes)
    known_interaction_partner = []
    gene_regulation = []

    num_cluster = 5

    # plot_string_image_for_each_gene()
    get_network_interactions()

