'''
Gene scoring description
'''

import random
import argparse as arg
from statistics import mean

'''
Added four arguments to parser: 
1. --GI parameter that takes gene interaction file path as a string 
2. --ef parameter that takes loci and genes at loci file path as a string
3. --nb parameter that takes the number of bins as an integer
4. --nt parameter that takes the number of trials as an integer
'''

parser = arg.ArgumentParser(description="Create two networks and determines whether their edge density is significantly "
                                        "different")
parser.add_argument("--GI", "-gene_interaction_file", type=str, help="gene interaction file path (default STRING.txt)",
                    default="STRING.txt")
parser.add_argument("--ef", "-experimental_file", type=str, help="experimental loci file path (default Input.gmt.txt)",
                    default="Input.gmt.txt")
parser.add_argument("--nb", "-n_bins", type=int, help="number of bins (default: 128)", default=128)
parser.add_argument("--nt", "-n_trials", type=int, help="number of trials (default: 1000)", default=1000)

args = parser.parse_args()


'''
Parameter f_experimental takes .gmt formatted file containing loci and genes at that loci.
Function opens file and reads each line creating a dictionary where keys are the loci number and the values are a list
of genes at that loci. 
Returns dictionary. 
'''
def dict_loci(f_experimental):
    f = open(f_experimental,"r")
    loci_gene_dict = {}
    loci = 0
    for line in f:
        line = line.strip("\n")
        line_list = line.split("\t")
        loci_gene_dict[loci] = line_list[2:]
        loci += 1
    f.close()
    return loci_gene_dict

'''
Parameter f_cofunction_net takes STRING database file containing gene interactions. 
dict_genes_interactions opens files and reads each line creating a dictionary of dictionaries where the keys are genes 
and the values are a list of genes they interact with. 
Returns dictionary.
'''
def dict_gene_interactions (f_cofunction_net):
    f = open(f_cofunction_net,"r")
    gene_interaction_dict = {}
    connected_gene_list = []
    for line in f:
        line = line.strip("\n")
        line = line.split("\t")
        if gene_interaction_dict.get(line[0]) is None:
            connected_gene_list = []
        connected_gene_list.append(line[1])
        gene_interaction_dict[line[0]] = connected_gene_list
    f.close()
    return gene_interaction_dict

'''
FA_subnetwork takes parameter f_experimental that string file path. 
Returns list of 12 strings. 
'''
def FA_subnetwork(f_experimental):
    loci_gene_dict = dict_loci(f_experimental)
    FAsubnet = []
    for loci in loci_gene_dict:
        rand_gene = random.choice(loci_gene_dict[loci])
        FAsubnet.append(rand_gene)
    return FAsubnet

'''
gene_score_dict takes parameter string file path
calls dict_loci function to create a dict of integers (keys) and list of strings (value), takes these lists and converts
them into keys for a new dictionary. 
returns a dict where keys are strings and values are empty lists
'''
def gene_score_dict(f_experimental):
    loci_gene_dict = dict_loci(f_experimental)
    gene_score_dict = {}
    for loci in loci_gene_dict:
        for gene in loci_gene_dict[loci]:
            gene_score_dict[gene] = []
    return gene_score_dict

'''
gene_num_edges takes two parameters: subnetwork a list of 12 strings, and stringdatabase a string file path 
returns a dictionary of string (keys) and a list of strings (values) of length 12 
'''
def subnet_connections(gene, subnetwork, gene_interactions):
    subnet_conn_dict = {}
    if gene in gene_interactions:
        total_list_conn = gene_interactions[gene]
        subnet_minusgene = subnetwork[:]
        subnet_minusgene.remove(gene)
        conn_list = []
        for conn_genes in subnet_minusgene:
            if conn_genes in total_list_conn:
                conn_list.append(conn_genes)
    else:
        conn_list = []
    subnet_conn_dict[gene] = conn_list
    return subnet_conn_dict


'''
subnet_edgecounts takes parameter of a dictionary of strings (keys) and list of strings (values)
Returns a dictionary of the same strings (keys) and 
'''
def subnet_edgecounts(subnet_conn_dict):
    subnet_counts = {}
    for key in subnet_conn_dict:
        subnet_counts[key] = len(subnet_conn_dict[key])
    return subnet_counts

'''
gene_scores takes two parameters both string file paths 
calls FA_subnetwork to generate list of 12 strings 
calculates gene score for every gene for that subnetwork 
'''
def gene_score(FAsubnetwork, gene_interactions):
    gene_scoring_dict = gene_score_dict(FAsubnetwork)
    OGsubnetwork = FA_subnetwork(FAsubnetwork)
    loci_genes = dict_loci(FAsubnetwork)
    for index in range(12):
        keep_gene = OGsubnetwork[index]
        for every_gene in loci_genes[index]:
            OGsubnetwork[index] = every_gene
            sc = subnet_connections(every_gene,OGsubnetwork, gene_interactions)
            sec = subnet_edgecounts(sc)
            gene_scoring_dict[every_gene].append(sec[every_gene])
        OGsubnetwork[index] = keep_gene
    return gene_scoring_dict


'''

'''
def n_trials_avgscore(num, FAsubnetwork, stringdatabase):
    gene_interactions = dict_gene_interactions (stringdatabase)
    nscores_dict = {}
    for trial in range(num):
        onetrial = gene_score(FAsubnetwork, gene_interactions)
        if trial == 0:
            nscores_dict = onetrial
        else:
            for key in onetrial:
                nscores_dict[key] += onetrial[key]
    for key in nscores_dict:
        avg = mean(nscores_dict[key])
        nscores_dict[key] = avg
    return nscores_dict



import time
FA_loci_filepath = "C:\\Users\\ascol\\OneDrive\\Desktop\\7711\\Module3Day3\\Input.gmt.txt"
stringdatabase = "C:\\Users\\ascol\\OneDrive\\Desktop\\7711\\Module3Day3\\STRING.txt"
t0 = time.time()
print(n_trials_avgscore(5000, FA_loci_filepath, stringdatabase))
t1 = time.time()
print(t1-t0)










