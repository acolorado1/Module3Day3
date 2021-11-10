'''
GeneScoring is a script that takes a .gmt formatted file and a STRING.txt database and calculates a score of each gene
found in the the .gmt formatted file. The score is calculated by taking the average number of connections that the gene
has in n_trials subnetworks. A .sif file is then written to the specified directory wherein the genes containing the
highest scores in each loci make up the subnetwork. Visualization is done in Cytoscape.
'''

import random
import argparse as arg
from statistics import mean

'''
Added four arguments to parser: 
1. --nt parameter that takes the number of trials as an integer
2. --gmt parameter that takes loci and genes at loci file path as a string
3. --sdb parameter that takes string database file path as a string
4. --sfd parameter that takes new file path as an string
'''

parser = arg.ArgumentParser(description="Calculate average gene scores and visualize networks of the genes with the "
                                        "highest scores from each loci" )
parser.add_argument("--nt", "-n_trials", type=int, help="number of trials to run",
                    default=5000)
parser.add_argument("--gmt", "-gmt_formated_file", type=str, help="experimental loci file path (default Input.gmt.txt)",
                    default="Input.gmt.txt")
parser.add_argument("--sdb", "-string_database", type=str, help="gene interaction file path (default STRING.txt)",
                    default="STRING.txt")
parser.add_argument("--sfd", "-SIF_file_directory", type=str, help="file path of new sif file (default subnetwork.sif",
                    default="subnetwork.sif")

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
gene_num_edges takes two parameters: subnetwork a list of 12 strings, gene the current gene of interest and 
gene_interactions a dictionary of genes (keys) and lists of genes that they interact with (values). 
returns a dictionary of one gene string (keys) and a list of strings (values) genes that they are connected to within 
the network. 
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
Returns a dictionary of the same string (keys) and a integer as the value. 
'''
def subnet_edgecounts(subnet_conn_dict):
    subnet_counts = {}
    for key in subnet_conn_dict:
        subnet_counts[key] = len(subnet_conn_dict[key])
    return subnet_counts

'''
gene_scores takes two parameters both string file paths 
calls gene_score_dict function to generate an empty dictionary of keys containing all the genes from the .gmt file 
calls FA_subnetwork to generate a list of 12 random genes (strings) containing one gene from each of the 12 loci 
calls dict_loci to generate a dictionary where keys are integers denoting loci and values are lists of genes contained 
in each loci. 
for every loci (there are 12), the gene at that loci in the subnetwork generated is kept 
for every_gene in the loci replace the gene for that loci with every_gene 
    get the list of genes every_gene interacts with in the subnetwork 
    calculate the number of genes every_gene is connected to 
    append to gene_scoreing_dict the number of connected genes (integer) as list value for every_gene
replace the gene at loci in the subnetwork with the original gene chosen 
return dictionary of genes (keys) and number of genes they were connected to in the subnetwork as a list (values)
'''
def gene_score(FAsubnetwork, gene_interactions, numL):
    gene_scoring_dict = gene_score_dict(FAsubnetwork)
    OGsubnetwork = FA_subnetwork(FAsubnetwork)
    loci_genes = dict_loci(FAsubnetwork)
    for index in range(numL):
        keep_gene = OGsubnetwork[index]
        for every_gene in loci_genes[index]:
            OGsubnetwork[index] = every_gene
            sc = subnet_connections(every_gene,OGsubnetwork, gene_interactions)
            sec = subnet_edgecounts(sc)
            gene_scoring_dict[every_gene].append(sec[every_gene])
        OGsubnetwork[index] = keep_gene
    return gene_scoring_dict

'''
n_trials_avgscore takes three parameters: num, an integer representing the number of trials to be run, FAsubnetwork and 
string database both which are string file paths 
they interact with (values) 
for every trial 
    a gene score is calculated for every gene located in all the loci provided 
    if it is the first trial key value is initializes as a list of one integer 
    if not the first trial following lists of integers are concatenated together
results in a dictionary of string genes (keys) with list of num integers (value) 
for every gene in the dictionary 
    average of their value a list is taken 
    value in dictionary is replaced with the average score 
dictionary of genes (keys) and average gene scores from num trials (values) is returened 
'''
def n_trials_avgscore(num, FAsubnetwork, gene_interactions ,numL):
    nscores_dict = {}
    for trial in range(num):
        onetrial = gene_score(FAsubnetwork, gene_interactions, numL)
        if trial == 0:
            nscores_dict = onetrial
        else:
            for key in onetrial:
                nscores_dict[key] += onetrial[key]
    for key in nscores_dict:
        avg = mean(nscores_dict[key])
        nscores_dict[key] = avg
    return nscores_dict

'''
Visualizing results portion of the code:

createssiffile takes two parameters: conn_dict a dictionary of genes (keys) and list of genes they're connected to
(values) and the filepath of the files being written. 
writes a .sif file to be used in cytoscpe to visualize subnetworks 
'''
def createsiffile(conn_dict, newdirectory):
    sif_file = open(newdirectory, 'w+')
    for key in conn_dict:
        listofstrings = [str(item) for item in conn_dict[key]]
        stringlist = '\t'.join(listofstrings)
        row = "{0}\tpp\t{1}\n".format(str(key), stringlist)
        sif_file.write(row)
    sif_file.close()

'''
VisualizeGeneScoring takes four parameters: 
num (integer) -> number of trials to run 
FAsubnetwork (string) -> .gmt file path 
stringdatabase (string) -> STRING database file path 
newdirectory (string) -> file path where .sif file will be written to 
calls dict_gene_interactions to generate dictionary of all genes in databased (keys) and list of connected genes (values)
calls n_trials_avgscore to create a dict of all genes in all loci (keys) and average gene score for num trials (value) 
calls loci_genes to create dictionary of loci (keys) and list of genes at loci (values) 
empty list is initialized 
for every loci in loci_genes
    gene with highest average gene score is found 
    gene is appended to list 
empty dictionary is initialized 
for every gene in in high_subnetwork list 
    list of connected genes is found 
    gene is added as a key and list of connected genes is its value to empty dictionary 
createsiffile is called on the dictionary and .sif is written to specified directory 
'''
def VisualizeGeneScoring (num, FAsubnetwork, stringdatabase, newdirectory):
    gene_interactions = dict_gene_interactions(stringdatabase)
    loci_genes = dict_loci(FAsubnetwork)
    gene_scores = n_trials_avgscore(num, FAsubnetwork, gene_interactions, len(loci_genes))
    high_subnetwork = []
    for loci in loci_genes:
        highest_gene_score = 0
        highest_gene = ''
        for gene in loci_genes[loci]:
            if gene_scores[gene] > highest_gene_score:
                highest_gene_score = gene_scores[gene]
                highest_gene = gene
        high_subnetwork.append(highest_gene)
    high_gene_conn_dict = {}
    for gene in high_subnetwork:
        gene_conn_dict = subnet_connections(gene, high_subnetwork, gene_interactions)
        high_gene_conn_dict[gene] = gene_conn_dict[gene]
    createsiffile(high_gene_conn_dict, newdirectory)


VisualizeGeneScoring(args.nt, args.gmt, args.sdb, args.sfd)
