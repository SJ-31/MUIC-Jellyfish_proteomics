#################################################
#### ClueGO cyREST Python 3 Example Workflow ####
#################################################

# In Windows install first python for Windows 3.7.2 from https://www.python.org/downloads/release/python-372/
#### Please make sure that Cytoscape v3.6+' is started and the Cytoscape Apps 'yFiles Layout Algorithms' and 'ClueGO v2.5.2' are installed before running this script! ####

# required libraries that need to be installed before starting the script
# In Windows install first pip with e.g. py get-pip.py  (https://bootstrap.pypa.io/get-pip.py)
import requests # to install use "pip install requests" or "sudo apt install python3-request" in ubuntu linux
import json
import os
import csv
from pathlib import Path
from urllib.parse import quote
import time

# Helper functions
def read_gene_list(data):
    gene_ids=[]
    with open(data,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for row in reader:
            gene_ids.append(row[0])
    return gene_ids

def writeLines(lines,out_file):
    file = open(out_file,'w')
    for line in lines:
        file.write(line)
    file.close()

def writeBin(raw,out_file):
    file = open(out_file,'wb')
    file.write(raw)
    file.close()



#### choose the example to run with one or two gene lists ####
EXAMPLE_SELECTION = "ClueGO Rest Example for one gene list"
#EXAMPLE_SELECTION = "ClueGO Rest Example for two gene lists"

#### Basic settings for cyREST ####
# home.folder = "/home/user"
HOME_FOLDER = str(Path.home())

# Windows specific command to set home folder to user-home
HOME_FOLDER = HOME_FOLDER.replace("\\Documents","")

SEP = "/"
OUTPUT_FOLDER = HOME_FOLDER+SEP+EXAMPLE_SELECTION
try:
    os.stat(OUTPUT_FOLDER)
except:
    os.mkdir(OUTPUT_FOLDER)

CLUEGO_HOME_FOLDER = HOME_FOLDER+SEP+"ClueGOConfiguration"+SEP+"v2.5.3"
PORT_NUMBER = "1234"
HOST_ADDRESS = "localhost"
HEADERS = {'Content-Type': 'application/json'}


# define base urls
CYTOSCAPE_BASE_URL = "http://"+HOST_ADDRESS+":"+PORT_NUMBER+SEP+"v1"
CLUEGO_BASE_URL = CYTOSCAPE_BASE_URL+SEP+"apps"+SEP+"cluego"+SEP+"cluego-manager"

print("Example analysis Parameters for "+EXAMPLE_SELECTION)
print("User Home Folder: "+HOME_FOLDER)
print("User Output Folder: "+OUTPUT_FOLDER)
print("Cytoscape Base URL: "+CYTOSCAPE_BASE_URL)
print("ClueGO Base URL: "+CLUEGO_BASE_URL)

#### 0.0 Start up ClueGO in case it is not running yet
response = requests.post(CYTOSCAPE_BASE_URL+SEP+"apps"+SEP+"cluego"+SEP+"start-up-cluego", headers=HEADERS)
# wait 2 seconds to make sure ClueGO is started
if(str(response)=="<Response [500]>"):
    print("wait 2 secs")
    time.sleep(2) 

#### 1.0 Select the ClueGO Organism to analyze ####
organism_name = "Homo Sapiens" # (run "1.1 Get all ClueGO organisms" to get all options)
print("1.0 Select the ClueGO Organism to analyze: "+organism_name)
response = requests.put(CLUEGO_BASE_URL+SEP+"organisms"+SEP+"set-organism"+SEP+quote(organism_name), headers=HEADERS)

## [optional functions and settings, un comment and modify if needed]
#  # 1.1 Get all ClueGO organisms
#response <- requests.get(CLUEGO_BASE_URL+SEP+"organisms"+SEP+"get-all-installed-organisms")
#print(response.json())

# # 1.2 Get all info for installed organisms
#response <- requests.get(CLUEGO_BASE_URL+SEP+"organisms"+SEP+"get-all-organism-info")
#print(response.json())

#### 2.0 Upload IDs for a specific Cluster ####
print("2.0 Upload IDs for a specific Cluster")

cluster1 = 1
cluster2 = 2

# upload gene list for cluster 1
FILE_LOCATION = CLUEGO_HOME_FOLDER+SEP+"ClueGOExampleFiles/GSE6887_Bcell_Healthy_top200UpRegulated.txt"
gene_list1 = json.dumps(read_gene_list(FILE_LOCATION))
response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"upload-ids-list"+SEP+quote(str(cluster1)), data=gene_list1, headers=HEADERS)


# 2.1 Set the number of Clusters
max_input_panel_number = 1
if(EXAMPLE_SELECTION == "ClueGO Rest Example for two gene lists"):
    max_input_panel_number = 2

response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"max-input-panel"+SEP+str(max_input_panel_number))

if(EXAMPLE_SELECTION == "ClueGO Rest Example for two gene lists"):
    # upload gene list for cluster 2
    
    FILE_LOCATION = CLUEGO_HOME_FOLDER+SEP+"ClueGOExampleFiles/GSE6887_NKcell_Healthy_top200UpRegulated.txt"
    gene_list2 = json.dumps(read_gene_list(FILE_LOCATION))
    response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"upload-ids-list"+SEP+quote(str(cluster2)), data=gene_list2, headers=HEADERS)

# To add here in the same way if you have more than 2 cluster to compare



## [optional functions and settings, un comment and modify if needed]
## 2.2 Select the ClueGO ID type
#id_type_name = "# Automatic #" # (run "2.3 Get ClueGO ID types" to get all options)
#response = requests.put(CLUEGO_BASE_URL+SEP+"ids"+SEP+"set-id-type"+SEP+id_type_name, headers=HEADERS)

## 2.3 Refresh ClueGO source files
#requests.post(CLUEGO_BASE_URL+SEP+"ids"+SEP+"refresh-cluego-id-files", headers=HEADERS)

## 2.4 Get ClueGO ID types
#response = requests.get(CLUEGO_BASE_URL+SEP+"ids"+SEP+"get-all-installed-id-types", headers=HEADERS)
#print(response.json())


## [optional functions and settings, un comment and modify if needed]
# 2.5 Set analysis properties for a Cluster, this is to repeat for each input cluster
input_panel_index = cluster1 # set here the cluster number e.g. "1"
node_shape = "Ellipse" # ("Ellipse","Diamond","Hexagon","Octagon","Parallelogram","Rectangle","Round Rectangle","Triangle","V")
cluster_color = "#ff0000" # The color in hex, e.g. #F3A455
min_number_of_genes_per_term = 3
min_percentage_of_genes_mapped = 4
no_restrictions = False # "True" for no restricions in number and percentage per term
response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"set-analysis-properties"+SEP+str(input_panel_index)+SEP+node_shape+SEP+quote(cluster_color)+SEP+str(min_number_of_genes_per_term)+SEP+str(min_percentage_of_genes_mapped)+SEP+str(no_restrictions), headers=HEADERS)

if(EXAMPLE_SELECTION == "ClueGO Rest Example for two gene lists"):
    input_panel_index = cluster2 # set here the cluster number e.g. "2"
    node_shape = "Ellipse" # ("Ellipse","Diamond","Hexagon","Octagon","Parallelogram","Rectangle","Round Rectangle","Triangle","V")
    cluster_color = "#0000ff" # The color in hex, e.g. #F3A455
    min_number_of_genes_per_term = 3
    min_percentage_of_genes_mapped = 4
    no_restrictions = False # "True" for no restricions in number and percentage per term
    response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"set-analysis-properties"+SEP+str(input_panel_index)+SEP+node_shape+SEP+quote(cluster_color)+SEP+str(min_number_of_genes_per_term)+SEP+str(min_percentage_of_genes_mapped)+SEP+str(no_restrictions), headers=HEADERS)

# 2.6 Select visual style
if(EXAMPLE_SELECTION == "ClueGO Rest Example for two gene lists"):
    visual_style = "ShowClusterDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)
else:
    visual_style = "ShowGroupDifference"  # (ShowGroupDifference, ShowSignificanceDifference, ShowClusterDifference)

response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"select-visual-style"+SEP+visual_style, headers=HEADERS)


####  3.0 Select Ontologies
print("3.0 Select Ontologies")

selected_ontologies = json.dumps(["3;Ellipse","8;Triangle","9;Rectangle"]) # (run "3.1 Get all available Ontologies" to get all options)
response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-ontologies", data=selected_ontologies, headers=HEADERS)

## [optional functions and settings, un comment and modify if needed]
## 3.1 Get all available Ontologies
#response = requests.get(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"get-ontology-info", headers=HEADERS)
#print(response.json())

## 3.2 Select Evidence Codes
#evidence_codes = json.dumps(["All"]) # (run "3.3 Get all available Evidence Codes" to get all options)
#response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-evidence-codes", data=evidence_codes, headers=HEADERS)

## 3.3 Get all available Evidence Codes
#response = requests.get(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"get-evidence-code-info", headers=HEADERS)
#print(response.json())

# 3.4 Set min, max GO tree level
min_level = 5
max_level = 6
all_levels = False
response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-min-max-levels"+SEP+str(min_level)+SEP+str(max_level)+SEP+str(all_levels), headers=HEADERS)

## 3.5 Use GO term significance cutoff
#p_value_cutoff = 0.05
#use_significance_cutoff = True
#response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+str(use_significance_cutoff)+SEP+str(p_value_cutoff), headers=HEADERS)

## 3.6 Use GO term fusion
#use_go_term_fusion = True
#response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+str(use_go_term_fusion), headers=HEADERS)

## 3.7 Set statistical parameters
#enrichment_type = "Enrichment/Depletion (Two-sided hypergeometric test)" # ("Enrichment (Right-sided hypergeometric test)", "Depletion (Left-sided hypergeometric test)", "Enrichment/Depletion (Two-sided hypergeometric test)")
#multiple_testing_correction = True
#use_mid_pvalues = False
#use_doubling = False
#response = requests.put(CLUEGO_BASE_URL+SEP+"stats"+SEP+enrichment_type+SEP+str(multiple_testing_correction)+SEP+str(use_mid_pvalues)+SEP+str(use_doubling), headers=HEADERS)

## 3.8 Set the Kappa Score level
#kappa_score = 0.4
#response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies","set-kappa-score-level"+SEP+str(kappa_score), headers=HEADERS)

## 3.9 Set grouping parameters
#do_grouping = True
#coloring_type = "Random" # ("Random","Fix")
#group_leading_term = "Highest Significance" # ("Highest Significance","#Genes / Term","%Genes / Term","%Genes / Term vs Cluster")
#grouping_type = "Kappa Score" # ("Kappa Score","Tree")
#init_group_size = 1
#perc_groups_for_merge = 50
#perc_terms_for_merge = 50
#response = requests.put(CLUEGO_BASE_URL+SEP+"grouping"+SEP+str(do_grouping)+SEP+coloring_type+SEP+group_leading_term+SEP+grouping_type+SEP+str(init_group_size)+SEP+str(perc_groups_for_merge)+SEP+str(perc_terms_for_merge), headers=HEADERS)


#### 4.0 Run ClueGO Analysis ####
print("4.0 Run ClueGO Analysis")

# Run the analysis an save log file
analysis_name = EXAMPLE_SELECTION
analysis_option = "Cancel and refine selection" # ("Continue analysis","Skip the grouping","Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
response = requests.get(CLUEGO_BASE_URL+SEP+quote(analysis_name)+SEP+quote(analysis_option), headers=HEADERS)
log_file_name = OUTPUT_FOLDER+SEP+"ClueGO-Example-Analysis-log.txt"
writeLines(response.text,log_file_name)
#print(response.text)

# 4.1 Get network id (SUID) (CyRest function from Cytoscape)
response = requests.get(CYTOSCAPE_BASE_URL+SEP+"networks"+SEP+"currentNetwork", headers=HEADERS)
current_network_suid = response.json()['data']['networkSUID']
#print(current_network_suid)

print("Save results")

# Get network graphics (CyRest function from Cytoscape)
image_type = "svg" # png, pdf
response = requests.get(CYTOSCAPE_BASE_URL+SEP+"networks"+SEP+str(current_network_suid)+SEP+"views"+SEP+"first."+image_type)
image_file_name = OUTPUT_FOLDER+SEP+"ClueGOExampleNetwork."+image_type
writeBin(response.content,image_file_name)

# 4.2 Get ClueGO result table
response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-cluego-table"+SEP+str(current_network_suid))
table_file_name = OUTPUT_FOLDER+SEP+"ClueGO-Example-Result-Table.txt"
writeLines(response.text,table_file_name)

# 4.3 Get ClueGO genes and main functions
number_of_functions_to_add = 3
response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-main-functions"+SEP+str(current_network_suid)+SEP+str(number_of_functions_to_add))
table_file_name = OUTPUT_FOLDER+SEP+"ClueGO-Genes-With-Main-Functions.txt"
writeLines(response.text,table_file_name)

# 4.4 Get genes result table
response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-gene-table"+SEP+str(current_network_suid))
table_file_name = OUTPUT_FOLDER+SEP+"ClueGO-Gene-Table.txt"
writeLines(response.text,table_file_name)

# 4.5 Get Kappascore Matrix
response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-kappascore-matrix"+SEP+str(current_network_suid))
table_file_name = OUTPUT_FOLDER+SEP+"ClueGO-Kappascore-Matrix.txt"
writeLines(response.text,table_file_name)

# 4.6 Get binary Gene-Term Matrix
response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-binary-gene-term-matrix"+SEP+str(current_network_suid))
table_file_name = OUTPUT_FOLDER+SEP+"ClueGO-Binary-Gene-Term-Matrix.txt"
writeLines(response.text,table_file_name)

# 4.7 ClueGO Result Chart
# Get result charts for both cluster as pie chart
chart_type = "PieChart" # ("PieChart","BarChart")
image_type = "svg" # ("svg","png","pdf")
response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-cluego-result-chart"+SEP+str(current_network_suid)+SEP+str(cluster1)+SEP+chart_type+SEP+image_type)
image_file_name = OUTPUT_FOLDER+SEP+"ClueGO-"+chart_type+"-For-Cluster"+str(cluster1)+"."+image_type
writeBin(response.content,image_file_name)

if(EXAMPLE_SELECTION == "ClueGO Rest Example for two gene lists"):
    response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-cluego-result-chart"+SEP+str(current_network_suid)+SEP+str(cluster2)+SEP+chart_type+SEP+image_type)
    image_file_name = OUTPUT_FOLDER+SEP+"ClueGO-"+chart_type+"-For-Cluster"+str(cluster2)+"."+image_type
    writeBin(response.content,image_file_name)


### [optional functions and settings, un comment and modify if needed]
## 4.8 Remove ClueGO analysis result
#print("Remove ClueGO Network")
## Remove analysis to reduce memory usage. This is important when using patch modes that create lots of analyses.
#response = requests.delete(CLUEGO_BASE_URL+SEP+"remove-cluego-analysis-result"+SEP+str(current_network_suid))

## 4.9 Remove all ClueGO analysis results
#print("Remove all ClueGO Networks")
#response = requests.delete(CLUEGO_BASE_URL+SEP+"remove-all-cluego-analysis-results")

print("Example analysis "+EXAMPLE_SELECTION+" done")
