#!/usr/bin/env python3

#Importing of librarys
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
from multiprocessing import cpu_count
from shutil import copyfile
import re
import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import datetime
import time
from scipy import stats
from scipy.cluster import hierarchy
from warnings import simplefilter

#Ignoring some warnings
simplefilter("ignore", FutureWarning)
simplefilter("ignore", UserWarning)
simplefilter("ignore", hierarchy.ClusterWarning)

#Setting version
version = "1.00"

#Counting how much cpu core are in PC
num_threads = cpu_count()

#Setting a help message
help = """
MPACT - Multimetric PAirwise Comparison Tool version """+version+""" - 24 jan 2025
This program perform all-against-all pairwise comparisons for biological sequences 
and protein structures with five diferent metrics. The results are presented in 
heatmap graph, frequency distribution, trees and matrices.
(c) 2024. Igor Custódio dos Santos & Arthur Gruber

Usage: MPACT -i <input_file> -s <structures> -p <protocols> -o <output_directory>

Mandatory parameters:

-i input <file name>       	Input file (FASTA sequence file)

-s PDBs <directory name>	3D-structures directory (PDBs directory)

-p protocol <integer> (default=all methods)
      	Metric options		0 - pairwise aa sequence identity
				1 - pairwise aa sequence similarity
				2 - maximum likelihood (ML)
				3 - pairwise structural alignment
				4 - pairwise 3Di-character alignment
	Multiple metrics - Examples:
				All methods: -m 0,1,2,3,4
				Combined ML and structural alignment: -m 2,3

-o output <directory name>	Output directory to store files

Additional parameters:

-conf <configuration file name> (default=none)
                                Configuration file in which parameters are specified

-c color			Color palette for graphs (matplotlib options)
				(default=RdBu)

-f mafft <parameters>		Only "--maxiterate" and "--thread" MAFFT parameters
				are accepted
				e.g. -f "--maxiterate 1000 --thread 20"

-g create subgroups <integer>,<lower limit>,<upper limit> (default=none)
      	Metric options		0 - pairwise aa sequence identity (% value)
				1 - pairwise aa sequence similarity (% value)
				2 - maximum likelihood (ML distance value)
				3 - pairwise structural alignment (TM-score value)
				4 - pairwise 3Di-character alignment (% value)
	Examples:
				ML 7.5-9,0: -s 2,7.5,9.0
				Similarity 70-85%: -s 2,70,85
				TM-score 0.75-1.0: -s 3,0.75,1

-h help                       	Print this page of help

-k <clustering method> (default=average)
                                The clustering method to be used on heatmap.
                                Use scipy.cluster.hierarchy.linkage method ("average",
                                "single","complete","weighted","centroid","median",
                                "ward") or neighbor-joining ("nj")

-l lines <yes|no>		Use cell dividing lines on heatmap (default=no)

-n names* <id|org|all>		Names for labels. Examples:
                                  id - YP_003289293.1 (default)
                                  org - Drosophila totivirus 1
                                  all - YP_003289293.1 - Drosophila totivirus 1

-q iqtree <parameters>		Only "-m", "-bb" and "-nt" IQTREE parameters are
				accepted
				e.g. -q "-m Q.pfam+F+R6 -bb 1000 -nt 20"

-t tags <prefixes>		Prefixes of already known groups to obtain their
                                intra- and inter-group value ranges. Separate
                                the prefixes by commas.
				e.g. -t "Group1,Group2,Group3"

-v version			Show the program's version.

*Only valid for NCBI naming format:
YP_003289293.1 RNA-dependent RNA polymerase [Drosophila totivirus 1]
For other naming formats, all terms are used
"""

#Setting parser and all parameters with argparse
parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument('-conf')
parser.add_argument('-s')
parser.add_argument('-p', default='0,1,2,3,4')
parser.add_argument('-o')
parser.add_argument('-c', default='RdBu')
parser.add_argument('-k', default='nj')
parser.add_argument('-l', default='no')
parser.add_argument('-n', default="id")
parser.add_argument('-g', default = None)
parser.add_argument('-t', default = None)
parser.add_argument('-h', '--help', action='store_true')
parser.add_argument('-v', '--version', action='store_true')
parser.add_argument('-f', default='--maxiterate 1000 --thread '+str(round(num_threads/2)))
parser.add_argument('-q', default='-bb 1000 -nt '+str(round(num_threads/2)))
args = parser.parse_args()

#Defining accessory functions.

def no_empty(x):
  """It removes empty strings of a list."""
  while '' in x:
    x.remove('')
  return x

def can_be_float(string):
  """It checks if a string can be converted into a float number."""
  try:
    float(string)
    return True
  except ValueError:
    return False

def isInt(value):
  """It checks if a string object can be converted into a integer number."""
  try:
    int(value)
    return True
  except:
    return False




def configuration_file(args):
  """From a configuration file, this function sets all the parameters it specifies."""

  with open(args.conf, "r") as read_conf:
    conf_lines = read_conf.readlines()

  input_fasta = args.i
  pdb_directory = args.s
  protocols = args.p
  output_directory = args.o
  heatmap_pallette = args.c
  heatmap_lines = args.l
  heatmap_labels = args.n
  subgroups_formation = args.g
  mafft_parameters = args.f
  iqtree_parameters = args.q
  prefixes = args.t
  clustering_method = args.k

  for line in conf_lines:
    if line.strip() in ["", "\n"]:
      continue
    letter = line.split("=")[0].split()[0]
    argument = line.split("=")[1].split()[0]

    #For input FASTA file
    if letter == "i":
      input_fasta = argument

    #For input structures in PDB diretory
    if letter == "s":
      pdb_directory = argument

    #For protocols
    if letter == "p":
      protocols = argument

    #For output directory
    if letter == "o":
      output_directory = argument

    #For colour pallette of heatmaps
    if letter == "c":
      heatmap_pallette = argument

    #For lines in heatmaps
    if letter == "l":
      heatmap_lines = argument

    #For names/labels in heatmaps
    if letter == "n":
      heatmap_labels = argument

    #For clustering in subgroups
    if letter == "g":
      subgroups_formation = argument

    #For MAFFT program parameters
    if letter == "f":
      mafft_parameters = argument

    #For IQ-TREE program parameters
    if letter == "q":
      iqtree_parameters = argument

    #For prefixes to obtain ranges
    if letter == "t":
      prefixes = argument

    #For clustering method
    if letter == "k":
      clustering_method = argument

  return (input_fasta, pdb_directory,
          protocols, output_directory,
          heatmap_pallette, heatmap_lines,
          heatmap_labels, subgroups_formation,
          mafft_parameters, iqtree_parameters, prefixes, clustering_method)


def check_fasta(fasta):
  """It verify that if the fasta file is valid."""
  if fasta == None:
    return "not_ncbi", "protein"
  elif os.path.isdir(fasta):
    print("Error: Input (-i) is a directory.")
    sys.exit()
  else:
    if os.path.isfile(fasta):
      with open(fasta, 'r') as file:
        lines = file.readlines()
      if not lines:
        print(f"Error: Input file ({fasta}) is empty.")
        sys.exit()
      if not lines[0].startswith(">"):
        print(f"Error: Input file ({fasta}) is not in FASTA format.")
        sys.exit()
    else:
      print(f"Error: Input file ({fasta}) not found.")
      sys.exit()

    #Determine that if the headers of the fasta are in the ncbi standard or not.
    type_header = "ncbi"
    for line in lines:
      if line.startswith(">") and not re.match(r'>.* .* \[.*\]', line):
        type_header= "not_ncbi"
        break

    #Determine that if the sequences of the fasta are nucleotide or protein.
    is_protein = None
    for line in lines:
      line = line.strip()
      if line.startswith(">"):
        continue
      elif any(base in "FEILfeil" for base in line.upper()):
        is_protein = True
        break
      else:
        is_protein = False

    #Return the result of this function.
    if is_protein:
      return type_header, "protein"
    elif is_protein is not None:
      print("NUCLEOTIDE")
      return type_header, "nucleotide"

def check_methods(args):
    """Check the methods/protocols parameter (-p)"""
    pattern = r"[0-4](\,[0-4])*"
    if re.search(pattern, args.p):
      metric_list = args.p.split(",")
      nrecognmethods = ""
      for m in metric_list:
        if m not in "01234":
          nrecognmethods += m
    else:
      print("Metrics sintax (-p parameter) is incorrect. Please correct.")

    if len(nrecognmethods) != 0:
      print(f"Methods '{nrecognmethods}' not recognized.")

    return metric_list

def check_clustering_method(cluster_method):
  """Check the clustering method parameter (-k)"""
  valid_clustering_methods = ["average","single","complete","weighted","centroid","median","ward", "nj"]
  if cluster_method not in valid_clustering_methods:
    print(f"ERROR: Clustering method '{args.k}' not supported.")
    sys.exit()

def check_prefixes(args):
  """Check the prefixes parameter (-t)"""
  prefixes = args.t.split(",")

  nrecogn = []
  if any(x in args.p for x in ["0", "1", "2"]):
    with open(args.i, "r") as read_fasta:
      fasta = read_fasta.readlines()
    headers = list(x for x in fasta if x.startswith(">"))

    for prefix in prefixes:
      b = 0
      for header in headers:
        if header.startswith(">"+prefix):
          b=1
          break
      if b==0:
         nrecogn.append(prefix)
      else:
        continue
    if len(nrecogn) != 0:
      string_nrecogn = ", ".join(nrecogn)
      print(f"Prefixes '{string_nrecogn}' not found on FASTA file.")
      sys.exit()

  if any(x in args.p for x in ["3", "4"]):
    pdbs_list = os.listdir(args.s)

    for prefix in prefixes:
      b = 0
      for pdb in pdbs_list:
        if pdb.startswith(prefix):
          b=1
          break
      if b==0:
         nrecogn.append(prefix)
      else:
        continue
    if len(nrecogn) != 0:
      string_nrecogn = ", ".join(nrecogn)
      print(f"Prefixes '{string_nrecogn}' not found on PBDs directory.")
      sys.exit()

  return prefixes


def check_pdb(pdb_dict):
  """Verify that if the PDB files directory exists"""
  if pdb_dict != None:
    if os.path.isdir(pdb_dict) == False:
      print(f"\nERROR: '{pdb_dict}' (-s) is not a directory.")
      sys.exit()

    file_list = os.listdir(pdb_dict)

    #If it is not empty;
    if len(file_list) == 0:
      print(f"\nERROR: The directory '{pdb_dict}' (-s) is empty.")
      sys.exit()

    #And if all files end with ".pdb":
    for pdb in file_list:
      if pdb.endswith(".pdb")==False:
        print(f'In "{pdb_dict}" directory, the file "{pdb}" might not be a pdb format file.\nPlease take it off from the input directory or rename it with pdb extention.')
        sys.exit()

def check_output_dir(output):
  """Add a new directory if the user don't provide one."""
  if output != None:
    return output
  else:
    i = 1
    while True:
      if os.path.isdir(f"output{str(i)}"):
        i +=1
        continue
      else:
        save_output_dir = f"output{str(i)}"
        break
    print(f"Output directory not specified. Saving files in {save_output_dir}")

  return save_output_dir


def check_heatmap_param(args):
  """Checks the heatmap related parameters."""

  #If...
  #Argument -l is valid;
  if args.l not in ['yes', 'no']:
    print("Error: Line parameter (-l) not recognized.")
    sys.exit()

  #Argument -c is valid;
  color = args.c
  if not sns.color_palette(color):
    print("Error: Color parameter (-c) not recognized.")
    sys.exit()

  #Argument -n is valid;
  if args.n not in ['id', 'org', 'all']:
    print("Error: Label parameter (-n) not recognized.")
    sys.exit()



def check_mafft(param_mafft):
  """Check the MAFFT program parameter (-f)"""
  num_threads = cpu_count()
  mafft_list = param_mafft.split(' ')
  if len(mafft_list) % 2 != 0:
    print("ERROR: MAFFT parameters not accepted. It has an odd number of parameters.")
    sys.exit()
  for i in range(0, len(mafft_list), 2):
    if mafft_list[i] == "--maxiterate":
      if not isInt(mafft_list[i+1]):
        print("ERROR: MAFFT parameter --maxiterate is not integer")
        sys.exit()
    elif mafft_list[i] == "--thread":
      if not isInt(mafft_list[i+1]):
        print("ERROR: MAFFT parameter --thread is not integer")
        sys.exit()
      else:
        if int(mafft_list[i+1]) > num_threads:
          print("Number of processors ("+mafft_list[i+1]+") exceeds the total CPUs available on the server. Using "+str(round(num_threads))+" CPUs.")
          mafft_list[i+1] = str(round(num_threads))
    else:
      print("ERROR: MAFFT parameter '"+mafft_list[i]+"' not accepted by 3A-DGT. Please correct your syntax.")
      sys.exit()
  if len(mafft_list) == 2:
    if "--maxiterate" not in mafft_list:
      mafft_list.append("--maxiterate")
      mafft_list.append("1000")
    if "--thread" not in mafft_list:
      mafft_list.append("--thread")
      mafft_list.append(str(round(num_threads/2)))

  mafft_cmd = " ".join(mafft_list)

  return mafft_cmd

def check_clustering(clust_cmd):
  """Check the data partitioning parameter (-g)"""
  if "," not in clust_cmd:
    print("ERROR: Clustering command not valid. Please review the syntax.")
    sys.exit()

  clust_cmd_list = clust_cmd.split(",")

  if len(clust_cmd_list) != 3:
    print("ERROR: Clustering command not valid. Please review the syntax.")
    sys.exit()

  if can_be_float(clust_cmd_list[1]) == False:
    print("ERROR: Specified lower limit is not a number. Please review the syntax.")
    sys.exit()

  if can_be_float(clust_cmd_list[2]) == False:
    print("ERROR: Specified upper limit is not a number. Please review the syntax.")
    sys.exit()

  lower_limit = float(clust_cmd_list[1])
  upper_limit = float(clust_cmd_list[2])

  if upper_limit < lower_limit:
    print("ERROR: Lower limit is higher than upper limit. Please review the syntax.")
    sys.exit()


  if clust_cmd_list[0] == "0":
    if lower_limit < 0 or lower_limit > 100:
      print("ERROR: Specified lower limit is out of identity range. Choose a number between 0 and 100 (%).")
      sys.exit()
    if upper_limit < 0 or upper_limit > 100:
      print("ERROR: Specified upper limit is out of identity range. Choose a number between 0 and 100 (%).")
      sys.exit()
  elif clust_cmd_list[0] == "1":
    if lower_limit < 0 or lower_limit > 100:
      print("ERROR: Specified lower limit is out of similarity range. Choose a number between 0 and 100 (%).")
      sys.exit()
    if upper_limit < 0 or upper_limit > 100:
      print("ERROR: Specified upper limit is out of similarity range. Choose a number between 0 and 100 (%).")
      sys.exit()

  elif clust_cmd_list[0] == "2":
    if lower_limit < 0 or lower_limit > 10:
      print("ERROR: Specified lower limit is out of maximum-likelihood distance range. Choose a number between 0 and 10.")
      sys.exit()
    if upper_limit < 0 or upper_limit > 10:
      print("ERROR: Specified upper limit is out of maximum-likelihood distance range. Choose a number between 0 and 10.")
      sys.exit()

  elif clust_cmd_list[0] == "3":
    if lower_limit < 0 or lower_limit > 1:
      print("ERROR: Specified lower limit is out of TM-scores range. Choose a number between 0 and 1.")
      sys.exit()
    if upper_limit < 0 or upper_limit > 1:
      print("ERROR: Specified upper limit is out of TM-scores range. Choose a number between 0 and 1.")
      sys.exit()

  elif clust_cmd_list[0] == "4":
    if lower_limit < 0 or lower_limit > 100:
      print("ERROR: Specified lower limit is out of 3Di characters similarity range. Choose a number between 0 and 100 (%).")
      sys.exit()
    if upper_limit < 0 or upper_limit > 100:
      print("ERROR: Specified upper limit is out of 3Di characters similarity range. Choose a number between 0 and 100 (%).")
      sys.exit()

  else:
    print("ERROR: Method of clustering command not recognized. Please review the syntax.")
    sys.exit()

  print("The sequence subset generation parameter is correct.")
  return clust_cmd_list


def check_iqtree(param_iqtree):
  """Check the IQ-TREE parameter (-n)"""
  num_threads = cpu_count()
  iqtree_list = param_iqtree.split(' ')
  if len(iqtree_list) % 2 != 0:
    print("ERROR: IQ-TREE parameters not accepted. It has an odd number of parameters.")
    sys.exit()
  for i in range(0, len(iqtree_list), 2):
    if iqtree_list[i] == "-bb":
      if not isInt(iqtree_list[i+1]):
        print("ERROR: IQ-TREE parameter -bb is not integer.")
        sys.exit()
      if isInt(iqtree_list[i+1]) and int(iqtree_list[i+1]) < 1000:
        print("ERROR: IQ-TREE parameter -bb must be greater than 1000.")
        sys.exit()
    elif iqtree_list[i] == "-nt":
      if not isInt(iqtree_list[i+1]):
        print("ERROR: IQ-TREE parameter -nt is not integer")
        sys.exit()
      else:
        if int(iqtree_list[i+1]) > num_threads:
          print("Number of processors ("+iqtree_list[i+1]+") exceeds the total CPUs available on the server. Using "+str(round(num_threads))+" CPUs.")
          iqtree_list[i+1] = str(round(num_threads))
    elif iqtree_list[i] == "-m":
      continue
    else:
      print("ERROR: IQ-TREE parameter '"+iqtree_list[i]+"' not accepted by 3A-DGT. Please correct your syntax.")
      sys.exit()

  if "-bb" not in iqtree_list:
    iqtree_list.append("-bb")
    iqtree_list.append("1000")
  if "-nt" not in iqtree_list:
    iqtree_list.append("-nt")
    iqtree_list.append(str(round(num_threads/2)))

  iqtree_cmd = " ".join(iqtree_list)
  return iqtree_cmd

def log_header(log, args, type_sequence, type_header, version):
  """Open the logfile and write a header with run informations"""
  #Setting date time for logfile
  now = datetime.datetime.now()
  date_format = "%m-%d-%Y %H:%M:%S"
  formatted_date = now.strftime(date_format)

  #Preparing configuration file informations for logfile
  if args.conf:
    with open(args.conf, "r") as read_conf_file:
      conf_file_content = "".join(read_conf_file.readlines())
    log_conf_file = f" {args.conf}\n\n*************** Configuration file content ***************\n"
    log_conf_file += f"{conf_file_content}\n**********************************************************"
  else:
    log_conf_file = "NOT SPECIFIED"

  #Preparing FASTA informations for logfile
  if args.i != None:
    fasta_name = args.i

    count_fasta = 0
    with open(args.i, 'r') as fasta:
      line = fasta.readline()
      while line:
        if line.startswith('>'):
          count_fasta += 1
        line = fasta.readline()

    if type_sequence == "protein":
      log_type_sequence = "Protein"
    elif type_sequence == "nucleotide":
      log_type_sequence = "Nucleotide"

    if type_header == "ncbi":
      log_type_header = "NCBI standard"
    elif type_header == "not_ncbi":
      log_type_header = "Not NCBI standard"

  else:
    fasta_name = "NOT SPECIFIED"
    count_fasta = " - "
    log_type_sequence = " - "
    log_type_header = " - "

  #Preparing PDB directory informations for logfile
  if args.s != None:
    pdb_name = args.s
    count_pdb = len(os.listdir(args.s))
  else:
    pdb_name = "NOT SPECIFIED"
    count_pdb = " - "


  #Preparing methods informations for logfile
  method_dict = {"0": "(0) pairwise sequence alignment identity",
                 "1": "(1) pairwise sequence alignment similarity",
                 "2": "(2) maximum likelihood distance",
                 "3": "(3) pairwise structural alignment - TM-score",
                 "4": "(4) pairwise 3di character alignment similarity",
                 }
  methods = ""
  for met in method_dict:
    if met in args.p:
      methods += f"\n  - {method_dict[met]}"

  #Preparing clustering methods informations for logfile
  if args.k == "nj":
    log_clust_method = "Neighbor-joining (nj)"
  else:
    log_clust_method = f"{args.k} (from scipy.cluster.hierarchy.linkage) and Neighbor-joining (nj)"

  #Preparing heatmap label informations for logfile
  header_dict ={"id": "Identification code", "org": "Organism", "all": "All the header"}
  if args.n in header_dict:
    log_names = f"{header_dict[args.n]} ({args.n})"

  #Preparing MAFFT parameters informations for logfile
  if args.f == '--maxiterate 1000 --thread '+str(round(num_threads/2)):
    log_mafft = args.f+" (default)"
  else:
    log_mafft = args.f

  #Preparing IQ-TREE parameters informations for logfile
  if args.q == '-bb 1000 -nt '+str(round(num_threads/2)):
    log_iqtree = args.q+" (default)"
  else:
    log_iqtree = args.q

  #Preparing data partitioning informations for logfile
  method_partit_dict = {"0": "Identity",
                 "1": "Similarity",
                 "2": "Maximum likelihood distance",
                 "3": "TM-score",
                 "4": "3Di similarity",
                 }
  if args.g:
    method_partitioning = args.g.split(",")[0]
    lower_limit = args.g.split(",")[1]
    upper_limit = args.g.split(",")[2]
    if method_partitioning in ["0", "1", "4"]:
      log_data_partitioning = f"{args.g} ({lower_limit}% to {upper_limit}% of {method_partit_dict[method_partitioning]})"
    elif method_partitioning in ["2", "3"]:
      log_data_partitioning = f"{args.g} ({lower_limit} to {upper_limit} of {method_partit_dict[method_partitioning]})"
  else:
    log_data_partitioning = "NOT SPECIFIED"

  #Preparing prefixes informations for logfile
  if args.t:
    log_prefixes = args.t
  else:
    log_prefixes = "NOT SPECIFIED"

  log.write(f"""MPACT - Multimetric PAirwise Comparison Tool version {version} - 24 jan 2025
(c) 2025. Igor Custódio dos Santos & Arthur Gruber

***Logfile (Date: {formatted_date})***

Main parameters:

- Configuration file: {log_conf_file}

- Input fasta file (-i): {fasta_name}
  - Number of sequences: {count_fasta}
  - Type of sequences: {log_type_sequence}
  - Type of headers: {log_type_header}
- Input PDB directory (-s): {pdb_name}
  - Number of files: {count_pdb}
- Method (-p): {methods}
- Output directory (-o): {args.o}

Additional parameters:

- Heatmaps clustering method (-k): {log_clust_method}
- Heatmaps color (-c): {args.c}
- Names in heatmaps (-n): {log_names}
- Heatmap lines (-l): {args.l}
- Mafft parameters (-f): {log_mafft}
- IQ-TREE parameters (-q): {log_iqtree}
- Subgroup parameter (Data partitioning) (-g): {log_data_partitioning}
- Sample subgroups prefixes (-t): {log_prefixes}

""")

def where_to_save_graphics(output):
  """Define what number of graphic directory to save"""
  if "/" in output:
    output_log = output.split("/")[-1]
  else:
    output_log = output
  if os.path.isdir(output)==False:
    os.mkdir(output)
  i = 1
  print(f"Output directory alrealdy exists: {output}")
  while True:
    if os.path.isdir(f"{output}/graphics_dir_{str(i)}") or os.path.isfile(f"{output}/{output_log}_{str(i)}.log"):
      i +=1
      continue
    else:
      save_graphics_dir = f"{output}/graphics_dir_{str(i)}"
      logfile = f"{output}/{output_log}_{str(i)}.log"
      break

  return save_graphics_dir, logfile

def entry_file(args, type_header):
  """Define the input output"""
  if "/" in args.o:
    last_name_input = args.o.split("/")[-1]
  else:
    last_name_input = args.o
  if type_header == 'not_ncbi':
    with open(args.i, 'r') as file:
      with open(f"{args.o}/renamed_{last_name_input}", 'w') as file2:
        for line in file:
          if line.startswith(">"):
            line = line.split(" ")[0]+"\n"
            line = line.replace('(', '')
            line = line.replace(')', '')
            line = line.replace(']', '')
            line = line.replace('[', '')
            line = line.replace('.', '')
            line = line.replace(',', '')
            line = line.replace(':', '')
            line = line.replace(' ', "_")
          file2.write(line)
    enter = f"{args.o}/renamed_{last_name_input}"

  elif type_header == "ncbi":
    enter = args.i

  return enter

def correct_label(data, args, type_header, type_file, order):
  """Correct the label to -n parameter option"""

  if type_header == "not_ncbi" and args.n in ["org", "all"]:
    label = "all"
  elif type_file == "pdb" and args.n in ["org", "all"]:
    label = "all"
  else:
    label = args.n

  if label == 'id':
    new_order = order
    return data, new_order

  elif (label == "org" or label == "all") and type_file == "fasta":
    new_data = []

    with open(args.i, "r") as read_fasta:
      for line in read_fasta:
        if line.startswith(">"):
          taxon = line.split("[")[1].rsplit("]")[0]
          line = line.strip()[1:]
          line = line.replace('(', '')
          line = line.replace(')', '')
          line = line.replace(']', '')
          line = line.replace('[', '')
          line = line.replace(',', '')
          line = line.replace(':', '')
          line = line.replace(' ', "_")
          new_data.append([taxon, line])

    for o in data:
      for t in o[:-1]:
        for n in new_data:
          if t in n[1]:
            if label == "org":
              data[data.index(o)][o.index(t)] = n[0]
            elif label == "all":
              data[data.index(o)][o.index(t)] = n[1]


    new_order = []
    for o in order:
      for seq in new_data:
        if o in seq[1]:
          if label == "org":
            new_order.append(seq[0])
          elif label == "all":
            new_order.append(seq[1])
          continue

    return data, new_order

  elif label == "all" and type_file == "pdb":
    new_data = []

    pdb_list = os.listdir(args.s)
    for pdb in pdb_list:
      pdb = pdb.replace('(', '')
      pdb = pdb.replace(')', '')
      pdb = pdb.replace(']', '')
      pdb = pdb.replace('[', '')
      pdb = pdb.replace(',', '')
      pdb = pdb.replace(':', '')
      pdb = pdb.replace(' ', "_")
      new_data.append(pdb[:-4])

    for o in data:
      for t in o[:-1]:
        for n in new_data:
          if t in n:
            data[data.index(o)][o.index(t)] = n

    new_order = []
    for o in order:
      for seq in new_data:
        if o in seq:
          new_order.append(seq)
        continue

    return data, new_order

def run_needle(fasta_file, output_dir, type_header, type_sequence, log):
  """Run needle program"""
  save_needle_dir = f"{output_dir}/needle_dir"
  if "/" in output_dir:
    output_dir = output_dir.split("/")[-1]
  log.write(f"- Saving in directory: {output_dir}/needle_dir\n")
  if os.path.isdir(save_needle_dir) == False:
    os.mkdir(save_needle_dir)
  sequencias = []
  sequencia_atual = []

  with open(fasta_file, 'r') as arquivo:
    for linha in arquivo:
      linha = linha.strip()
      if linha.startswith(">"):
        if sequencia_atual:
          sequencias.append(sequencia_atual)
        sequencia_atual = [linha]
      else:
        sequencia_atual.append(linha)

    if sequencia_atual:
      sequencias.append(sequencia_atual)

  needle_file = save_needle_dir+"/"+output_dir+".needle"
  stdout_needle_file = save_needle_dir+"/"+output_dir+"_stdout.txt"
  with open(needle_file, 'w') as open_needle:
    open_needle.write("")
  with open(stdout_needle_file, 'w') as open_needle:
    open_needle.write("")

  if type_sequence == 'protein':
    needle_modelcmd_aa = f"needle -sprotein1 sequence.fasta -sprotein2 restofsequences.fasta -gapopen 10.0 -gapextend 0.5 -datafile EBLOSUM62 -outfile temp_out_file.needle >stdout_needle_file"
    log.write(f"- Command model of Needle for protein sequences: {needle_modelcmd_aa}\n")
  if type_sequence == 'nucleotide':
    needle_modelcmd_nucleotide = f"needle -snucleotide1 sequence.fasta -snucleotide2 restofsequences.fasta -gapopen 10.0 -gapextend 0.5 -datafile EDNAFULL -outfile temp_out_file.needle >stdout_needle_file"
    log.write(f"- Command model of Needle for nucleotide sequences: {needle_modelcmd_nucleotide}\n")

  for seq in sequencias:
    in_seq = seq[0][1:].split(' ')[0]
    in_seq_file = seq[0][1:].split(' ')[0]+'.fasta'
    out_seq_file = seq[0][1:].split(' ')[0]+'.needle'
    with open(in_seq_file, 'w') as arquivo:
      arquivo.write('\n'.join(seq))
    with open('comparated_sequences.fasta', 'w') as arquivo:
      for l in sequencias[sequencias.index(seq):]:
        arquivo.write('\n'.join(l)+'\n')
    try:
      if type_sequence == 'protein':
        needle_cmd_aa = f"needle -sprotein1 {in_seq_file} -sprotein2 comparated_sequences.fasta -gapopen 10.0 -gapextend 0.5 -datafile EBLOSUM62 -outfile {out_seq_file} >{stdout_needle_file}"
        subprocess.call(needle_cmd_aa, shell = True)

      if type_sequence == 'nucleotide':
        needle_cmd_nucleotide = f"needle -snucleotide1 {in_seq_file} -snucleotide2 comparated_sequences.fasta -gapopen 10.0 -gapextend 0.5 -datafile EDNAFULL -outfile {out_seq_file} >{stdout_needle_file}"
        subprocess.call(needle_cmd_nucleotide, shell = True)

    except Exception as erro_needle:
      log.write(f"\nError at running Needle:\n{erro_needle}")
      print(f"Error at running Needle:\n{erro_needle}")
      sys.exit
    os.remove(in_seq_file)
    os.remove('comparated_sequences.fasta')
    with open(out_seq_file, 'r') as open_moment_needle:
      with open(needle_file, 'a') as open_needle:
        for linha in open_moment_needle:
          open_needle.write(linha)
    os.remove(out_seq_file)

def data_from_needle(file, type_sequence, metric):
  """Extract the Identity and Similarity values from needle output file and returns a list of comparisons"""
  data_ident = []
  data_sim = []
  with open(file, 'r') as open_needle:
    pair_sim = []
    pair_ident = []
    for needle_line in open_needle:
      if needle_line.startswith("# 1:"):
        pair_sim.append(needle_line[5:-1])
        pair_ident.append(needle_line[5:-1])
      elif needle_line.startswith("# 2:"):
        pair_sim.append(needle_line[5:-1])
        pair_ident.append(needle_line[5:-1])
      elif needle_line.startswith("# Identity:"):
        pair_ident.append(round(float(needle_line.split("(")[1].split("%")[0]), 2))
      elif needle_line.startswith("# Similarity:"):
        pair_sim.append(round(float(needle_line.split("(")[1].split("%")[0]), 2))
      else:
        continue
      if len(pair_sim) == 3:
        data_sim.append(pair_sim)
        pair_sim = []
      if len(pair_ident) == 3:
        data_ident.append(pair_ident)
        pair_ident = []

    if metric == "ident":
      return data_ident
    elif metric == "simil":
      return data_sim

def save_matrix(perc_list, save_path, log):
  """Save a csv format matrix of all-against-all comparisons"""
  matrix_name = save_path

  if save_path.endswith("_ident_matrix.csv"):
    matrix_title_log = "Identity"
  elif save_path.endswith("_simil_matrix.csv"):
    matrix_title_log = "Similarity"
  elif save_path.endswith('_tmscores_matrix.csv'):
    matrix_title_log = "TM-scores"
  elif save_path.endswith('_mldist_matrix.csv'):
    matrix_title_log = "Maximum-likelihood Distance"
  elif save_path.endswith('_3Di_simil_matrix.csv'):
    matrix_title_log = "3Di Characters Similarity"

  order = []
  for l in perc_list:
    if l[0] not in order:
      order.append(l[0])
  matriz = []
  for k in range(len(order)):
    linha = [None] * len(order)
    matriz.append(linha)

  for element in perc_list:
    i = order.index(element[0])
    j = order.index(element[1])
    if i < j:
      h = i
      i = j
      j = h
    matriz[j][i] = matriz[i][j] = element[2]

  table = {}
  table[''] = order
  for l in range(0, len(matriz)):
    table[order[l]] = matriz[l]
  final_table = pd.DataFrame(table)
  final_table.set_index('', inplace=True)
  final_csv = final_table.to_csv(index = True)
  with open(matrix_name, "w") as file:
    file.write(final_csv)

  #Calculate the matrix statistics.
  matrix_inline = []
  matrix_for_stats = matriz
  for l in range(len(matrix_for_stats)):
    matrix_for_stats[l].pop(l)
    matrix_inline += matrix_for_stats[l]

  mean_log = np.mean(matrix_inline)
  std_log = np.std(matrix_inline)
  median_log = np.median(matrix_inline)
  mode_log = stats.mode(matrix_inline)
  range_log = max(matrix_inline) - min(matrix_inline)

  #Write the results on the logfile.
  log.write(f"""
{matrix_title_log} matrix in csv format ({matrix_name}):
{final_csv}

Report of {matrix_title_log} matrix (not considering comparisons of samples against themselves):
- Mean: {mean_log}
- Standard Deviation (STD): {std_log}
- Median: {median_log}
- Mode: {mode_log}
- Range (Max - Min): {range_log}
""")

  return round(mean_log, 4), round(std_log, 4)

def read_matrix_csv(matrix_name):
  """Read csv matrix and return the list of comparisons"""
  with open(matrix_name, "r") as open_matrix:
    extracted_matrix = []
    for line in open_matrix.readlines():
      line = line.strip()
      extracted_matrix.append(line.split(","))

    data = []
    columns = extracted_matrix[0]

    for j, liney in enumerate(extracted_matrix):
      if j == 0:
        continue
      for i, coordx in enumerate(liney):
        if i == 0:
          continue
        data.append([columns[i], columns[j], float(coordx)])

  return data

def nj_tree(perc_list, type_metric, args, log):
  """Generate a neighbor-joining tree from a distance matrix"""
  if "/" in args.o:
    last_name_dir = args.o.split("/")[-1]
  else:
    last_name_dir = args.o
  #Name and create neighbor-joining directory if it doesn't exists.
  save_nj_dir = f"{args.o}/neighbor-joining_dir"
  if os.path.isdir(save_nj_dir) == False:
    os.mkdir(save_nj_dir)

  #Create the name of the files where tree and distance matrix will be save.
  name_njtree = "{}/{}_{}.tree".format(save_nj_dir, last_name_dir, type_metric)
  name_njdistmat = "{}/{}_{}_njdistmat.csv".format(save_nj_dir, last_name_dir, type_metric)

  #Create a list of all sequences names.
  order = []
  longest = ""
  for l in perc_list:
    if l[0] not in order:
      order.append(l[0])
    if len(l[0]) > len(longest):
      longest = l[0]

  #Create empty matrix.
  matriz = []
  for k in range(len(order)):
    linha = [None] * len(order)
    matriz.append(linha)

  #Put every value in its respective place.
  for element in perc_list:
    i = order.index(element[0])
    j = order.index(element[1])
    if i < j:
      h = i
      i = j
      j = h
    if type_metric == "simil" or type_metric == "ident" or type_metric == "3Di":
      matriz[j][i] = matriz[i][j] = round(1-(element[2]/100), 4)
    elif type_metric == "mldist" or type_metric == "2.5Di":
      matriz[j][i] = matriz[i][j] = round(element[2], 4)
    elif type_metric == "tmscores":
      matriz[j][i] = matriz[i][j] = round(1-element[2], 4)

  for n in matriz:
    for o in n:
      if o == None:
        print(order[matriz.index(n)]+" X "+order[n.index(o)])

  #Transform distance matrix in csv and save it in neighbor-joining directory:
  table = {}
  table[''] = order
  for l in range(0, len(matriz)):
    table[order[l]] = matriz[l]
  final_table = pd.DataFrame(table)
  final_table.set_index('', inplace=True)
  final_csv = final_table.to_csv(index = True)
  with open(name_njdistmat, "w") as file:
    file.write(final_csv)

  for_matrix = []
  for l in range(0, len(matriz)):
    for_matrix.append(matriz[l][0:l+1])

  distance_matrix = DistanceMatrix(order, for_matrix)
  constructor = DistanceTreeConstructor()
  njtree = constructor.nj(distance_matrix)
  njtree.root_at_midpoint()
  njtree.ladderize(reverse=True)

  Phylo.write(njtree, name_njtree, "newick")
  name_njtree = "{}/{}_{}.tree".format(save_nj_dir, last_name_dir, type_metric)
  read_tree=open(name_njtree, "r").read()
  l = read_tree.split(",")
  nj_order=[]
  for x in range(len(order)):
    l[x]=l[x].replace(")","")
    l[x]=l[x].replace("(","")
    l[x]=l[x].replace(" ","")
    l[x]=l[x].replace("\t","")
    l[x]=l[x].replace("\n","")
    l[x]=l[x].replace("'","")
    seq_name = l[x].split(":")[0]

    nj_order.append(seq_name)

  return nj_order


def make_mltree_order(perc_list, treefile, args, name_ordered_mltree):
  """Sort and root the IQ-TREE phylogenetic tree"""
  order = []
  longest = ""
  for l in perc_list:
    if l[0] not in order:
      order.append(l[0])
    if len(l[0]) > len(longest):
      longest = l[0]

  ml_tree = Phylo.read(treefile, "newick")
  ml_tree.root_at_midpoint()
  ml_tree.ladderize(reverse=True)

  Phylo.write(ml_tree, name_ordered_mltree, "newick")

  read_tree=open(name_ordered_mltree, "r").read()
  l = read_tree.split(",")
  mltree_order=[]
  for x in range(len(order)):
    l[x]=l[x].replace(")","")
    l[x]=l[x].replace("(","")
    l[x]=l[x].replace(" ","")
    l[x]=l[x].replace("\t","")
    l[x]=l[x].replace("\n","")
    l[x]=l[x].replace("'","")
    seq_name = l[x].split(":")[0]

    mltree_order.append(seq_name)

  return mltree_order

def run_mafft(fasta_file, output, cmd_mafft, log):
  """Run multiple sequence alignment with MAFFT program"""
  if "/" in args.o:
    last_name_dir = output.split("/")[-1]
  else:
    last_name_dir = output

  save_align_dir = f"{output}/mafft_dir"
  if "/" in output:
    output = output.split("/")[-1]
  log.write(f"- Saving in directory: {save_align_dir}\n")
  if os.path.isdir(save_align_dir) == False:
    os.mkdir(save_align_dir)
  real_path_fasta = os.path.realpath(fasta_file)
  real_path_output = os.path.realpath(save_align_dir+"/"+last_name_dir+".align")
  error_out = os.path.realpath(save_align_dir+"/error_"+last_name_dir+"_mafft")

  cmd = 'mafft '+cmd_mafft+' '+real_path_fasta+' > '+real_path_output+' 2> '+error_out
  log.write(f"- MAFFT command: {cmd}\n")
  subprocess.call(cmd, shell=True)

def run_iqtree(align_file, output, cmd_iqtree, align, log):
  """Run phylogenetic reconstruction with IQ-TREE program"""
  if "/" in args.o:
    last_name_dir = output.split("/")[-1]
  else:
    last_name_dir = output

  if align == "3di":
    save_dir = "3diiqtree_dir"
    method = "_3Di"
  elif align == "aa":
    save_dir = "iqtree_dir"
    method = ""
  elif align == "mega":
    save_dir = "megaiqtree_dir"
    method = "_mega"

  save_iqtree_dir = f"{output}/{save_dir}"
  log.write(f"- Saving in directory: {save_iqtree_dir}\n")
  if os.path.isdir(save_iqtree_dir) == False:
    os.mkdir(save_iqtree_dir)

  copyfile(align_file, f"{save_iqtree_dir}/{last_name_dir}{method}")


  real_path_align = os.path.realpath(f"{save_iqtree_dir}/{last_name_dir}{method}")
  error_out = os.path.realpath(f"{save_iqtree_dir}/error1_{last_name_dir}_iqtree")
  error_out2 = os.path.realpath(f"{save_iqtree_dir}/error2_{last_name_dir}_iqtree")

  cmd = f"iqtree2 -s {real_path_align} {cmd_iqtree} 1>{error_out} 2>{error_out2}"

  log.write(f"- IQ-TREE command: {cmd}\n")
  try:
    subprocess.call(cmd, shell=True)
  except Exception as err:
    print("IQTREE reported an error:\n\n"+str(err))
  else:
    print("Running IQ-TREE...")

  if os.path.isfile(f"{save_iqtree_dir}/{last_name_dir}{method}.mldist") == False:
    print("Parameter -m of IQTREE is not valid. Please correct your sintaxe.")
    sys.exit()

def data_from_mldist(mldist_file, log):
  """Extract the maximum likelihood values from IQ-TREE results"""
  data = []
  with open(mldist_file, 'r') as arquivo:
    lines = arquivo.readlines()
    for i in range(1, len(lines)):
      line = lines[i]
      for j in range(1, len(line.split())):
        distance = line.split()[j]
        data.append([line.split()[0], lines[j].split()[0], round(float(distance), 4)])
  return data

def run_tmalign(input_dir, output, log):
  """Run pairwise structural alignment with TM-align program"""
  if "/" in args.o:
    last_name_dir = output.split("/")[-1]
  else:
    last_name_dir = output

  output_tmalign = f"{output}/TMalign_dir"
  log.write(f"- Saving in directory: {output_tmalign}\n")
  if os.path.isdir(output_tmalign) == False:
    os.mkdir(output_tmalign)

  file_list = os.listdir(input_dir)

  for pdb in file_list:
    if pdb.endswith(".pdb")==False:
      print(f'In "{input_dir}" directory, the file "{pdb}" might not be a pdb format file.\nPlease take it off from the input directory or rename it with pdb extention.')
      sys.exit()

  combinations = []
  for yfile in file_list:
    for xfile in file_list:
      if xfile==yfile:
        continue
      if [yfile, xfile] in combinations or [xfile, yfile] in combinations:
        continue
      combinations.append([yfile, xfile])

  log.write("- TM-align command model: TMalign pdb_file1 pdb_file2 -o directory_output -a T 1>tmalign.log\n")

  for pair in combinations:
    in_filex = pair[0]
    in_filey = pair[1]

    out1 = pair[0][:-4].split(' ')[0]
    out2 = pair[1][:-4].split(' ')[0]

    if ":" in pair[0]:
      out1 = out1.replace(":", "")
    if " " in pair[0]:
      out1 = out1.replace(" ", "_")
    if ":" in pair[1]:
      out2 = out2.replace(":", "")
    if " " in pair[1]:
      out2 = out2.replace(" ", "_")
    if os.path.isdir(f"{output_tmalign}/{out1}_x_{out2}")==False and os.path.isdir(f"{output_tmalign}/{out2}_x_{out1}")==False:
      os.mkdir(f"{output_tmalign}/{out1}_x_{out2}")
    elif os.path.isdir(f"{output_tmalign}/{out1}_x_{out2}") and os.path.isfile(f"{output_tmalign}/{out1}_x_{out2}/{out1}_x_{out2}.log"):
      continue
    elif os.path.isdir(f"{output_tmalign}/{out2}_x_{out1}") and os.path.isfile(f"{output_tmalign}/{out2}_x_{out1}/{out2}_x_{out1}.log"):
      continue

    tmalign_cmd = f"TMalign {input_dir}/'{in_filex}' {input_dir}/'{in_filey}' "
    tmalign_cmd += f"-o {output_tmalign}/{out1}_x_{out2}/{out1}_x_{out2} "
    tmalign_cmd += f"-a T 1>{output_tmalign}/{out1}_x_{out2}/{out1}_x_{out2}.log"

    print(f"Running TM-align with {pair[0]} vs {pair[1]}.")

    try:
      subprocess.call(tmalign_cmd, shell=True)
    except:
      print(f"Error at running TM-align with {pair[0]} vs {pair[1]}.")
      sys.exit()

    for tmalign_dirs in os.listdir(f"{output}/TMalign_dir"):
      if tmalign_dirs.endswith("_all"):
        try:
          os.remove(os.path.realpath(f"{output}/TMalign_dir/{tmalign_dirs}"))
        except:
          pass
        try:
          os.remove(os.path.realpath(f"{output}/TMalign_dir/{tmalign_dirs[:-4]}"))
        except:
          pass

    if len(pair[0]) >= 90:
      os.rename(f"{input_dir}/'{in_filex}'", f"{input_dir}/'{pair[0]}'")
    if len(pair[1]) >= 90:
      os.rename(f"{input_dir}/'{in_filey}'", f"{input_dir}/'{pair[1]}'")

def data_from_tmalign(output):
  """Extract TM-score values from TM-align program directory results and return a list of comparisons"""
  if "/" in args.o:
    last_name_dir = output.split("/")[-1]
  else:
    last_name_dir = output

  data = []
  all = []
  for tmalign_dirs in os.listdir(f"{output}/TMalign_dir"):
    if "_x_" not in tmalign_dirs:
      continue
    if tmalign_dirs.split("_x_")[0] not in all:
      all.append(tmalign_dirs.split("_x_")[0])
    if tmalign_dirs.split("_x_")[1] not in all:
      all.append(tmalign_dirs.split("_x_")[1])
    with open(f"{output}/TMalign_dir/{tmalign_dirs}/{tmalign_dirs}.log", 'r') as open_moment_tmalign:
      for line in open_moment_tmalign:
        if "TM-score=" not in line:
          continue
        if "Aligned length=" in line:
          data.append([tmalign_dirs.split("_x_")[0], tmalign_dirs.split("_x_")[1], float(line.split("TM-score=")[1].split(", ")[0])])
        elif "if normalized by average length of two structures" in line:
          data.append([tmalign_dirs.split("_x_")[0], tmalign_dirs.split("_x_")[1], float(line.split("score= ")[1].split(" (if")[0])])
        else:
          continue
  for same in all:
    if [same, same, 1.0] not in data:
      data.append([same, same, 1.0])
  return data


def run_foldseek(args, log):
  """Run 3Di alphabet prediction with Foldseek program and organize the output sequences"""
  if "/" in args.o:
    last_name_dir = args.o.split("/")[-1]
  else:
    last_name_dir = args.o

  if "/" in args.o:
    entry_last_name_dir = args.s.split("/")[-1]
  else:
    entry_last_name_dir = args.s

  if args.s != None and args.s.endswith("/"):
    args.s = args.s[:-1]
  if os.path.isdir(f"{args.o}/foldseek_dir") == False:
    os.mkdir(f"{args.o}/foldseek_dir")
  log.write(f"- Saving in: {args.o}/foldseek_dir\n")
  if args.s != None and args.s.endswith("/"):
    entry_pdbs = f'{args.s}*'
  else:
    entry_pdbs = f'{args.s}/*'
  num_threads = str(round(cpu_count()/2))
  temp_out = f'{last_name_dir}_foldseek_output'
  erro1 = f'{args.o}/foldseek_dir/error1'
  erro2 = f'{args.o}/foldseek_dir/error2'
  try:
    foldseek_cmd = f'foldseek structureto3didescriptor -v 0 --threads {num_threads} --chain-name-mode 1 {entry_pdbs} {temp_out} >{erro1} 2>{erro2}'
    subprocess.call(foldseek_cmd, shell = True)
    log.write(f'- Foldseek command: {foldseek_cmd}\n')
  except Exception as err:
    print("Foldseek reported an error:\n\n"+str(err))
    log.write("\nFoldseek reported an error:\n\n"+str(err)+"\n")
    sys.exit

  from_foldseek_out = temp_out
  to_foldseek_out = f'{args.o}/foldseek_dir/{temp_out}.txt'

  try:
    os.rename(from_foldseek_out, to_foldseek_out)
  except:
    print(f"There was an error in the Foldseek running. Please check out:\n{erro1}\n{erro2}")
    sys.exit()

  from_foldseekdb_out = f'{temp_out}.dbtype'
  to_foldseekdb_out = f'{args.o}/foldseek_dir/{temp_out}.dbtype'

  try:
    os.rename(from_foldseekdb_out, to_foldseekdb_out)
  except:
    print(f"There was an error in the Foldseek running. Please check out:\n{erro1}\n{erro2}")
    sys.exit()

  entry = f'{args.o}/foldseek_dir/{temp_out}.txt'
  file_3Di = f"{args.o}/{entry_last_name_dir}_3Di.fasta"

  log.write(f"""- FASTA file with 3Di characters: {file_3Di}
""")

  with open(entry, "r") as foldseek_file:
    lines = foldseek_file.readlines()

  with open(file_3Di, "w") as open_3Di:
    for line in lines:
      if ".pdb" in line.split('\t')[0]:
        seq_name = line.split('\t')[0].split('.pdb')[0]
      else:
        seq_name = line.split('\t')[0]
      seq_AA = line.split("\t")[1]
      seq_3Di = line.split("\t")[2]
      open_3Di.write(f">{seq_name}\n{seq_3Di}\n")

  os.remove(f"{args.o}/foldseek_dir/{last_name_dir}_foldseek_output.dbtype")


def run_needle3Di(fasta_3Difile, output_dir):
  """Run 3Di alphabet sequences with needle program"""
  if "/" in args.o:
    last_name_dir = output_dir.split("/")[-1]
  else:
    last_name_dir = output_dir

  save_needle3Di_dir = f"{output_dir}/needle3Di_dir"
  log.write(f"- Saving in: {save_needle3Di_dir}\n")

  if os.path.isdir(save_needle3Di_dir) == False:
    os.mkdir(save_needle3Di_dir)

  seqs = []
  actual_seq = []
  with open(fasta_3Difile, 'r') as file_3Di:
    for line in file_3Di:
      line = line.strip()
      if line.startswith(">"):
        if actual_seq:
          seqs.append(actual_seq)
        actual_seq = [line]
      else:
        actual_seq.append(line)

    if actual_seq:
      seqs.append(actual_seq)

  needle_file = save_needle3Di_dir+"/"+last_name_dir+"_3Di_simil.needle"
  stdout_needle_file = save_needle3Di_dir+"/"+last_name_dir+"_stdout.txt"
  with open(stdout_needle_file, 'w') as open_needle:
    open_needle.write("")
  with open(needle_file, 'w') as open_needle:
    open_needle.write("")

  log.write("- Command model of Needle (3Di characters alignment): ")
  log.write(f"needle -asequence sequence.fasta -bsequence restofsequences.fasta -gapopen 8.0 -gapextend 2.0 -datafile 3Di_matrix.txt -aformat3 pair -outfile temp_out_file.needle\n")

  for seq in seqs:
    in_seq = seq[0][1:].split(' ')[0]
    in_seq_file = seq[0][1:].split(' ')[0]+'.fasta'
    out_seq_file = seq[0][1:].split(' ')[0]+'.needle'
    with open(in_seq_file, 'w') as file_3Di:
      file_3Di.write('\n'.join(seq))
    with open('comparated_sequences.fasta', 'w') as file_3Di:
      for l in seqs[seqs.index(seq):]:
        file_3Di.write('\n'.join(l)+'\n')
    try:
      subprocess.call(f"needle -asequence {in_seq_file} -bsequence comparated_sequences.fasta -gapopen 8.0 -gapextend 2.0 -datafile /home/geninfo/isantos/programação/3A-DGT/3Di_matrix.txt -aformat3 pair -outfile {out_seq_file} >{stdout_needle_file}", shell = True)
    except:
      print("Error at running Needle.")
      sys.exit
    else:
      print(f"Running Needle with {in_seq}...")
      time.sleep(1)
    os.remove(in_seq_file)
    os.remove('comparated_sequences.fasta')
    with open(out_seq_file, 'r') as open_moment_needle:
      with open(needle_file, 'a') as open_needle:
        for line in open_moment_needle:
          open_needle.write(line)
    os.remove(out_seq_file)

def data_from_needle3Di(file, type_sequence):
  """Extract 3Di character similarity values and return a list of comparisons"""
  data_score = []
  with open(file, 'r') as open_needle:
    pair_score = []
    for needle_line in open_needle:
      if needle_line.startswith("# 1:"):
        pair_score.append(needle_line[5:-1])
      elif needle_line.startswith("# 2:"):
        pair_score.append(needle_line[5:-1])
      if type_sequence == "3di":
        if needle_line.startswith("# Similarity:"):
          pair_score.append(round(float(needle_line.split("(")[1].split("%")[0]), 2))
      else:
        continue
      if len(pair_score) == 3:
        data_score.append(pair_score)
        pair_score = []
  return data_score


def score_to_dist25di(data):
  dist_data = []
  for AB in data:
    for XX in data:
      if [AB[0], AB[0]] == XX[:2]:
        AA = XX[2]
      if [AB[1], AB[1]] == XX[:2]:
        BB = XX[2]

    dist_data.append([AB[0], AB[1], 1.0-((2.0*AB[2])/(AA+BB))])

  return dist_data

def create_clustermap(all_data, args, type_header, log):
  """Create heatmep with clustermap function and save the dendrogram of the method of cluster."""
  for data in all_data:
    print("Generating clustermap file: "+data+"_"+args.k)
    log.write("- Clustermap file: "+data+"_"+args.k+"\n")

    pdf_file_name = data + '_'+args.k+'.pdf'
    jpg_file_name = data + '_'+args.k+'.jpg'
    svg_file_name = data + '_'+args.k+'.svg'
    tree_file_name = data + '_dend_'+args.k+'.tree'

    metric = ""
    if jpg_file_name.endswith(f"_mldist_{args.k}.jpg"):
      metric = "Maximum Likelihood Distance"
      tick_locs = [0, 3, 6, 9]
      tick_labels = ['0', '3', '6', '9']
      type_file = "fasta"

    elif jpg_file_name.endswith(f"_simil_{args.k}.jpg"):
      metric = "Similarity (%)"
      tick_locs = [0, 25, 50, 75, 100]
      tick_labels = ['0', '25', '50', '75', '100']
      type_file = "fasta"

    elif jpg_file_name.endswith(f"_ident_{args.k}.jpg"):
      metric = "Identity (%)"
      tick_locs = [0, 25, 50, 75, 100]
      tick_labels = ['0', '25', '50', '75', '100']
      type_file = "fasta"

    elif jpg_file_name.endswith(f"_3Di_{args.k}.jpg"):
      metric = "3Di Similarity (%)"
      tick_locs = [0, 25, 50, 75, 100]
      tick_labels = ['0', '25', '50', '75', '100']
      type_file = "pdb"

    elif jpg_file_name.endswith(f"_tmscores_{args.k}.jpg"):
      metric = "TM-score"
      tick_locs = [0, 0.25, 0.5, 0.75, 1]
      tick_labels = ['0', '0.25', '0.5', '0.75', '1']
      type_file = "pdb"

    order = []
    for l in all_data[data]:
      if l[0] not in order:
        order.append(l[0])

    all_data[data], order = correct_label(all_data[data], args, type_header, type_file, order)

    matriz = []
    for k in range(len(order)):
        linha = [None] * len(order)
        matriz.append(linha)

    for element in all_data[data]:
        i = order.index(element[0])
        j = order.index(element[1])
        if i < j:
            h = i
            i = j
            j = h
        matriz[j][i] = matriz[i][j] = element[2]

    table = {}
    table[''] = order
    for l in range(0, len(matriz)):
      table[order[l]] = matriz[l]

    final_table = pd.DataFrame(table)
    final_table.set_index('', inplace=True)

    color = args.c
    if sns.color_palette(color):
      color_msa = color
      if color.endswith("_r"):
        color_pa = color[:-2]
      else:
        color_pa = color+"_r"
    else:
      print("Error: Invalid color (-c).")
      sys.exit()

    if args.l == "yes":
      linew = 0.15
    elif args.l == "no":
      linew = 0

    font_size_cbar = 21
    font_size= 35 / np.sqrt(len(order))


    if jpg_file_name.endswith(f"_mldist_{args.k}.jpg"):
      clustermap = sns.clustermap(final_table, method=args.k,cmap=color_msa,
                                  xticklabels=True, yticklabels=True,
                                  linewidths=linew,
                                  vmin=0, vmax=9,
                                  figsize=(25, 25),
                                  annot_kws={"size": font_size},
                                  cbar_kws={"ticks":tick_locs},
                                  cbar_pos=(.07, .39, .05, .36)
                                  )


    elif jpg_file_name.endswith(f"_simil_{args.k}.jpg") or jpg_file_name.endswith(f"_ident_{args.k}.jpg") or jpg_file_name.endswith(f"_3Di_{args.k}.jpg"):
      clustermap = sns.clustermap(final_table, method = args.k, cmap=color_pa,
                                  xticklabels=True, yticklabels=True,
                                  linewidths=linew,
                                  vmin=0, vmax=100,
                                  figsize=(25, 25),
                                  annot_kws={"size": font_size},
                                  cbar_kws={"ticks":tick_locs},
                                  cbar_pos=(.07, .39, .05, .36)
                                  )

    elif jpg_file_name.endswith(f"_tmscores_{args.k}.jpg"):
      clustermap = sns.clustermap(final_table, method = args.k, cmap=color_pa,
                                  xticklabels=True, yticklabels=True,
                                  linewidths=linew,
                                  vmin=0, vmax=1,
                                  figsize=(25, 25),
                                  annot_kws={"size": font_size},
                                  cbar_kws={"ticks":tick_locs},
                                  cbar_pos=(.07, .39, .05, .36)
                                  )

    cbar = clustermap.cax
    cbar.tick_params(labelsize=font_size_cbar)
    cbar.set_ylabel(metric, rotation=90, va='center', fontsize=font_size_cbar, labelpad=24)
    clustermap.ax_row_dendrogram.set_visible(False)

    plt.savefig(jpg_file_name, dpi=400, bbox_inches="tight")
    plt.savefig(svg_file_name, dpi=400, format='svg', bbox_inches="tight")


    #Save the cluster dendrogram as tree file.
    labels = [t.get_text() for t in clustermap.ax_heatmap.yaxis.get_majorticklabels()]
    link = clustermap.dendrogram_col.linkage

    dict_tree = {numero: 0 for numero in range(len(order)*2-1)}

    for code in order:
      dict_tree[order.index(code)] = code

    node = len(order)
    for coord in link:
      label1 = dict_tree[int(coord[0])]
      label2 = dict_tree[int(coord[1])]

      if int(coord[0]) >= len(order):
        tree_dist_parts = dict_tree[int(coord[0])].rsplit(',', 1)[1].split(':')
        total_branch_distance = 0
        for part in tree_dist_parts:
          if ')' in part:
            total_branch_distance += float(part.split(')')[0])
        dist1 = str(round(coord[2]/2, 4)-total_branch_distance)
      else:
        dist1 = str(round(coord[2]/2, 4))

      if int(coord[1]) >= len(order):
        tree_dist_parts = dict_tree[int(coord[1])].rsplit(',', 1)[1].split(':')
        total_branch_distance = 0
        for part in tree_dist_parts:
          if ')' in part:
            total_branch_distance += float(part.split(')')[0])
        dist2 = str(round(coord[2]/2, 4)-total_branch_distance)
      else:
        dist2 = str(round(coord[2]/2, 4))

      dict_tree[node] = f'({label1}:{dist1},{label2}:{dist2})'
      node += 1

    with open(tree_file_name, 'w') as write_tree:
      write_tree.write(dict_tree[len(order)*2-2])

    plt.clf()

    return labels


def create_njheatmap(all_data, nj_order, args, order_type, type_header, log):
  """Create heatmap with heatmap function"""
  for data in all_data:
    print("Generating Neighbor-Joining heatmap file: "+data+"_"+order_type)
    log.write("- Neighbor-joining heatmap file: "+data+"_"+order_type+"\n")

    jpg_file_name = data + '_'+order_type+'.jpg'
    svg_file_name = data + '_'+order_type+'.svg'

    metric = ""
    if jpg_file_name.endswith("_mldist_"+order_type+".jpg"):
      metric = "Maximum Likelihood Distance"
      tick_locs = [0, 3, 6, 9]
      tick_labels = ['0', '3', '6', '9']
      type_file = "fasta"

    elif jpg_file_name.endswith("_simil_"+order_type+".jpg"):
      metric = "Similarity (%)"
      tick_locs = [0, 25, 50, 75, 100]
      tick_labels = ['0', '25', '50', '75', '100']
      type_file = "fasta"

    elif jpg_file_name.endswith("_ident_"+order_type+".jpg"):
      metric = "Identity (%)"
      tick_locs = [0, 25, 50, 75, 100]
      tick_labels = ['0', '25', '50', '75', '100']
      type_file = "fasta"

    elif jpg_file_name.endswith("_3Di_"+order_type+".jpg"):
      metric = "3Di Similarity (%)"
      tick_locs = [0, 25, 50, 75, 100]
      tick_labels = ['0', '25', '50', '75', '100']
      type_file = "pdb"

    elif jpg_file_name.endswith("_tmscores_"+order_type+".jpg"):
      metric = "TM-score"
      tick_locs = [0, 0.25, 0.5, 0.75, 1]
      tick_labels = ['0', '0.25', '0.5', '0.75', '1']
      type_file = "pdb"

    all_data[data], nj_order = correct_label(all_data[data], args, type_header, type_file, nj_order)

    matriz = []
    for k in range(len(nj_order)):
        linha = [None] * len(nj_order)
        matriz.append(linha)

    for element in all_data[data]:
        i = nj_order.index(element[0])
        j = nj_order.index(element[1])
        if i < j:
            h = i
            i = j
            j = h
        matriz[j][i] = matriz[i][j] = element[2]

    table = {}
    table[''] = nj_order
    for l in range(0, len(matriz)):
      table[nj_order[l]] = matriz[l]

    final_table = pd.DataFrame(table)
    final_table.set_index('', inplace=True)

    color = args.c
    if sns.color_palette(color):
      color_msa = color
      if color.endswith("_r"):
        color_pa = color[:-2]
      else:
        color_pa = color+"_r"
    else:
      print("Error: Invalid color (-c).")
      sys.exit()

    if args.l == "yes":
      linew = 0.15
    elif args.l == "no":
      linew = 0

#    plt.subplots(figsize=(25, 25))

    plt.rcParams['svg.fonttype'] = 'none'

    sns.set (rc = {'figure.figsize':(25, 25)})

    font_size = 25 / np.sqrt(len(nj_order))
    font_size_cbar = 50 / np.sqrt(len(nj_order))




    if jpg_file_name.endswith("_mldist_"+order_type+".jpg"):
      heatmap = sns.heatmap(final_table, cmap=color_msa,
                   xticklabels=True, yticklabels = True, linewidths=linew,
                   vmin=0, vmax=9,
                   cbar=False,
                   annot_kws={"size": font_size})

    elif jpg_file_name.endswith("_simil_"+order_type+".jpg") or jpg_file_name.endswith("_ident_"+order_type+".jpg") or jpg_file_name.endswith("_3Di_"+order_type+".jpg"):
      heatmap = sns.heatmap(final_table, square=True, cmap=color_pa,
                   xticklabels=True, yticklabels = True,
                   linewidths=linew,
                   vmin=0, vmax=100,
                   cbar=False,
                   annot_kws={"size": font_size})

    elif jpg_file_name.endswith("_tmscores_"+order_type+".jpg"):
      heatmap = sns.heatmap(final_table, cmap=color_pa,
                   xticklabels=True, yticklabels = True,
                   linewidths=linew,
                   vmin=0, vmax=1,
                   cbar=False,
                   annot_kws={"size": font_size})

    heatmap.set_aspect("equal")
    heatmap.tick_params(left=False, right=True, bottom=True, labelleft=False, labelright=True, labelbottom=True, rotation=0)
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90, fontsize=font_size)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=font_size)
    plt.xlabel('', fontsize=font_size)

    heatmap.tick_params(bottom=True, top=False, axis='x', length = 3, width = .2)
    heatmap.tick_params(right=True, left=False, axis='y', length = 3, width = .2)

    plt.rcParams['axes.grid'] = False
    cbar = plt.colorbar(heatmap.collections[0], ax=heatmap, location='left', shrink=0.4, aspect=6, extend='neither')
    cbar.ax.yaxis.set_ticks_position('right')
    cbar.ax.yaxis.set_label_position('right')
    cbar.ax.set_ylabel(metric, rotation=90, va='center', fontsize=font_size_cbar)
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.set_yticklabels(tick_labels, rotation=0, fontsize=font_size_cbar)
    cbar.ax.tick_params(right=True, left=False, axis='y', length = 5, width = .2)
    cbar.outline.set_linewidth(0)


    plt.tight_layout()
    plt.savefig(jpg_file_name, dpi=400, bbox_inches="tight")
    plt.savefig(svg_file_name, dpi=400, format='svg', bbox_inches="tight")
    plt.clf()

def cluster_in_subdata(clust_cmd_list, subdata, order, args, align_type, log):
  """Execute the data partitioning, according to user-defined thresolds"""

  print("Data partitioning in execution...")
  method_partit_dict = {"0": "Identity",
                 "1": "Similarity",
                 "2": "Maximum likelihood distance",
                 "3": "TM-score",
                 "4": "3Di similarity",
                 }

  method_partitioning = clust_cmd_list[0]
  lower_limit = float(clust_cmd_list[1])
  upper_limit = float(clust_cmd_list[2])
  if method_partitioning in ["0", "1", "4"]:
    log.write(f"- Partitioning argument (-g): {args.g} ({lower_limit}% to {upper_limit}% of {method_partit_dict[method_partitioning]})\n")
    print(f"- Partitioning argument (-g): {args.g} ({lower_limit}% to {upper_limit}% of {method_partit_dict[method_partitioning]})")
  elif method_partitioning in ["2", "3"]:
    log.write(f"- Partitioning argument (-g): {args.g} ({lower_limit} to {upper_limit} of {method_partit_dict[method_partitioning]})\n")
    print(f"- Partitioning argument (-g): {args.g} ({lower_limit} to {upper_limit} of {method_partit_dict[method_partitioning]})")

  log.write("\nOrder of comparison: {';'.join(order)}\n")

  new_subdata = []
  for comp in subdata:
    if comp[0] == comp[1]:
      continue
    if [comp[1], comp[0], comp[2]] in new_subdata:
      continue
    new_subdata.append(comp)
  subdata = new_subdata

  outgroups = {}
  already_in = []
  sequences_in = []
  cluster_report = []

  i=1
  for query in order:
    if query in sequences_in:
      continue
    log.write(f"Cluster {i} (HEAD: {query}):\n")
    l = 1
    subgroup = [query]
    for comp in subdata:
      if query in comp:
        if comp not in already_in:
          already_in.append(comp)
          if comp[2] > lower_limit and comp[2] < upper_limit:
            if comp[0] == query and comp[1] not in subgroup:
              log.write(f"- {comp[2]} with {comp[1]} (In range... ACCEPTED IN CLUSTER {i})\n")
              subgroup.append(comp[1])
              sequences_in.append(comp[1])
              l+=1
            else:
              log.write(f"- {comp[2]} with {comp[0]} (In range... ACCEPTED IN CLUSTER {i})\n")
              subgroup.append(comp[0])
              sequences_in.append(comp[0])
              l+=1
          else:
            if comp[0] == query and comp[1] not in subgroup:
              log.write(f"- {comp[2]} with {comp[1]} (Out of range...)\n")
            else:
              log.write(f"- {comp[2]} with {comp[0]} (Out of range...)\n")
    for seq in subgroup:
      for comp in subdata:
        if seq in comp and comp not in already_in:
          already_in.append(comp)
    outgroups["{}cl{}s".format(i, l)] = subgroup
    print(f"- Cluster {i} ({i}cl{l}s): {','.join(subgroup)}.")
    cluster_report.append(f"  * Cluster {i} ({i}cl{l}s): {','.join(subgroup)}.\n")
    log.write("\n")
    i+=1


  method_for_name = {"0": "ident", "1": "simil", "2": "mldist", "3": "tm", "4": "3di", "5": "3dimldist", "6": "megamldist"}
  method_name = method_for_name[clust_cmd_list[0]]

  dict_clstr = args.o+"/clstr_"+method_name+"_"+str(clust_cmd_list[1])+"_"+str(clust_cmd_list[2])
  if os.path.isdir(dict_clstr) == False:
    os.mkdir(dict_clstr)

  log.write(f"\nData partitioning report (Directory: {dict_clstr}):\n")
  print("\n")
  for line_report in cluster_report:
    log.write(line_report)

  if clust_cmd_list[0] == "0" or clust_cmd_list[0] == "1" or clust_cmd_list[0] == "2":
    sequences = []
    actual_sequence = []
    with open(args.i, 'r') as fasta_file:
      lines = fasta_file.readlines()
      l=1
      for line in lines:
        linep = line.strip()
        if linep.startswith(">"):
          if actual_sequence:
            sequences.append(actual_sequence)
          actual_sequence = [linep]
        else:
          actual_sequence.append(linep)
        if l == len(lines):
          sequences.append(actual_sequence)
        l+=1

    for name in outgroups:
      with open(dict_clstr+"/"+name+".fasta", "w") as write_fasta:
        for code in outgroups[name]:
          for sequence in sequences:
            if code in sequence[0]:
              for line in sequence:
                write_fasta.write(line+'\n')

  elif clust_cmd_list[0] == "3" or clust_cmd_list[0] == "4":

    aa3to1={
       'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
       'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
       'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
       'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
       'MSE':'M',
    }

    pdb_list = os.listdir(args.s)

    pdbs = []
    for pdb in pdb_list:
      with open(args.s+"/"+pdb, "r") as read_pdb:
        lines = read_pdb.readlines()

      current_seq = ""
      current_pos = 1
      for line in lines:
        if line.startswith("ATOM"):
          list_attributes = no_empty(line.split(" "))
          if isInt(list_attributes[4]):
            if list_attributes[4] == str(current_pos):
              current_seq += aa3to1[list_attributes[3]]
              current_pos += 1
          else:
            if list_attributes[5] == str(current_pos):
              current_seq += aa3to1[list_attributes[3]]
              current_pos += 1
      pdbs.append(">{}\n{}\n".format(pdb[:-4], current_seq))

    for name in outgroups:
      with open(dict_clstr+"/"+name+".fasta", "w") as write_fasta:
        for code in outgroups[name]:
          for pdb in pdbs:
            if code in pdb:
              write_fasta.write(pdb)


def create_freqplot(all_data):
  """Create the frequency distribuition plot."""
  for data in all_data:
    print(f"Generating frequency plot file: {data}_freqplot")
    log.write("- Frequency plot file: "+data+"_freqplot\n")

    jpg_file_name = data + '_freqplot.jpg'
    svg_file_name = data + '_freqplot.svg'

    if data.endswith("_ident") or data.endswith("_simil") or data.endswith("_3Di"):
      freqplot_data = {numero: 0 for numero in range(101)}
      for n in all_data[data]:
        if n[0] != n[1]:
          if round(n[2]) in freqplot_data:
            freqplot_data[round(n[2])] += 1

    elif data.endswith("_tmscores") or data.endswith("_2.5Di"):
      freqplot_data = {numero/100: 0 for numero in range(101)}
      for n in all_data[data]:
        if n[0] != n[1]:
          if round(n[2], 2) in freqplot_data:
            freqplot_data[round(n[2], 2)] += 1

    elif data.endswith("_mldist") or data.endswith("_3dimldist") or data.endswith("_megamldist"):
      freqplot_data = {numero: 0 for numero in [round(x * 0.1, 1) for x in range(91)]}
      for n in all_data[data]:
        if n[0] != n[1]:
          if round(n[2], 1) in freqplot_data:
            freqplot_data[round(n[2], 1)] += 0.5

    freqplot_data = pd.DataFrame.from_dict(freqplot_data, orient = 'index')
    freqplot_data = freqplot_data.rename(columns={0: ""})

    sns.set_style("white")
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['axes.linewidth'] = 1.2
    plt.subplots(figsize=(6.4, 4.8))

    freqplot = sns.lineplot(data=freqplot_data, palette='Greens_d')

    for spine in freqplot.spines.values():
      spine.set_color('black')

    freqplot.tick_params(axis='both', colors='black')

    if data.endswith("_ident"):
      freqplot.set_xlabel('Pairwise identity (%)')
      freqplot.set_title('Identity distribution')
    elif data.endswith("_simil"):
      freqplot.set_xlabel('Pairwise similarity (%)')
      freqplot.set_title('Similarity distribution')
    elif data.endswith("_mldist"):
      freqplot.set_xlabel('Maximum-likelihood distance')
      freqplot.set_title('Maximum-likelihood distance distribution')
    elif data.endswith('_tmscores'):
      freqplot.set_xlabel('TM-score')
      freqplot.set_title('TM-score distribution')
    elif data.endswith('_3Di'):
      freqplot.set_xlabel('3Di similarity (%)')
      freqplot.set_title('3Di similarity distribution')
    freqplot.set_ylabel('Absolute frequency')

    def num2freq(x):
      return 100*(x/(len(all_data[data])/2))
    def freq2num(x):
      return (x/100)*(len(all_data[data])/2)
    secax = freqplot.secondary_yaxis('right', functions=(num2freq, freq2num))
    secax.set_ylabel('Relative frequency (%)')

    plt.savefig(jpg_file_name, dpi=200)
    plt.savefig(svg_file_name, dpi=200, format='svg')

    plt.clf()

######################################### BLOCK OF FUNCTIONS (TABLE OF RANGES) ######################################

def create_range_table(csv_file, prefixes, mean, std):

  """Read the comparison matrix from a CSV file and return index order and values."""
  if csv_file.endswith("mldist"):
    with open(csv_file, 'r') as file:
      size = int(file.readline().strip())

      data = file.readlines()

    index_order = []
    values_matrix = []

    for line in data:
      parts = line.split()
      index_order.append(parts[0])
      values_matrix.append([float(value) for value in parts[1:]])

  else:
    df = pd.read_csv(csv_file, index_col=0)
    index_order = df.index.tolist()
    values_matrix = df.values.tolist()


  """Get the min and max comparison values for each prefix."""
  results = {}
  index_to_position = {index: pos for pos, index in enumerate(index_order)}

  for prefix in prefixes:
    intragroup_min = 100
    intragroup_max = -100
    intergroup_min = 100
    intergroup_max = -100
    for i, sample1 in enumerate(index_order):
      for j, sample2 in enumerate(index_order):
        if i!=j:
          if sample1.startswith(prefix) and sample2.startswith(prefix):
            if values_matrix[i][j] > intragroup_max:
              intragroup_max = values_matrix[i][j]
            if values_matrix[i][j] < intragroup_min:
              intragroup_min = values_matrix[i][j]
          elif sample1.startswith(prefix) and not sample2.startswith(prefix):
            if values_matrix[i][j] > intergroup_max:
              intergroup_max = values_matrix[i][j]
            if values_matrix[i][j] < intergroup_min:
              intergroup_min = values_matrix[i][j]
          else:
            continue

    results[prefix] = {
            'intragroup': {'min': intragroup_min, 'max': intragroup_max},
            'intergroup': {'min': intergroup_min, 'max': intergroup_max}
    }

  results["Mean+-STD"] = f"{mean}+-{std}"

  return results

def save_results_to_csv(results, prefixes, output_file, log):
  """Save the results dictionary to a CSV file."""
  range_table = {}

  first_column = []
  for prefix in prefixes:
    first_column += [prefix, "intragroup","intergroup"]
  first_column += ["Mean+-STD"]
  range_table["Prefix"] = first_column

  for metric in results:
    column = []
    for prefix in results[metric]:
      if prefix == "Mean+-STD":
        column.append(results[metric][prefix])
        continue
      column.append("")
      for type_comparison in results[metric][prefix]:
        maximum  = results[metric][prefix][type_comparison]['max']
        minimum  = results[metric][prefix][type_comparison]['min']
        column.append(f"{minimum}-{maximum}")
    range_table[metric]=column

  results_df = pd.DataFrame(range_table)
  ranges_csv = results_df.to_csv(index = False)
  results_df.to_csv(output_file, sep=';', encoding='utf-8', index=False)

  log.write(f"""

Range table based on tags/prefixes:

{ranges_csv}""")

################################################################################################################################


if __name__ == '__main__':
  #Print help message if no parameter was specified or if it's called (-h)
  if not len(sys.argv)>1:
    print(help)
  elif args.help == True:
    print(help)

  #Print version message if it's called (-v)
  elif args.version == True:
    print("""
MPACT - Multimetric PAirwise Comparison Tool version """+version+""" - 24 jan 2025
This program perform all-against-all pairwise comparisons for biological sequences 
and protein structures with five diferent metrics. The results are presented in 
heatmap graph, frequency distribution, trees and matrices.
(c) 2024. Igor Custódio dos Santos & Arthur Gruber
""")
  else:
    #Set parameters according to configuration file if it's specified (-conf)
    if args.conf != None:
      args.i,args.s,args.p,args.o,args.c,args.l,args.n,args.g,args.f,args.q,args.t,args.k = configuration_file(args)

    #Check the parameters individually.
    type_header, type_sequence = check_fasta(args.i)
    check_pdb(args.s)
    args.p = check_methods(args)
    check_heatmap_param(args)
    cmd_iqtree = check_iqtree(args.q)
    cmd_mafft = check_mafft(args.f)
    check_clustering_method(args.k)

    args.o = check_output_dir(args.o)
    if "/" in args.o:
      last_name_dir = args.o.split("/")[-1]
    else:
      last_name_dir = args.o

    if args.g:
      clust_cmd_list = check_clustering(args.g)
    else:
      clust_cmd_list = None

    if args.t:
      prefixes = check_prefixes(args)
    else:
      prefixes = None

  #  else:
   #   type_header = 'not_ncbi'
    #  type_sequence = "PDB"

    #Defines the paths to all main output directorys.
    save_graphics_dir, logfile = where_to_save_graphics(args.o)
    save_needle_dir = f"{args.o}/needle_dir"
    save_align_dir = f"{args.o}/mafft_dir"
    save_iqtree_dir = f"{args.o}/iqtree_dir"
    save_tmalign_dir = f"{args.o}/TMalign_dir"
    save_foldmason_dir = f"{args.o}/foldmason_dir"
    save_3Dineedle_dir = f"{args.o}/needle3Di_dir"


#    if args.i and os.path.isdir(args.i) == False:
 #     args.i = entry_file(args, type_header)
    with open(logfile, "w") as log:
      log_header(log, args, type_sequence, type_header, version)

      all_data = {}

      data_for_ranges_table = {}

      part1 = "no"
      if "0" in args.p or "1" in args.p or "2" in args.p:
        part1 = "yes"
      if args.i == None and part1 == "yes":
        print("""
WARNING:
Not possible to run pairwise aa/nucleotide sequence alignment ('0' and '1') 
or maximum likelihood ('2').
FASTA file (-i) not specified.
""")
        log.write("""
WARNING:
Not possible to run pairwise aa/nucleotide sequence alignment ('0' and '1') 
or maximum likelihood ('2').
FASTA file (-i) not specified.
""")
      else:

        if "0" in args.p or "1" in args.p:
          log.write("""
***** Pairwise Alignment *****
""")
          print("\nMethod: pairwise sequence alignment with NW")
          if os.path.isfile(save_needle_dir+"/"+last_name_dir+".needle"):
            log.write("""
- Needle file already exists. Skipping Needle run.
- Using data from: """+save_needle_dir+"/"+args.o+".needle\n")
            print("Needle file already exists. Skipping Needle run.")
          else:
            log.write("\nNeedle\n")
            print("Running pairwise sequence alignments with needle…")
            run_needle(args.i, args.o, type_header, type_sequence, log)
            print("Pairwise sequence alignments completed.")
                                                                                                


          #Identity protocol execution
          if "0" in args.p:
            log.write("""
********** Identity Protocol **********
""")
            print("\nIdentity protocol")
            all_data = {}
            matrix_name = save_needle_dir+"/"+last_name_dir+"_ident_matrix.csv"
            if os.path.isfile(matrix_name):
              data = read_matrix_csv(matrix_name)
            else:
              print(save_needle_dir+"/"+last_name_dir+".needle", type_sequence, 'ident')
              data = data_from_needle(save_needle_dir+"/"+last_name_dir+".needle", type_sequence, 'ident')
#            print(data)
            mean, std = save_matrix(data, matrix_name, log)
            if prefixes != None:
              data_for_ranges_table["AA identity %"] = create_range_table(matrix_name, prefixes, mean, std)

            all_data[save_graphics_dir+"/"+last_name_dir+"_ident"] = data

            log.write("\nIdentity outputs:\n")

            nj_order = nj_tree(data, "ident", args, log)
            clstr_order = nj_order

            if clust_cmd_list and clust_cmd_list[0] == "0":
              log.write("\nData partitioning based on Identity values:\n")
              
              for subdata in all_data:
                cluster_in_subdata(clust_cmd_list, all_data[subdata], clstr_order, args, "ident", log)

            else:
              if os.path.isdir(save_graphics_dir) == False:
                os.mkdir(save_graphics_dir)

              print("Generating frequency distribution plot of pairwise sequence alignment distance…")
              create_freqplot(all_data)
              print("Done.")

              print("Generating heatmap plot of pairwise sequence alignment distance…")
              create_njheatmap(all_data, nj_order, args, "nj", type_header, log)
              if args.k != None and args.k != "nj":
                create_clustermap(all_data, args, type_header, log)
              print("Done.")


            all_data = {}


          #Similarity protocol execution
          if "1" in args.p:
            log.write("""
********** Similarity Protocol **********
""")
            print("\nSimilarity protocol")
            all_data = {}
            matrix_name = save_needle_dir+"/"+last_name_dir+"_simil_matrix.csv"
            if os.path.isfile(matrix_name):
              data = read_matrix_csv(matrix_name)
            else:
              data = data_from_needle(save_needle_dir+"/"+last_name_dir+".needle", type_sequence, 'simil')
            mean, std = save_matrix(data, matrix_name, log)
            if prefixes != None:
              data_for_ranges_table["AA similarity %"] = create_range_table(matrix_name, prefixes, mean, std)

            all_data[save_graphics_dir+"/"+last_name_dir+"_simil"] = data

            log.write("\nSimilarity outputs:\n")

            nj_order = nj_tree(data, "simil", args, log)
            clstr_order = nj_order

            if clust_cmd_list:
              for subdata in all_data:
                cluster_in_subdata(clust_cmd_list, all_data[subdata], clstr_order, args, "simil", log)

            else:
              if os.path.isdir(save_graphics_dir) == False:
                os.mkdir(save_graphics_dir)

              print("Generating frequency distribution plot of pairwise sequence alignment distance…")
              create_freqplot(all_data)
              print("Done.")

              print("Generating heatmap plot of pairwise sequence alignment distance…")
              create_njheatmap(all_data, nj_order, args, "nj", type_header, log)
              if args.k != None and args.k != "nj":
                create_clustermap(all_data, args, type_header, log)
              print("Done.")

            all_data = {}



	#Maximum-likelihood protocol execution
        if "2" in args.p:
          log.write("""
********** Maximum Likelihood Distance **********
""")
          print("\nMethod: Maximum-likelihood distance")

          log.write("\nMAFFT:\n")

          if os.path.isfile(save_align_dir+"/"+last_name_dir+".align"):
            log.write("""- Multiple sequences alignment file already exists. Skipping MAFFT run.
- Using data from: {save_align_dir}/{last_name_dir}.align\n""")
            print("MAFFT aligned file already exists. Skipping MAFFT step.")

          else:
            print("Running multiple sequence alignment with MAFFT…")
            run_mafft(args.i, args.o, cmd_mafft, log)
            print("Multiple sequence alignment completed.")



          log.write("\nMAFFT:\n")
          if os.path.isfile(save_iqtree_dir+"/"+last_name_dir+".mldist"):
            log.write(f"""- IQ-TREE maximum-likelihood distance matrix file already exists. Skipping IQ-TREE run.
- Using data from: {save_iqtree_dir}/{last_name_dir}.mldist\n""")
            print("IQ-TREE maximum-likelihood distance matrix file already exists. Skipping IQ-TREE step.")

          else:
            print("Running phylogenetic analysis IQ-TREE…")
            run_iqtree(save_align_dir+"/"+last_name_dir+".align", args.o, cmd_iqtree, "aa", log)
            print("Phylogenetic analysis completed.")

          if prefixes != None:
            data_for_ranges_table["ML distance"] = create_range_table(save_iqtree_dir+"/"+last_name_dir+".mldist", prefixes, mean, std)
          data = data_from_mldist(save_iqtree_dir+"/"+last_name_dir+".mldist", log)

          matrix_name = f"{save_iqtree_dir}/{last_name_dir}_mldist_matrix.csv"
          mean, std = save_matrix(data, matrix_name, log)

          all_data[save_graphics_dir+"/"+last_name_dir+"_mldist"] = data

          log.write("\nMaximum-likelihood distance outputs:\n")

          nj_order = nj_tree(data, "mldist", args, log)

          name_mltree_ordered = save_iqtree_dir+"/"+last_name_dir+"_ordered.treefile"
          mltree_order = make_mltree_order(data, save_iqtree_dir+"/"+last_name_dir+".treefile", args, save_iqtree_dir+"/"+last_name_dir+"_ordered.treefile")
          clstr_order = mltree_order

          if clust_cmd_list:
            for subdata in all_data:
              cluster_in_subdata(clust_cmd_list, all_data[subdata], clstr_order, args, "simil", log)

          else:
            if os.path.isdir(save_graphics_dir) == False:
              os.mkdir(save_graphics_dir)

            print("Generating frequency distribution plot of pairwise sequence alignment distance…")
            create_freqplot(all_data)
            print("Done.")

            print("Generating heatmap plot of pairwise sequence alignment distance…")
            create_njheatmap(all_data, nj_order, args, "nj", type_header, log)
            create_njheatmap(all_data, mltree_order, args, "phylogeny", type_header, log)
            if args.k != None and args.k != "nj":
              create_clustermap(all_data, args, type_header, log)
            print("Done.")

          all_data = {}


      part2 = "no"
      if "3" in args.p or "4" in args.p:
        part2 = "yes"
      if args.s == None and part2 == "yes":
        print("""
WARNING:
Not possible to run pairwise structural alignment ('3'),
or pairwise 3di character alignment ('4').
PDB files directory (-s) not specified.
""")
        log.write("""
WARNING:
Not possible to run pairwise structural alignment ('3'),
or pairwise 3di character alignment ('4').
PDB files directory (-s) not specified.
""")
      else:
        if args.s != None and args.s.endswith("/"):
          args.s = args.s[:-1]

        if "3" in args.p:
          print("\nMethod: Structural alignment (TM-score)")
          log.write("""
********** Structural alignment (TM-score) **********
""")
          log.write("\nTM-align:\n")
          print("Running pairwise structural alignments with TM-align")
          run_tmalign(args.s, args.o, log)
          print("Pairwise structural alignments completed.")

          matrix_name = save_tmalign_dir+"/"+last_name_dir+"_tmscores_matrix.csv"
          if os.path.isfile(matrix_name):
            data = read_matrix_csv(matrix_name)
          else:
            data = data_from_tmalign(args.o)
          save_matrix(data, matrix_name, log)
          if prefixes != None:
            data_for_ranges_table["TM-score"] = create_range_table(matrix_name, prefixes, mean, std)

          mean, std = save_matrix(data, matrix_name, log)
          all_data[save_graphics_dir+"/"+last_name_dir+"_tmscores"] = data

          log.write("\nTM-score outputs:\n")

          nj_order = nj_tree(data,  "tmscores", args, log)
          clstr_order = nj_order

          if clust_cmd_list:
            for subdata in all_data:
              cluster_in_subdata(clust_cmd_list, all_data[subdata], clstr_order, args, "simil", log)

          else:
            if os.path.isdir(save_graphics_dir) == False:
              os.mkdir(save_graphics_dir)

            print("Generating frequency distribution plot of pairwise sequence alignment distance…")
            create_freqplot(all_data)
            print("Done.")

            print("Generating heatmap plot of pairwise sequence alignment distance…")
            create_njheatmap(all_data, nj_order, args, "nj", type_header, log)
            if args.k != None and args.k != "nj":
              create_clustermap(all_data, args, type_header, log)
            print("Done.")

          all_data = {}


        if "4" in args.p:
          if "/" in args.o:
            entry_last_name_dir = args.s.split("/")[-1]
          else:
            entry_last_name_dir = args.s

          print("\nMethod: pairwise 3Di-character alignment")
          log.write("""
********** 3Di Similarity **********
""")
          file_AA_3Di = f"{args.o}/{entry_last_name_dir}_AA_3Di.fasta"
          file_3Di = f"{args.o}/{entry_last_name_dir}_3Di.fasta"

          log.write("\nFoldseek:\n")
          if os.path.isfile(file_AA_3Di) == False and os.path.isfile(file_3Di) == False:
            print("Converting PDB to 3Di-character sequences with Foldseek")
            log.write("\nRunning Foldseek:\n")
            run_foldseek(args, log)
            print("Conversion completed.\n")
          else:
            log.write(f"""- Skipping Foldseek run. Output files already exists:
  - {file_3Di}
  - {file_AA_3Di}\n""")
            print("3Di and aa/3Di FASTA files already exist. Skipping Foldseek run.")


          needle_file = f'{args.o}/needle3Di_dir/{last_name_dir}_3Di_simil.needle'
          log.write("\nNeedle (3Di characters alignment):\n")
          if os.path.isfile(needle_file) == False:
            print("Running pairwise 3Di-character sequence alignments with needle…")
            run_needle3Di(file_3Di, args.o)
            print("Pairwise 3Di-character sequence alignments completed.")
          else:
            log.write("- Needle 3Di similarity file already exists. Skipping Needle run.\n")
            print("Needle 3Di similarity file already exists. Skipping Needle run.")

          matrix_name = save_3Dineedle_dir+"/"+last_name_dir+"_3Di_simil_matrix.csv"
          if os.path.isfile(matrix_name):
            data = read_matrix_csv(matrix_name)
          else:
            data = data_from_needle3Di(needle_file, "3di")
          mean, std = save_matrix(data, matrix_name, log)
          if prefixes != None:
            data_for_ranges_table["3Di similarity %"] = create_range_table(matrix_name, prefixes, mean, std)

#          nj_order = nj_tree(data, "3Di", args, log)

          all_data[save_graphics_dir+"/"+last_name_dir+"_3Di"] = data

          log.write("\n3Di Similarity outputs:\n")

          nj_order = nj_tree(data, "3Di", args, log)
          clstr_order = nj_order

          if clust_cmd_list:
            for subdata in all_data:
              cluster_in_subdata(clust_cmd_list, all_data[subdata], clstr_order, args, "simil", log)

          else:
            if os.path.isdir(save_graphics_dir) == False:
              os.mkdir(save_graphics_dir)

            print("Generating frequency distribution plot of pairwise sequence alignment distance…")
            create_freqplot(all_data)
            print("Done.")

            print("Generating heatmap plot of pairwise sequence alignment distance…")
            create_njheatmap(all_data, nj_order, args, "nj", type_header, log)
            if args.k != None and args.k != "nj":
              create_clustermap(all_data, args, type_header, log)
            print("Done.")

          all_data = {}

      if prefixes:
#        print(data_for_ranges_table)
        save_results_to_csv(data_for_ranges_table, prefixes, args.o+"/metrics_ranges_table.csv", log)


