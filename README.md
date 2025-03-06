# MPACT

MPACT, the Multimetric PAirwise Comparison Tool, is an integrated toolbox for pairwise all-against-all comparison of primary nucleotide or amino acid sequences, and 3D protein structures. 

## Instalation

MPACT does not need to be installed. The user should only download one of the Docker files or the mpact.py script.

## Requirements

### Using the docker file

In this case, you must only install the Docker program (https://www.docker.com/) on your computer.

### Using the MPACT python script

MPACT uses the following third-party programs to perform the sequence and structural comparisons:

* Needle
* MAFFT
* IQ-TREE
* Foldseed
* TM-align

## Usage
### Docker file
In this repository, two docker images are available. The first image (mpact.tar.gz) contains all third-party programs and was generated on a computer with an Intel® Core™ i7-5500U 2.40 GHz processor and a Linux operating system (Ubuntu 24.04.1 LTS). Due to Foldseek limitations, this image is only compatible with machines with an Intel® Core™ i7-5500U 2.40 GHz processor or higher. The second image (mpact_no_fs.tar.gz) does not contain the Foldseek program and therefore does not have the processor limitations.

To use both images, first you need to decompress the file with gzip program:

```
gzip -d mpact.tar.gz
```
Then you must upload the docker image to your local repository:

```
docker load < mpact.tar
```

To run the docker image, you need to set your local path in the run command:

```
docker run -v <your_path>:/mnt mpact -i /mnt/input.fasta -s /mnt/input_dir -o /mnt/output_dir
```
### Python script

```
python3 mpact.py -i input.fasta -s input_dir -o output_dir
```
### Mandatory parameters
-i | input <file name>:      Input fasta file containing biological sequences. 

-s | structure <directory>:  Directory containing 3D protein structure files in PDB format. 

### Optional parameters
-conf <file>:                Configuration file.

-c | color <parameter>:      Color palette for graphs (default = RdBu).

-f | MAFFT <parameters>:     Allows the user to change parameters –maxitarate and –threads to be used by MAFFT.

-l (lines) <yes|no>:         Use cells dividing lines on heatmap (default = no).

-n | names <id|org|all>:     Names for labels.                       
			                          This parameter will work only if the sequence headers have the NCBI standard (e.g.: ">YP_003289293.1 RNA-dependent RNA polymerase 
                                [Drosophila totivirus 1]). Otherwise, sequence headers appended with “_” in whitespace will be used as labels. For the name showed 
                                in the example: 
                                
                                id: YP_003289293.1
									              org: Drosophila totivirus 1
									              all: YP_003289293.1 Drosophila totivirus 1

-o | output <name>:          Output directory.

-p | protocols <als3>:       Protocols to be used in the execution. By default, all protocols will be used.

-q | IQ-TREE <parameters>:   IQTREE parameters. Only "-m", "-bb" e "-nt" are accepted. 
				                        e.g.: -iqtree "-m Q.pfam+F+R6 -bb 1000 -nt 20”

-v | version:                Program version.






-p | protocols <als3>:       Protocols to be used in the execution. By default, all protocols will be used. 

## Tutorial

Follow the instructions in the MPACT_manual.pdf file to learn how to use MPACT program and interpret the results.

## Reference

If you use this program for your publication, please cite:



