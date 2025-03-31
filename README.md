# MPACT

MPACT, the Multimetric PAirwise Comparison Tool, is an integrated toolbox for pairwise all-against-all comparison of primary nucleotide or amino acid sequences, and 3D protein structures. 

## Instalation

We provide docker files for an easy installation. Alternatively, 
the user can download the MPACT program and perform a manual installation 
of the required third-party applications. 

## Requirements

### Using the docker file

In this case, you must install the Docker program (https://www.docker.com/) on your computer.

### Using the MPACT python script

MPACT uses the following third-party programs to perform the sequence and structural comparisons:

* Needle (https://emboss.sourceforge.net/apps/release/6.0/emboss/apps/needle.html)
* MAFFT (https://mafft.cbrc.jp/alignment/software/)
* IQ-TREE (https://github.com/iqtree/)
* Foldseek (https://github.com/steineggerlab/foldseek)
* TM-align (https://zhanggroup.org/TM-align/)

## Usage
### Docker file
Two MPACT docker images are available on the Docker repository. These images come pre-installed with the third-party programs mentioned in the previous section, eliminating the need for manual installation. The first image, named mpact, includes all third-party programs and was created on a system with an Intel® Core™ i7-5500U 2.40 GHz processor running Linux (Ubuntu 24.04.1 LTS). Due to Foldseek limitations, this image is only compatible with machines equipped with an Intel® Core™ i7-5500U 2.40 GHz processor or higher. The second image, mpact_no_fs, does not include the Foldseek program and, therefore, is not subject to the processor restrictions mentioned above.

To use them, the user only needs to install Docker (https://docs.docker.com/) and download one of them. To get the complete image, use the following command line:
s
```
docker push lilianesantana/mpact:tagname
```

To obtain the image without the Foldseek program, use the command below:

```
docker push lilianesantana/mpact_no_fs:tagname
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
```
-i | input <file name>:      Input fasta file containing biological sequences. 

-s | structure <directory>:  Directory containing 3D protein structure files in PDB format.
```

### Optional parameters
```
-conf <file>:                Configuration file.

-c | color <parameter>:      Color palette for graphs (default = RdBu).

-f | MAFFT <parameters>:     Allows the user to change parameters –maxitarate and –threads to be used by MAFFT.

-l (lines) <yes|no>:         Use cells dividing lines on heatmap (default = no).

-n | names <id|org|all>:     Names for labels.                       
			        This parameter will work only if the sequence headers have the NCBI standard (e.g.: ">YP_003289293.1 RNA-dependent RNA polymerase 
                                [Drosophila totivirus 1]). Otherwise, sequence headers appended with “_” in whitespace will be used as labels. For the name shown 
                                in the example: 
                                
                                id: YP_003289293.
				org: Drosophila totivirus 1
				all: YP_003289293.1 Drosophila totivirus 1

-o | output <name>:          Output directory.

-p | protocols <als3>:       Protocols to be used in the execution. By default, all protocols will be used.

-q | IQ-TREE <parameters>:   IQTREE parameters. Only "-m", "-bb" e "-nt" are accepted. 
			        e.g.: -iqtree "-m Q.pfam+F+R6 -bb 1000 -nt 20”

-v | version:                Program version.

```

## Tutorial

Follow the instructions in the MPACT_manual.pdf file to learn how to use MPACT program and interpret the results.

## Reference

If you use this program for your publication, please cite:

MPACT program (developed by Igor Custódio dos Santos and Arthur Gruber, University of São Paulo, Brazil, manuscript in preparation).

## Contact

To report bugs, to ask for help and to give any feedback, please contact Arthur Gruber (argruber@usp.br) or Igor Custódio dos Santos (igor223tgec@usp.br).

