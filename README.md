# CRISPR-Cas-Docker
In silico docking of CRISPR RNA and Cas protein

## Description

The CCD server docks a CRISPR RNA (crRNA) onto its Cas protein in silico and generates Top 10 docking models using HDOCK. This bioinformatic tool will be useful when your prokaryotic genomes have multiple CRISPR arrays and Cas systems, so optimal crRNA-Cas protein pairs are not clear. As a preliminary study, CCD can predict the RNA-protein interaction of a given crRNA-Cas pair through in silico experiments, before conducting any time-consuming and expensive in vitro and in vivo experiments. You can either provide experimental 3D structures of your crRNA and Cas protein directly or use 3D-predicted crRNA and AlphaFold-predicted Cas proteins.

## Getting Started

### When you have no Cas protein and crRNA experimental structures:

#### 1. Have your Cas protein sequence ready

You may have a Cas gene cluster from your sequencing experiment or from a public database like CRISPRCasDB (https://crisprcas.i2bc.paris-saclay.fr/MainDb/StrainList). For now, the CCD can only handle the 3D structure prediction of Class II effector proteins such as Cas9, Cas12, and Cas13 that are more widely used as genomic tools. There are plans to upgrade the CCD to handle multimeric Cas proteins of Class I.

After extracting the DNA sequence of your Cas gene, have it translated into a protein sequence (using a bioinformatic translation tool, such as https://www.bioinformatics.org/sms2/translate.html). Please pay attention to the reading frame and the strand direction. Also, please remove the * at the end of the protein, if present.

#### 2. Have your crRNA sequence ready

Next, extract your crRNA sequence (either in RNA or DNA) from your sequencing experiment or from a public database like CRISPRCasDB (https://crisprcas.i2bc.paris-saclay.fr/MainDb/StrainList).

#### 3. Provide Cas protein sequence as an input to CCD

The next step requires you to provide an amino acid sequence of your Cas protein to predict its 3D structure using AlphaFold. You can either directly input an amino acid sequence or upload a FASTA file containing a single amino acid sequence in the standard format.

#### 4. Provide crRNA sequence as an input to CCD

The final step requires you to provide an RNA or DNA sequence of your crRNA to predict its 3D structure using RNAFold and RoseTTAFold. You can either directly input an RNA or DNA sequence or upload a FASTA file containing a single RNA or DNA sequence in the standard format.

#### 5. Fill in your query information and press [send]

For your Cas protein, it may take 30 minutes to run AlphaFold to predict its 3D structure. For your crRNA, it may take 15 minutes to run RNAFold and RoseTTAFold to predict its 3D structure. For your crRNA and Cas protein pair, it may take 15 minutes to generate the Top 10 docking models from HDOCK docking experiments.

### When you have Cas protein and crRNA experimental structures:

#### 1. Provide Cas protein structure as an input to CCD

The first step requires you to provide the 3D structure of your Cas protein. You can either directly input a PDB ID of your Cas protein from the RCSB Protein Data Bank (https://www.rcsb.org/) or upload a PDB file containing a single protein structure in the standard format.

#### 2. Provide crRNA structure as an input to CCD

The next step requires you to provide the 3D structure of your crRNA in the standard PDB format.

#### 3. Fill in your query information and press [send]

For your crRNA and Cas protein pair, it may take 15 minutes to generate the Top 10 docking models from HDOCK docking experiments.

## More Help

http://www.crisprcasdocker.org/help

CRISPR-Cas-Docker: a web-based in-silico docking of crRNA to Cas protein with machine learning classification (in preparation).

## Authors

* Ho-min Park
* Jongbum Won 
* Yunseol Park 
* Esla Timothy Anzaku
* Joris Vankerschaver
* Arnout Van Messem
* Wesley De Neve
* Hyunjin Shim

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgements

2D RNA structures and 3D RNA structures were predicted with ViennaRNA v.2.5.1 and RoseTTAFold v.2.0.0, respectively. In silico docking experiments were modeled with HDOCK v.1.1.0. Protein structures were predicted with AlphaFold2, available under an open-source license at https://github.com/deepmind/alphafold. For protein structure similarity metrics, we used TM-align (https://zhanggroup.org/TM-align). 3-D Structure visualizations were created with 3Dmol.js (https://3dmol.csb.pitt.edu/doc/tutorial-embeddable.html). For data analysis, Python .3.8.13 (https://www.python.org), NumPy v.1.23.4 (https://github.com/numpy/numpy), seaborn v.0.12.0 (https://github.com/mwaskom/seaborn), Matplotlib v.3.5.3 (https://github.com/matplotlib/matplotlib), pandas v.1.4.3 (https://github.com/pandas-dev/pandas) were used.

## References

H. Park, Y. Park, J. Vankerschaver, A. Van Messem, W. De Neve, H. Shim. Rethinking protein drug design with highly accurate structure prediction of anti-CRISPR proteins. Pharmaceuticals 15:3, 310 (2022).

H. Park*, Y. Park*, U. Berani, E. Bang, J. Vankerschaver, A. Van Messem, W. De Neve, H. Shim. In silico optimization of RNA-protein interactions for CRISPR-Cas13-based antimicrobials. Biology Direct 17:1, 1-16 (2022).
