## Installation

1. Clone the repository
2. Create a conda [environment](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)
```
conda create --name precog python=2.7 scipy numpy pandas scikit-learn=0.19.0
```
3. Activate the environment
```
conda activate precog
```
4. Install the [HMMER](http://hmmer.org/download.html) package
>**Tip:** Use version **v3.1b2** for output comparable to that of the [webserver](precog.russelllab.org)/[original paper](https://pubmed.ncbi.nlm.nih.gov/31143927/)

5. See **help** to know the command-line parameters
 ```
./precog.py --help
```

>**Note:** This version does not include structure-based predictions( using *InterPreTS*). It has been tested on Ubuntu-based environments.

## Input
1. The input file must be [**FASTA formatted**](https://en.wikipedia.org/wiki/FASTA_format). For example:
>**\>sp|P30518|V2R_HUMAN Vasopressin V2 receptor OS=Homo sapiens OX=9606 GN=AVPR2 PE=1 SV=1**
MLMASTTSAVPGHPSLPSLPSNSSQERPLDTRDPLLARAELALLSIVFVAVALSNGLVLA
ALARRGRRGHWAPIHVFIGHLCLADLAVALFQVLPQLAWKATDRFRGPDALCRAVKYLQM
VGMYASSYMILAMTLDRHRAICRPMLAYRHGSGAHWNRPVLVAWAFSLLLSLPQLFIFAQ
RNVEGGSGVTDCWACFAEPWGRRTYVTWIALMVFVAPTLGIAACQVLIFREIHASLVPGP
SERPGGRRRGRRTGSPGEGAHVSAAVAKTVRMTLVIVVVYVLCWAPFFLVQLWAAWDPEA
PLEGAPFVLLMLLASLNSCTNPWIYASFSSSVSSELRSLLCCARGRTPPSLGPQDESCTT
ASSSLAKDTSS
>**\>Q14330**
MITRNNQDQPVPFNSSHPDEYKIAALVFYSCIFIIGLFVNITALWVFSCTTKKRTTVTIY
MMNVALVDLIFIMTLPFRMFYYAKDEWPFGEYFCQILGALTVFYPSIALWLLAFISADRY
MAIVQPKYAKELKNTCKAVLACVGVWIMTLTTTTPLLLLYKDPDKDSTPATCLKISDIIY
LKAVNVLNLTRLTFFFLIPLFIMIGCYLVIIHNLLHGRTSKLKPKVKEKSIRIIITLLVQ
VLVCFMPFHICFAFLMLGTGENSYNPWGAFTTFLMNLSTCLDVILYYIVSKQFQARVISV
MLYRNYLRSMRRKSFRSGSLRSLSNINSEML

2. Mutations can be given in the following format:
>\>Q14330/**Y60R**
>MITRNNQDQPVPFNSSHPDEYKIAALVFYSCIFIIGLFVNITALWVFSCTTKKRTTVTI**R**  
MMNVALVDLIFIMTLPFRMFYYAKDEWPFGEYFCQILGALTVFYPSIALWLLAFISADRY  
MAIVQPKYAKELKNTCKAVLACVGVWIMTLTTTTPLLLLYKDPDKDSTPATCLKISDIIY  
LKAVNVLNLTRLTFFFLIPLFIMIGCYLVIIHNLLHGRTSKLKPKVKEKSIRIIITLLVQ  
VLVCFMPFHICFAFLMLGTGENSYNPWGAFTTFLMNLSTCLDVILYYIVSKQFQARVISV  
MLYRNYLRSMRRKSFRSGSLRSLSNINSEML

3. However, the mutated sequence of a GPCR may not be provided again if it's wild type sequence has already been mentioned before.
>\>Q14330
>MITRNNQDQPVPFNSSHPDEYKIAALVFYSCIFIIGLFVNITALWVFSCTTKKRTTVTI**Y**  
MMNVALVDLIFIMTLPFRMFYYAKDEWPFGEYFCQILGALTVFYPSIALWLLAFISADRY  
MAIVQPKYAKELKNTCKAVLACVGVWIMTLTTTTPLLLLYKDPDKDSTPATCLKISDIIY  
LKAVNVLNLTRLTFFFLIPLFIMIGCYLVIIHNLLHGRTSKLKPKVKEKSIRIIITLLVQ  
VLVCFMPFHICFAFLMLGTGENSYNPWGAFTTFLMNLSTCLDVILYYIVSKQFQARVISV  
MLYRNYLRSMRRKSFRSGSLRSLSNINSEML
>**\>Q14330/Y60R**

## Ouput
The output file contains following headers:
1. **GPCR/MUT**:
Name of the input. Must be alphanumeric.
2. **GNAI3 - GNAL**:
Couping values predicted by PRECOG (probabilities) or known from IUPHAR (PC: Primary Coupling, SC: Secondary Coupling) or the shedding assay experiment (LogRAi >= -1.0 is coupled, otherwise uncoupled). If the input sequence cannot be searched against a G-protein model or is unavaliable in IUPHAR or the data from the shedding assay experiment, it is shown with **-**.
3. **7TM1_POS/BW/ALN_POS**:
Denotes the Pfam 7tm1 (7TM1_POS) position, Ballesteros-Weinstein numbering (BW) or alignment position (ALN_POS) of the positions affected in the given input GPCR.
4. **Mutation_Info**:
Information about the input mutation. Not applicable for wild type input.

## Citation
[PubMed](https://pubmed.ncbi.nlm.nih.gov/31143927/)
>Singh G, Inoue A, Gutkind JS, Russell RB, Raimondi F. PRECOG: PREdicting COupling probabilities of G-protein coupled receptors. Nucleic Acids Res. 2019 Jul 2;47(W1):W395-W401. doi: 10.1093/nar/gkz392. PMID: 31143927; PMCID: PMC6602504.

## Contact
gurdeep[at]bioquant[.]uni-heidelberg[.]de
