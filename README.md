# ish-312-ec
# Repository for Research on AP3D1 gene and ARPC gene family 
Methods have been taught by Dr. Rest and should not be posted anywhere else without permission.

# 1. BLAST AP3D1 against the database to identify AP3D1 homologs.

You will be conducting a blast searches for the ARPC protein.

The gene code can be found on websites such as the NIH. You can download it directly from the website using the following command:

    ncbi-acc-download -F fasta -m protein NP_001248755.1
 Now, perform a blast search using the query protein:

     blastp -db ../allprotein.fas -query NP_001248755.1.fa -outfmt 0 -max_hsps 1 -out AP3D1.blastp.typical.out
Type in the following command there after: 

    blastp -db ../allprotein.fas -query NP_001248755.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out AP3D1.blastp.detail.out


This will produce our blastp file to use to find putative homologs 

  
## Filtering the BLAST output for putative homologs

Next we need to get our putative homologs to include for future trees:

For our specific research, we require the e-value to be less than 1e-13.

Here is a command you can use to filter our output file to satisfy this requirement:

    awk '{if ($6<0.0000000000001)print $1 }' arpc.blastp.detail.out > arpc.blastp.detail.filtered.out

This will produce a filtered.out file we can use and have a set of putative homologs to the query that is needed for the research.

    wc -l AP3D1.blastp.detail.filtered.out

# 2. Obtain and align the FMR1 gene family sequences

now that proteomes have been downloaded, we need to use seqkit to produce a .fas file of our homologs. 

     seqkit grep --pattern-file AP3D1.blastp.detail.filtered.out ../allprotein.fas > AP3D1.homologs.fas

## Perform a global multiple sequence alignment in mafft: 
Now we are going to use MAFFT which is a multiple sequence alignment program. We need to install it first before we use it. 

	wwget -O ~/tools/mafft-7.490-linux.tgz https://mafft.cbrc.jp/alignment/software/mafft-7.490-linux.tgz && tar xfvz ~/tools/mafft-7.490-linux.tgz -C ~/tools
```Now, use the following command to make a multiple sequence alignment using mafft:

```
	~/tools/mafft-linux64/mafft.bat --localpair --maxiterate 100  AP3D1.homologs.fas   > AP3D1.homologs.al.fas 
Provide some statistics about the alignment using [t_coffee](https://www.tcoffee.org/Projects/tcoffee/):

     t_coffee -other_pg seq_reformat -in AP3D1.homologs.al.fas -output sim

# . Constructing a Phylogenetic Tree for AP3D1 Homologs from Sequence Data

The data for your tree will be the alignment for ARPC homologs that you made in lab3:  `AP3D1.homologs.al.fas`.

Use IQ-TREE to find the maximum likehood tree estimate, specifying a particular model of substitution (VT+F+R10). 
```bash
iqtree -s AP3D1.homologs.al.fas -m VT+F+R10
```

## Root the optimal phylogeny

Let's root the tree, so that our view of it is more meaningful. By rooting, we will specify the divergence event that is the oldest.
 
### Midpoint rooting
 For now, let's use a type of rooting called midpoint - we'll hope that the root is halfway along the longest branch on the tree. Note: this may not be true, in which case the root will be wrong.
 
    gotree reroot midpoint -i AP3D1.homologs.al.fas.treefile -o AP3D1.homologs.al.mid.treefile

Now, we can look at the rooted tree at the command line:
    
	    nw_order -c n AP3D1.homologs.al.mid.treefile  | nw_display -

Note: the `nw_order` first ladderizes (orders) the clades by numbers of descendents. This can make large trees easier to look at. (But can also contribute to some biases in interpretation.)

Or, we can draw it graphically:
``` 
nw_order -c n AP3D1.homologs.al.mid.treefile  |   nw_display  -w 1000 -b 'opacity:0' -s  >  AP3D1.homologs.al.mid.svg -
```
# 3. Reconciling the FMR1 gene family with the species tree

## The Species Tree
At this point we created the species tree as part of trees to be used for bootstrap and domains ``species.tre``.

## The Gene Tree

We will be using the midpoint rooted gene tree was done previously. 
`answers.AP3D1.homologs.al.fas.treefile` 


`ls ~ myusername/AP3D1/AP3D1.homologs.al.mid.treefile`  

The resulting output should be: `/home/ec2-user myusername/AP3D1/AP3D1.homologs.al.mid.treefile`


## Reconcile the gene and species tree using Notung
Use the following command to perform the reconciliation:

```bash
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ../speciesTreeBilateriaCnidaria.tre -g AP3D1.homologs.al.mid.treefile --reconcile --speciestag prefix --savepng --events
```
Examine the file ``AP3D1.homologs.al.mid.treefile.reconciled.events.txt``. (Hint: use ``less`` with the `-S` command flag, or push this to your remote repository.)

Transfer and then look at the graphic of the reconciled gene tree: ``AP3D1.homologs.al.mid.treefile.reconciled.png``

This is one of the trees we will use! 

### Generate a RecPhyloXML object to view the gene-within-species tree
Use the following command:
```bash
    python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g AP3D1.homologs.al.mid.treefile.reconciled --include.species
```

The resulting file is called ``AP3D1.homologs.al.mid.treefile.reconciled.xml``
 
To create a gene-reconciliation-with species tree reconciliation, you may use thirdkind:
 ```bash
thirdkind -f arpc.homologs.al.mid.treefile.reconciled.xml -o  arpc.homologs.al.mid.treefile.reconciled.svg
 ```
 
Note, thirdkind cuts off the events at the root of the tree. 



## Rooting via Minimization of Duplications and Deletions in Notung



In notung, we can replace the --reconcile command with the --root flag to re-root the tree during using the reconciliation process. The command is: 

```bash
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ../speciesTreeBilateriaCnidaria.tre -g AP3D1.homologs.al.mid.treefile --root --speciestag prefix --savepng --events
```

The output file will be the same as before, but with  "rooting.0" added to the file name. For example: ``AP3D1.homologs.al.mid.treefile.rooting.0.events.txt``

 
# 4. Evaluating node (branch) support for the ARPC gene tree using the bootstrap 

In lecture, we discussed bootstrap support.
Let's re-run IQ-TREE, this time generating bootstrap support values. 


To obtain bootstrap support, its necessary iqtree tree search, and ask it to estimate bootstrap support using the ultrafast bootstrap. This will result in a bootstrapped tree for each bootstrap replicate we request.  

First, copy the alignment and tree files to the lab6 ARPC directory. 
```bash
cp ~ myusername/AP3D1.homologs.al.fas ~myusername/.
cp AP3D1.homologs.al.fas.treefile ~myusername/.
```

Then run this command in order to repeat the bootstrap command 

```bash=1
iqtree -s AP3D1.homologs.al.fas -nt AUTO -bb 1000 -m VT+F+R10 -t AP3D1.homologs.al.fas.treefile -pre AP3D1.bs -wbt
```
 This file with the 1000 bootstrapped trees is called `AP3D1.bs.ufboot`. (But you will need to run the command for your own gene).


Next, we'll use the set of bootstrapped trees to add support values to our optimal tree: (run this command and all following commands)
```bash
iqtree -t AP3D1.bs.ufboot -sup AP3D1.homologs.al.fas.treefile 
```
The resulting tree that contains bootstrap values is called ``AP3D1.bs.ufboot.suptree``.

Now, midpoint root the optimal tree again (the tree will be the same, but now it has bootstrap support values associated with the nodes.):
```bash
gotree reroot midpoint -i AP3D1.bs.ufboot.suptree -o AP3D1.bs.mid.suptree
```

Finally, print out the tree to a graphics file with the bootstrap support:
```bash
nw_display -s  AP3D1.bs.mid.suptree -w 1000 -b 'opacity:0' >  AP3D1.bs.mid.suptree.svg
```
Transfer this graphics file to your remote repository so you can view it.

Recall that you can also display this as a ladderized cladogram (i.e. without branch lengths):
```bash
nw_order -c n AP3D1.bs.mid.suptree | nw_topology - | nw_display -s -w 1000 > AP3D1.bs.mid.suptree.Cl.svg -
```
This will produce our gene trees with bootstrap support!
# Domain Identification for ARPC Proteins 

Now we are going to use RPS-BLAST to obtain the gene tree with domains within the protein sequences.

## Input Protein Sequences
We need to use unaligned protein sequence. We also need the sequence names to match the names we used in our phylogenetic analysis.
We will use `AP3D1.homologs.fas`, which were our original, unaligned sequences from lab 3. 

Make sure you are in the ARPC directory:

```bash    
cd ~ myusername/AP3D1
```
Let's make a copy of our raw unaligned sequence, removing the asterisk (stop codon) in the process. 
```bash
sed 's/*//' ~ myusername/AP3D1/AP3D1.homologs.fas > ~ myusername/AP3D1/AP3D1.homologs.fas
```
## Download the Pfam database, run RPS-BLAST, and reformat

**Note that the following command (line 2 by the green bar) only needs to be run one time on your instance!**
```bash=2
  wget -O ~/data/Pfam_LE.tar.gz ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz && tar xfvz ~/data/Pfam_LE.tar.gz  -C ~/data
```
Now, you can run RPS-BLAST:
```bash
rpsblast -query myusername/AP3D1/AP3D1.homologs.fas -db /home/ec2-user/data/Pfam -out ~myusername/AP3D1/AP3D1.rps-blast.out  -outfmt "6 qseqid qlen qstart qend stitle"
```
>``
Now reformat the output so that [EvolView](https://www.evolgenius.info/evolview-v3/) can use it:
```bash
awk 'BEGIN{FS="\t"} {print $1"\t"$2"\t"$3"@"$4"@"$5}' ~ myusername/AP3D1/AP3D1.rps-blast.out | sed 's/,/_/g' | sed 's/ /_/g' | sed 's/__/_/g' | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > ~ myusername/AP3D1.rps-blast.evol.out
```
## Plot the predicted Pfam domains on the phylogeny.

Instead of useing Evolview, we are using the command that Dr. Rest gave us that is ggtree!

```
awk 'BEGIN{FS="\t"} {gsub(/ |\.|,/, "_", $5)1} {print $1"\t"$2"\t"$3"@"$4"@"$5}' ~myusername/AP3D1/AP3D1.rps-blast.out | awk -F "@" 'BEGIN { OFS=FS }; { print $1,$2, substr($3, 1, 20); }' | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g'   > ~AP3D1/AP3D1.rps-blast.evol.out
```
This should get us our gene tree with domain support!

```
awk 'BEGIN{FS="\t"} {gsub(/ |\.|,/, "_", $5)1} {print $1"\t"$2"\t"$3"@"$4"@"$5}' ~myusername/AP3D1/AP3D1.rps-blast.out | awk -F "@" 'BEGIN { OFS=FS }; { print $1,$2, substr($3, 1, 20); }' | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g'   > ~AP3D1/AP3D1.rps-blast.evol.out
```
This should get us our gene tree with domain support!
