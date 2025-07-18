# Project study
Project study of my Master of science.

To run this tool please execute the main.py from the folder it is located in.   
Also you need to  download the sql dump of the cripsr cas database:  
https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=ccpp_db.zip  
Save the sql dump in "crispr_cas_db/db_parser/Crispr_Cas_Database_SQL_Dump.sql".

You also need to download RNAmotiFold from here:  
https://github.com/RNABioInfo/RNAmotiFold  
Save it in the same folder the main.py is located in. 

Be aware, that it is possible that the RNAmotiFold algorith might not work due to changes in the tool.

Many of the classes in this project are generic, abstract, or serve as parent classes for others. This design choice was made to promote code reuse. Much of the code is expected to be refactored over time to be more generic, rather than hardcoded. The only part that will remain implementation-specific is the database processing logic, which needs to be rewritten for each distinct database, and a small part of the analysis.  
Additionally, this project currently lacks testing and exception handling. Since it was primarily developed for internal data analysis and not intended as a widely used tool, these aspects were initially deprioritized. However, they may be introduced in a future version of the project. The documentation is also kept to a minimum.


Following is the underlying theory of the project study:

**RNA folding and architecture**

Ribonucleic acid (RNA) plays an important role in regulation processes, catalysis and gene
expression. Since the three-dimensional structure of RNA dictates its function, understanding
RNA folding is essential for deciphering its biological functions and for developing potential
RNA-targeted applications [1]. RNA has multiple levels of structural organization. The primary
structure refers to its nucleotide sequence. The secondary structure is defined as the set of
paired and unpaired bases within RNA. This structure is stabilized by hydrogen bonding and
base-stacking interactions, which result in the formation of stems, loops and bulges. The
tertiary structure represents the three-dimensional conformation of the RNA molecule. It is
formed by the interaction between secondary structure elements, including long-range
contacts [2]. The RNA structure formation is hierarchical: The primary sequence determines
the secondary structure through local base-pairing interactions, and subsequently folds into
the tertiary structure through long-range interactions [1].

RNA 3D motifs are found in hairpin loops, internal loops and multihex-junctions within
structured RNA [3]. RNA 3D motifs consist of recurrent sets of nucleotides that adopt specific
three-dimensional structures stabilized by characteristic non-canonical base pairs [4]. RNA 3D
motifs are conserved across unrelated RNA molecules, reflecting their structural and
functional importance [3], [5]. RNA 3D motifs serve a variety of functions, they can mediate
long-range interactions with other RNAs, within the same RNA, as well as with proteins or
other ligands. Additionally, some RNA 3D motifs play critical roles in shaping RNA architecture.
For example, C-loops increase the helical twist of the RNA helix. Another important function
of RNA 3D motifs is their stabilizing role during RNA folding: UNCG-tetraloops can act as
nucleation sites for the formation of hairpin loop stems [3].

**RNA structure prediction**

Accurate prediction of RNA secondary structure is crucial for understanding RNA, as its
structure dictates its function. Traditional secondary structure prediction methods aim to
identify the thermodynamically most stable structure, known as the minimum free energy
(mfe) structure [6], [7]. However, these methods often fail to predict native RNA structures
accurately, because they exclude modified bases, rely on an imperfect energy model and do
not account for the complex cellular environment or interactions with other molecules that
alter RNA conformation. Focusing solely on the mfe structure can be limiting, as it overlooks
the diversity of alternative structures that may also be functionally relevant. To address this,
classified dynamic programming (cDP) approaches, that are based on algebraic dynamic
programming (ADP), enable the grouping of predicted RNA structures based on shared
structural features [8]. Such a cDP approach is employed in RNAshapes, which utilizes shape
abstraction. This method groups RNA structures based on their composition and
arrangements of loops, resulting in abstract shape classes, each with a representative mfe
structure [9]. A similar strategy is utilized by RNAHeliCes, which uses positional abstraction,
where each helix is represented through its central position and associated loop type [10].

RNAmotiFold employs an alternative cDP approach by incorporating RNA 3D motifs into the
loop regions of secondary structures by matching loop sequences against known motif
sequences. Thus, predicted structures can be grouped based on their motifs, resulting in
classes, each represented by a mfe structure [11]. For this study, the terms “RNA 3D motif”,
“3D motif” and “motif” will be used synonymously, unless otherwise specified.

**CRISPR RNA**

The CRISPR-Cas system is an adaptive immune mechanism that protects prokaryotes from
mobile genetic elements such as phages. CRISPR stands for “Clustered Regularly Interspaced
Palindromic Repeats,” while Cas refers to “CRISPR-associated” proteins [12]. CRISPR-Cas
systems include a CRISPR array, which contains a leader sequence followed by short,
conserved repeat sequences interspersed with unique spacers [13]. Spacers derived from
foreign genetic elements enable CRISPR-Cas systems to recognize and target complementary
nucleic acid sequences [14]. The repeats are necessary for the recognition and binding of
CRISPR RNAs (crRNAs) by Cas proteins [15]. The genes encoding Cas proteins are typically
organized into an operon or gene cluster and are responsible for processing precursor CRISPR
RNA (pre-crRNA) into mature crRNA, which then guides Cas proteins to complementary
sequences in invading nucleic acids, facilitating their sequence-specific cleavage. This renders
the foreign genetic material harmless [14].

CRISPR-Cas systems are classified into multiple types and subtypes based on the identity and
organization of their cas genes. These subtypes differ in how they process crRNA and in the
protein complexes they express, which are responsible for the interference mechanism [12].
Depending on the subtype, the crRNA is either largely unstructured or adopt a specific
secondary structure, which serves as a structural motif that is recognized by Cas proteins [14]. 
Consequently, the repeats and crRNA represent ideal targets for secondary structure and RNA 3D motif analysis.