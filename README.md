# Project study
This repository contains the code for my Master's project study.

**The aim of this project was to examine repeat RNA and CRISPR RNA structures and their RNA 3D motifs from the CRISPR-Cas database.**

To run this tool please execute the main.py from the folder it is located in.   
You need to  download the sql dump of the cripsr cas database:  
https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=ccpp_db.zip  
Save the sql dump in "crispr_cas_db/db_parser/Crispr_Cas_Database_SQL_Dump.sql".

Furthermore, you need to download RNAmotiFold from here:  
https://github.com/RNABioInfo/RNAmotiFold  
Save it in the same folder the main.py is located in. 

Be aware, that it is possible that the RNAmotiFold algorith might not work due to changes in the tool. 
In that case try following commit: 33f9e8c131daabd151af74753c2dbbcee5751e32

After running the main.py a plot folder will be created, containing the visualizations of the results.

Many of the classes in this project are generic, abstract, or serve as parent classes for others. This design choice was made to promote code reuse. Much of the code is expected to be refactored over time to be more generic, rather than hardcoded. The only part that will remain implementation-specific is the database processing logic, which needs to be rewritten for each distinct database, and a small part of the analysis.  
Additionally, this project currently lacks testing and exception handling. Since it was primarily developed for internal data analysis and not intended as a widely used tool, these aspects were initially deprioritized. For the same reason the documentation is also kept to a minimum.


Following is the underlying theory of the project study:

### RNA folding and architecture

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


### RNA structure prediction

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


### CRISPR RNA

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


### Bibliography

[1] I. Tinoco und C. Bustamante, „How RNA folds“, J. Mol. Biol., Bd. 293, Nr. 2, S. 271–281,
Okt. 1999, doi: 10.1006/jmbi.1999.3001.  
[2] W. A. Haseltine, K. Hazel, und R. Patarca, „RNA Structure: Past, Future, and Gene Therapy
Applications“, Int. J. Mol. Sci., Bd. 26, Nr. 1, S. 110, Dez. 2024, doi: 10.3390/ijms26010110.  
[3] L. G. Parlea, B. A. Sweeney, M. Hosseini-Asanjan, C. L. Zirbel, und N. B. Leontis, „The RNA
3D Motif Atlas: Computational methods for extraction, organization and evaluation of
RNA motifs“, Methods San Diego Calif, Bd. 103, S. 99–119, Juli 2016, doi:
10.1016/j.ymeth.2016.04.025.  
[4] A. Lescoute, N. B. Leontis, C. Massire, und E. Westhof, „Recurrent structural RNA motifs,
Isostericity Matrices and sequence alignments“, Nucleic Acids Res., Bd. 33, Nr. 8, S. 2395–
2409, 2005, doi: 10.1093/nar/gki535.  
[5] L. Nasalean, J. Stombaugh, C. L. Zirbel, und N. B. Leontis, „RNA 3D Structural Motifs:
Definition, Identification, Annotation, and Database Searching“, in Non-Protein Coding
RNAs, N. G. Walter, S. A. Woodson, und R. T. Batey, Hrsg., Berlin, Heidelberg: Springer,
2009, S. 1–26. doi: 10.1007/978-3-540-70840-7_1.  
[6] I. L. Hofacker, W. Fontana, P. F. Stadler, L. S. Bonhoeffer, M. Tacker, und P. Schuster, „Fast
folding and comparison of RNA secondary structures“, Monatshefte Für Chem. Chem.
Mon., Bd. 125, Nr. 2, S. 167–188, Feb. 1994, doi: 10.1007/BF00818163.  
[7] M. Zuker und P. Stiegler, „Optimal computer folding of large RNA sequences using
thermodynamics and auxiliary information“, Nucleic Acids Res., Bd. 9, Nr. 1, S. 133–148,
Jan. 1981, doi: 10.1093/nar/9.1.133.  
[8] B. Voß, „Classified Dynamic Programming in RNA Structure Analysis“, Methods Mol. Biol.
Clifton NJ, Bd. 2726, S. 125–141, 2024, doi: 10.1007/978-1-0716-3519-3_6.  
[9] R. Giegerich, B. Voss, und M. Rehmsmeier, „Abstract shapes of RNA“, Nucleic Acids Res.,
Bd. 32, Nr. 16, S. 4843–4851, 2004, doi: 10.1093/nar/gkh779.  
[10] J. Huang, R. Backofen, und B. Voß, „Abstract folding space analysis based on helices“,
RNA, Bd. 18, Nr. 12, S. 2135–2147, Dez. 2012, doi: 10.1261/rna.033548.112.  
[11] M. Sebeke und B. Voß, „Incorporating RNA 3D Motifs into RNA secondary structure
prediction“.  
[12] J. Reeks, J. H. Naismith, und M. F. White, „CRISPR interference: a structural
perspective“, Biochem. J., Bd. 453, Nr. Pt 2, S. 155–166, Juli 2013, doi:
10.1042/BJ20130316.  
[13] K. Chylinski, A. Le Rhun, und E. Charpentier, „The tracrRNA and Cas9 families of type II
CRISPR-Cas immunity systems“, RNA Biol., Bd. 10, Nr. 5, S. 726–737, Mai 2013, doi:
10.4161/rna.24321.  
[14] E. Charpentier, H. Richter, J. van der Oost, und M. F. White, „Biogenesis pathways of
RNA guides in archaeal and bacterial CRISPR-Cas adaptive immunity“, FEMS Microbiol.
Rev., Bd. 39, Nr. 3, S. 428–441, Mai 2015, doi: 10.1093/femsre/fuv023.  
[15] V. Kunin, R. Sorek, und P. Hugenholtz, „Evolutionary conservation of sequence and
secondary structures in CRISPR repeats“, Genome Biol., Bd. 8, Nr. 4, S. R61, Apr. 2007, doi:
10.1186/gb-2007-8-4-r61.