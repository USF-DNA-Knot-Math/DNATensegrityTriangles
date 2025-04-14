This program is intended for use in calculating geometric and topological parameters of DNA chains from PDB files. 

There are two .m files for use, the first is a script wrapper, geoScript.m, and the second, GaussFrenet.m is the core function. Examples are included in geoScript.m for how to calculate on all PDB files from the Engineering Chirality manuscript.

GaussFrenet is a function file with several inputs: 

'filename' should correspond to a PDB file in the PDBs directory. All PDB files from the manuscript are included in the PDBs folder by default. For example: "2T3_merged.pdb"

'specialChainID' refers to which chain the user wishes to analyze. This can be determined by opening the PDB file in a graphic user interface, i.e. COOT, ChimeraX or the Web3X DNA server and to determine which chain is of interest. Typically, this would be the center strand of a DNA triangle. 

'circBool' governs whether there should be virtual closure of the strand. This would apply for center strands of tensegrity triangles. Unless actual ligation is carried out to circularize the strand, this is a theoretical consideration in order to determine the integer linking number of a motif. 

'allPlotBool' governs how many graphs to print while analyzing. 

'targetAtom' dictates which atom is pulled for the primary chain, or 'K' in Appendix B. This is, by default, the O3' atom in the phosphate backbone. It can be changed to others for non-standard polymers, etc. The second chain, 'L' from Appendix B, is set by default to the primary nitrogen connecting the sugar to the nucleotide. If non-standard residues are used, this will fail without modifying the code. Currently the code accepts all D- and L-DNA residue names.

'plotBaseNum' is an output parameter that governs the first "n" number of residues on specialChainID that are plotted to generate the XY plots of curvature, torsion, link, twist and writhe. This can be set to 0 to enable all residues to display. For the graphs in the manuscript, this was set to the length of the center strand in the asymmetric unit, i.e. 13, 14, 17, 24 nucleotides, etc.

Computation data is saved in the LinkData folder, and graphs are saved in the Figs folder, using the prefix included in the PDB name. These do not include a timestamp, and will overwrite if run on the same file.

Example usages of GaussFrenet.m include:

GaussFrenet("3T14_merged.pdb",'D',0,1, {'O3'''},14) 
% runs without virtual closure of the nick, plots all graphs with 14 bp steps

GaussFrenet("2T7_merged.pdb",'B',1,0, {'O3'''},7) 
% closes the nick for integer linking, does not display some graphs, uses 7 bp xrange

GaussFrenet("2T7_merged.pdb",'B',1,1, {'O3'''},7) 
% closes the nick for integer linking, displays all graphs, uses full 21 bp xrange
