# pairwise_sequence calculator 

Running sequence comparision at position i and j based on physiochemical score: <br>
 1)  Structure of proteins: packing of alpha-helices and pleated sheets. PNAS (1977), Vol. 74, No 10. p. 4130-4134. <br>
     C. Chothia, M. Levitt, D. Richardson <br>
 2)  Based on the review by Cyrus Chothia pioneer of structural bioinformatics <br>
     Average Volume per residue taken from Ann. Rev. Biochem. (1984) 53, 537-72.<br>
 3)  Conformation of amino acid side-chains in proteins <br>
     J Janin, S. Wodak, M. Levitt, B. Maigret from JMB (1978) 124, 357-386 <br>
 3)  Code is based on this great paper: Dmitry A. Afonnikov and Nikolay A. Kolchanov. CRASP: a program for analysis of   
     coordinated substitutions in multiple alignments of protein sequence. Nucleic acids res.(2004) W64-W67.<br>
 
     
 <div>
    <div></div>
    <hr class="styled-hr" />
    <div></div>
 </div>

 (Step1) align_pairwise('seq1','ChotVol.txt',0.8,0,-0.8,0) <br>
         [matlab will automatically generate a heatmap] <br>
   
 (Step2) run heatmap_raindbow.r to generate Heatmap using ggplot2 <br>

