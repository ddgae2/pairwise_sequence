# pairwise_sequence calculator 

Running sequence comparision at position i and j based on physiochemical score: <br>

 1)  Based on the review by Cyrus Chothia <br>
     Average Volume per residue taken from Ann. Rev. Biochem. (1984) 53, 537-72.<br>
 2)  Code is based on this great paper: Dmitry A. Afonnikov and Nikolay A. Kolchanov. CRASP: a program for analysis of   
     coordinated substitutions in multiple alignments of protein sequence. Nucleic acids res.(2004) W64-W67.<br>
     
 <div>
    <div></div>
    <hr class="styled-hr" />
    <div></div>
 </div>

 (Step1) align_pairwise('seq1','ChotVol.txt',0.8,0,-0.8,0) <br>
         [matlab will automatically generate a heatmap] <br>
   
 (Step2) run heatmap_raindbow.r to generate Heatmap using ggplot2 <br>

