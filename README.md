# Virus_Classification_FFP
Real-time analytics for large-scale classification of viruses using Feature Frequency Profile method

1. lists for 3905 viruses:

2. jellyfish commands: -m: k length, -s memory, -t threads and -o output, the target file is filename.txt
     
    ./jellyfish/bin/jellyfish count -m 7 -s 300M -t 200 input.fasta -o filename.jf
    ./jellyfish/bin/jellyfish dump -c filename.jf > filename.txt

    ( and the tutorial is here: http://www.genome.umd.edu/docs/JellyfishUserGuide.pdf) 
3. CRE: Cumulative Relative Entropy
4. FFP python file
5. Tree distance calculation in R: RF_dist.R

Please cite: Qian Zhang, Se-Ran Jun, Michael Leuze, David Ussery, Intawat Nookaew. Viral phylogenomics using an alignment-free method: A three-step approach to determine optimal length of k-mer. Scientific reports 7, 40712.https://www.nature.com/articles/srep40712  
