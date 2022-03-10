## TF-seq-matching
Python script for the following works:
* read Transcription Factor (TF) data files "1_sorted.csv" ~ "20_sorted.csv" and make a TF sequence database
* read a promoter sequence "SIRT3_promoter.txt" and find all TF sequences and their positions in the promoter (fast string matching algorithm is used)
* write matching results to "result.csv"

Modify the data file paths in "seq.py". 
