## TF-seq-matching
Python script for the following workflow:
* read Transcription Factor (TF) data files "1_sorted.csv" ~ "20_sorted.csv" to make a TF sequence database
* read a promoter sequence "SIRT3_promoter.txt" and find all TF sequence and corresponding positions (fast string matching algorithm is used)
* write results to "result.csv"

To reproduce the results, modify the file paths used in the script "seq.py". 
