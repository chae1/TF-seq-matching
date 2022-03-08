## TF-seq-matching
Python script for following workflow:
* read Transcription Factor (TF) data "1_sorted.csv" ~ "20_sorted.csv" to make a TF sequence database
* read a promoter sequence "SIRT3_promoter.txt" and find all TF and corresponding positions (fast string matching algorithm is used)
* write results to "result.csv"

To reproduce the results, modify the file paths used in the script "seq.py".
