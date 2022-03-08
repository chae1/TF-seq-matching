## TF-seq-matching
Workflow
* read "1_sorted.csv" ~ "20_sorted.csv" to make a Transcription Factor (TF) sequence database
* read a promoter sequence "SIRT3_promoter.txt" and find all TF and their positions by referencing database (fast string matching algorithm is used)
* write results to "result.csv"

Note: Modify file paths before use in "seq.py"
