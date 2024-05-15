from Bio.Align.Applications import ClustalOmegaCommandline

input_file = "input.txt"
output_file = "aligned.txt"

clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
clustalomega_cline()