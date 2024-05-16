from Bio.Align.Applications import ClustalOmegaCommandline

input_file = "io/input_HLAA.txt"
output_file = "io/aligned_HLAA.txt"

clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
clustalomega_cline()

input_file = "io/input_HLAB.txt"
output_file = "io/aligned_HLAB.txt"

clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
clustalomega_cline()

input_file = "io/input_HLAC.txt"
output_file = "io/aligned_HLAC.txt"

clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
clustalomega_cline()