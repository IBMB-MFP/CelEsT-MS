control=$(awk -F '\t' {print $9} input/benchmark_SRA_colourspace)
treatment=$(awk -F '\t' {print $10} input/benchmark_SRA_colourspace)