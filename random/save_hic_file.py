import hicstraw
import pandas as pd

HIC_FILE = "https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic"
CHROMOSOME = "chr22"
RESOLUTION = 5000

result = hicstraw.straw("observed", "SCALE", HIC_FILE, CHROMOSOME, CHROMOSOME, "BP", RESOLUTION)
bin_data = [[r.binX, r.binY, r.counts] for r in result]
df = pd.DataFrame(bin_data)
df.to_csv(f"{CHROMOSOME}_SCALE.tsv.gz", sep='\t', index=False, header=False, compression='gzip')