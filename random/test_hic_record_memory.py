import time

import hicstraw

hic_resolution = 5000
hic_file = "https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic"

hic = hicstraw.HiCFile(hic_file)
chromosome = "chr1"
matrix_object = hic.getMatrixZoomData(
    chromosome, chromosome, "observed", "SCALE", "BP", hic_resolution
)
# records = matrix_object.getRecords(0, 0, 0, 0)
records = matrix_object.getRecords(0, 75000, 0, 248956422)
import pdb

pdb.set_trace()
print(len(records))
