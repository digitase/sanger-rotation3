'''Read in and annotate known homoplasies in st239

Reference:
    Harris, S.R., Feil, E.J., Holden, M.T., Quail, M.A., Nickerson, E.K.,
    Chantratita, N., Gardete, S., Tavares, A., Day, N., Lindsay, J.A. and
    Edgeworth, J.D., 2010. Evolution of MRSA during hospital transmission and
    intercontinental spread. Science, 327(5964), pp.469-474.
'''

import pandas as pd
import sys
import os
sys.path.append(os.path.abspath("/nfs/users/nfs_b/bb9/workspace/rotation3/src/2_homoplasy/"))
import _summarise_homoplasies

known_hps = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/misc/homoplasies_st239_harris_et_al_2010_science_table1.csv', index_col=0)
known_hps = known_hps.set_index('SNP position')

annot = _summarise_homoplasies.annotate_homoplasies(known_hps, "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl")

annot.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/misc/homoplasies_st239_harris_et_al_2010_science_table1.annotated.csv')

