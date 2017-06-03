
import pandas as pd

def flatten_list_or_tuple(xs):
    if isinstance(xs, (list, tuple)):
        for x in xs:
            yield from flatten_list_or_tuple(x)
    else:
        yield xs

def tag_known_loci(candidates, known_amr_loci):
    for gs in candidates:
        gs = flatten_list_or_tuple(gs)
        gs_result = set()
        for g in gs:
            if g:
                for known in known_amr_loci:
                    if known.startswith(g.lower()[:3]):
                        gs_result.add(known)
                # print(process.extract(g.lower(), known_amr_loci, limit=2))
        yield gs_result

with open('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_annotate_candidates/known_amr_loci.txt', 'r') as infile:
    known_amr_loci = [x.strip().lower() for x in infile.readlines()]

hp_files = {
    'st22': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/st22_homoplasies.csv', 
    'st239': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/st239_homoplasies.csv', 
    #  'st30': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST30_BSAC_Pfizer/st30_homoplasies.csv', 
    'st8': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/st8_homoplasies.csv' 
}

for prefix in hp_files.keys():
    #  prefix = 'st22'
    hp = pd.read_csv( 
        hp_files[prefix], 
        converters={'gene': eval}
    )
    hp['known_amr_loci'] = list(tag_known_loci(hp['gene'], known_amr_loci))
    hp.loc[hp['known_amr_loci'].map(len) > 0].to_csv(
        '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/{}_homoplasies.known_amr_filtered.csv'.format(prefix)
    )

