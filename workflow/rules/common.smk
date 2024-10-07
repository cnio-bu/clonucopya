import pandas as pd
from pathlib import Path

##### Wildcard constraints #####
# wildcard_constraints:
#    sample="|".join(SAMPLES.index)

##### Helper functions #####

def get_resource(rule,resource):
    try:
        return config["resources"][rule][resource]
    except KeyError:
        return config["resources"]["default"][resource]


def get_cluster_ids():
    if Path(f"/home/cleon/PycharmProjects/MSc-Project/out/pyclone/output/tables/cluster.tsv").exists():
        clusters = pd.read_csv(f"/home/cleon/PycharmProjects/MSc-Project/out/pyclone/output/tables/cluster.tsv", sep="\t").set_index("cluster_id", drop=False)
        return [getattr(row, 'cluster_id') for row in clusters.itertuples()]
    else:
        return ['0']
