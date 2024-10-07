import pandas as pd
import numpy as np
import os
import tables


if __name__ == '__main__':
    file_name = '/home/cleon/PycharmProjects/MSc-Project/out/mapscape/933124/results.h5'

    # os.chmod(file_name, 0o777)

    # h5 file
    f = tables.open_file(file_name, mode='r')

    table1 = f.root.results.optimal

    # select the best tree
    best = list(np.array(table1.index))[list(np.array(table1.values)).index(True)]
    print(best)

    # array = f['/results']
    results = f.root.results
    print(list(results))
    print(list(f.root.results.num_mutations.values))

    # bic = f.root.results[best].bic
    # print(bic)

    # error_rate = f.root.results.error_rate
    # print(error_rate)

    # likelihood = f.root.results.likelihood
    # print(likelihood)

    # num_mutations = f.root.results.num_mutations
    # print(num_mutations)

    # num_nodes = f.root.results.num_nodes
    # print(num_nodes)

    # num_samples = f.root.results.num_samples
    # print(num_samples)

    # objective_value = f.root.results.objective_value
    # print(objective_value)

    # optimal = f.root.results.optimal
    # print(optimal)

    # tree_id = f.root.results.tree_id
    # print(tree_id)

    # tree_index = f.root.results.tree_index
    # print(tree_index)

    # tree_string = f.root.results.tree_string
    # print(tree_string)

    f.close()
