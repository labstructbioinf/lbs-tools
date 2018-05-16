import numpy as np


def prep_okeanos_jobs(joblist, n, resources='cores', n_res_per_job=1, account='GA67-18', sh_dir='_SH'):
    """
    Prepares the job files for the parallel run on Okeanos cluster
    :param joblist: list of jobs to calculate
    :param n: number of resources (cores or nodes) for the run
    :param resources: types of resources (cores or nodes)
    :param n_res_per_job: number of dedicated resources per job
    :param account: Okeanos account
    :param sh_dir: directory with the output .sh files
    :return:
    """
    assert(resources == 'cores' or resources == 'nodes')
    if resources == 'cores':
        max_res = n*24/n_res_per_job
    elif resources == 'nodes':
        max_res = n/n_res_per_job
    k = np.array(joblist)
    chunks = np.array_split(k, max_res)
    try:
        chunks = list(filter(None, chunks))
    except ValueError:
        pass
    for i in range(0, len(chunks)):
        f = open('%s/%s.sh' % (sh_dir, i), 'w')
        f.write('#!/bin/bash\n')
        for chunk in chunks[i]:
            f.write('%s\n' % (chunk))
        f.close()
