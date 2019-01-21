import subprocess
import sys
import os
import concurrent
import concurrent.futures
import ftplib
import pandas as pd
import re

def get_current_remote_version():
    """
    :return: fetches newest PDB db version from remote repository
    """
    ftp = ftplib.FTP("ftp.wwpdb.org")
    ftp.login("anonymous", "")
    raw_data = []
    ftp.dir('pub/pdb/data/status/', raw_data.append)
    dates = [entry.split(' ')[-1] for entry in raw_data if entry.split(' ')[-1].isnumeric()]
    latest_version = sorted(dates)[-1]
    return latest_version


def is_nucl_seq(seq):
    nuc_re = re.compile(r'[^ATGCUXNIF.]')
    res = nuc_re.search(seq)
    return not bool(res)


def is_x_seq(seq):
    x_re = re.compile(r'[^X.]')
    res = x_re.search(seq)
    return not bool(res)


def get_current_local_version(db_path):
    """
    :param db_path: location of the PDB db
    :return: newest version of PDB in the local repository
    """
    return sorted(os.listdir('{}/data/'.format(db_path)))[-1]


def get_previous_local_version(db_path, curr_version):
    """

    :param db_path:
    :param curr_version:
    :return:
    """
    dirs = [ver for ver in os.listdir('{}/data/'.format(db_path)) if ver < curr_version]
    return sorted(dirs)[-1]

def get_first_local_version(db_path):
    """

    :param db_path:
    :return:
    """
    return sorted(os.listdir('{}/data/'.format(db_path)))[0]


def fetch_cluster_data(dn_path):
    results = set()
    fns = ['bc-100.out', 'bc-95.out', 'bc-90.out', 'bc-70.out', 'bc-50.out', 'bc-40.out', 'bc-30.out',
           'clusters95.txt', 'clusters90.txt', 'clusters70.txt', 'clusters50.txt',
           'current-mmseqs/clusters-by-entity-100.txt', 'current-mmseqs/clusters-by-entity-95.txt',
           'current-mmseqs/clusters-by-entity-90.txt', 'current-mmseqs/clusters-by-entity-70.txt',
           'current-mmseqs/clusters-by-entity-50.txt', 'current-mmseqs/clusters-by-entity-40.txt',
           'current-mmseqs/clusters-by-entity-30.txt']
    for fn in fns:
        wget_cmd = 'wget -P {} ftp://resources.rcsb.org/sequence/clusters/{}'.format(dn_path, fn)
        res = os_cmd(wget_cmd)
        results.add(res)
    if len(results) == 1 and 0 in results:
        # print('\033[92mOK!\033[0m')
        return 0
    else:
        # print('\033[91mFailed!\033[0m')
        return 1


def parse_cluster_data(fn):
    f = open(fn, 'r')
    data = [line.rstrip() for line in f.readlines()]
    f.close()
    cluster_data = {}
    for c in range(1, len(data)+1):
        for entry in data[c-1].split(' '):
            pdb, chain = entry.split('_')
            cluster_data['{}_{}'.format(pdb.lower(), chain)] = str(c)
    return cluster_data


def os_cmd(cmd, check=False):
    result = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if check:
        if result == 0:
            pass
            #print('\033[92mOK!\033[0m')
        else:
            pass
            sys.exit(1)
    return result


def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass


def check_socket_output(fn):
    try:
        f = open(fn, 'rb')
        f.seek(-1024, 2)
        last = f.readlines()[-2].decode()
        f.close()
        if last.endswith('NO COILED COILS\n'):
            return False
        else:
            return True
    except (TypeError, KeyError, OSError):
        return False


def prep_paths(main_path, ver=None, pdbs=None):
    mkdir("{}/structures".format(main_path))
    mkdir("{}/biounits".format(main_path))
    mkdir("{}/data".format(main_path))
    if ver:
        mkdir("{}/data/{}".format(main_path, ver))
        mkdir("{}/data/{}/clusters".format(main_path, ver))
    #mkdir("{}/dssp_pdb".format(main_path))
    #mkdir("{}/dssp_pdb_biounit".format(main_path))
    #mkdir("{}/socket_pdb_biounit".format(main_path))
    #mkdir("{}/data/init/failed".format(main_path))
    #mkdir("{}/master_pdb".format(main_path))
    #mkdir("{}/master_pdb_biounit".format(main_path))
    if pdbs:
        for pdb in pdbs:
            mkdir("{}/structures/{}".format(main_path, pdb[1:3]))
            mkdir("{}/biounits/{}".format(main_path, pdb[1:3]))
            #mkdir("{}/dssp_pdb/{}".format(main_path, pdb[1:3]))
            ##mkdir("{}/dssp_pdb_biounit/{}".format(main_path, pdb[1:3]))
            #mkdir("{}/socket_pdb_biounit/{}".format(main_path, pdb[1:3]))
            #mkdir("{}/master_pdb/{}".format(main_path, pdb[1:3]))
            #mkdir("{}/master_pdb_biounit/{}".format(main_path, pdb[1:3]))


def multiprocess(func, cmd_dict, np=20, return_type='', print_progress=True, ok_status=0):
    """
    :param func: function to parallelize
    :param cmd_dict: dictionary with job id's (keys) and 'func' inputs as values
    :param np: number of processes
    :param return_type: 'failed' to return only failed job ids (according to ok_status' or 'all' for all results
    :param print_progress: print simple progress bar
    :param ok_status: value returned by 'func' indicating everything went ok
    :return: set of failed job ids or dict with ids and 'func' outputs
    """
    failed_jobs = set()
    job_count = len(cmd_dict)
    count = 0
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=np) as executor:
        futures = {executor.submit(func, cmd): ident for ident, cmd in cmd_dict.items()}
        for future in concurrent.futures.as_completed(futures):
            ident = futures[future]
            if return_type == 'failed':
                if future.result() != ok_status:
                    failed_jobs.add(ident)
            elif return_type == 'all':
                results[ident] = future.results()
            count += 1
            if print_progress:
                sys.stdout.write("\r%s/%s" % (count, job_count))
                sys.stdout.flush()
    if print_progress:
        print("")
    if return_type == 'failed':
        return failed_jobs
    elif return_type == 'all':
        return results
