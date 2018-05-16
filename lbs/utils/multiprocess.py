import subprocess
import multiprocessing
import os


def os_cmd(cmd):
    subprocess.call([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #os.system(cmd)


def multiprocess(func, tasks, n_cores, tasks_per_core=1):
    stdout_queue = multiprocessing.Queue()
    pool = multiprocessing.Pool(
        processes=n_cores, initargs=[stdout_queue], maxtasksperchild=tasks_per_core)
    for i, data in enumerate(pool.imap(func, tasks), 1):
        pass
    pool.close()
    pool.join()
