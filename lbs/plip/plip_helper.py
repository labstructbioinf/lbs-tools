import multiprocessing

# for multiprocess
def run_multiprocess(func, tasks, n_cores, tasks_per_core=1):
    stdout_queue = multiprocessing.Queue()
    pool = multiprocessing.Pool(
        processes=n_cores, initargs=[stdout_queue], maxtasksperchild=tasks_per_core)
    for i, data in enumerate(pool.starmap(func, tasks), 1):
        yield data
    pool.close()
    pool.join()
