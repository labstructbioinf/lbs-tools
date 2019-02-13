import logging
import os
def logger_setter(module,file_name=''):
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)
    logger = logging.Logger(os.path.basename(module))
    if file_name:
        fh = logging.FileHandler()
        fh.setLevel(logging.ERROR)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    logger.addHandler(handler)
    return logger


