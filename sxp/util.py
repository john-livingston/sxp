import os
import sys
import shutil


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if shutil.fnmatch.fnmatch(basename, pattern):
                filepath = os.path.join(root, basename)
                yield filepath


def ignore_oserror(f):
    """
    rather silly decorator for ignoring OSError exceptions
    """
    def wrapper(*args):
        try:
            f(*args)
        except OSError as e:
            print(e)
    return wrapper


@ignore_oserror
def mkdir(d):
    os.makedirs(d)


def mkdirs(*args):
    for d in args:
        mkdir(d)
