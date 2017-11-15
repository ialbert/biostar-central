import glob
import os


def collect_paths(base_dir, target_dir, file_patt):
    """
    inputs are base url, target directory and file pattern
    return a list of tuples of the form [(complete file path, filename),..]

    """
    patt = os.path.join(target_dir, file_patt)

    store = []

    for fname in glob.glob(patt):

        name = os.path.basename(fname)
        path = os.path.join(base_dir, name)

        store.append((path, name))

    return store
