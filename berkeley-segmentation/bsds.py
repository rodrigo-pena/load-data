# -*- coding: utf-8 -*-

r"""
This module implements functions used to load and manipulate data from the
Berkeley segmentation dataset:

[Link] https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/

"""

import numpy as np
import pandas as pd


def mask_from_human_annotation(dir, id):
    r"""
    Make segmentation mask from human annotation file.

    Parameters
    ----------
    dir : string
        Full path of the directory of the annotations.
    id : string
        Name of the annotation file

    Returns
    -------
    mask : numpy.ndarray
        Mask generated from the annotation file
    n_seg : int
        Number of segments
    """
    try:
        full_path = dir + id + str('.seg')
    except TypeError:
        print('Both dir and id should be strings.')

    header = pd.read_table(full_path, names=[0, 1], sep=' ', nrows=11)
    width = int(header[0]['width'])
    height = int(header[0]['height'])
    n_seg = int(header[0]['segments'])

    data = pd.read_table(full_path, skiprows=11, names=[0, 1, 2, 3], sep=' ')
    data = data.values

    count = 0
    mask = np.zeros((height, width))
    for i in data[:, 1]:
        mask[i, data[count, 2]:data[count, 3]] = data[count, 0]
        count += 1

    return mask, n_seg
