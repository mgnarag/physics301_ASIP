#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 19:07:24 2022

@author: markjeremynarag
"""

from PIL import Image
import numpy as np
from sklearn.decomposition import PCA

#...

img = Image.open('narag_1.png')
width, height = image.size
#img.show()

pca = PCA(n_components=(width*height))
