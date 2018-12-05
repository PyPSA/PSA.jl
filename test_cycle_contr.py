#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:48:31 2018

@author: fabian
"""

import pypsa
import pandas as pd

name = "elmod_8760h"

n = pypsa.Network("/home/vres/data/pypsa_models/elmod/{}.nc"
                  .format(name))
n.set_snapshots(n.snapshots[1:2])

n.lopf()

