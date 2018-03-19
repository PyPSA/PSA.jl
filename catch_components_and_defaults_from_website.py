#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 18:32:30 2018

@author: fabian
"""

import pandas as pd
import re
url = 'https://www.pypsa.org/doc/components.html'

def to_numeric(x): 
    s = pd.Series(x)
    return pd.to_numeric(s, errors='ignore')[0]

pypsadoc = pd.read_html(url, converters={'default':to_numeric})



def to_initializer(df):
    attr = (":" + df.attribute.astype(str)).tolist()
    return re.sub("'", "" , "DataFrame(repeat([Bool[]], outer=%s), \n %s)"%(len(df), attr))

component_map = pd.Series({
        'buses':3,
        'carriers':4,
        'global_constraints':5,
        'generators':6,
        'storage_units':7,
        'stores':8,
        'loads':9,
        'shunt_impendances':10,
        'lines':11,
        'line_types':12,
        'transformers':14,
        'transformer_types':15,
        'links':17,
        })

for comp in component_map.index:
    print comp + ' = ' + to_initializer(pypsadoc[component_map[comp]]) + '\n'


import pypsa
network = pypsa.Network()

#%%
timedependents = pd.Series(dir(network))[lambda df:df.str.find('_t', -2) > 0]

for comp in timedependents:
    descr = str(zip(getattr(network, comp), getattr(network, comp).values()))
    descr = re.sub('Empty DataFrame\n', 'DataFrame()', descr)
    descr = re.sub('Columns: \[\]\n', '', descr)
    descr = re.sub('Index: \[now\]', '', descr)
    descr = re.sub("'", '"' ,descr)
    print comp + '= Dict{String,DataFrame}( ' + descr + ')'
