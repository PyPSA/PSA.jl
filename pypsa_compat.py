
#Create a library with its own Network class that has a new CHP
#component and a new LOPF function for the CHP constraints

#NB: This only works with Python 3 because of super()

import pypsa, pandas as pd, numpy as np

from pypsa.descriptors import Dict




override_components = pypsa.components.components.copy()
override_component_attrs = Dict({k : v.copy() for k,v in pypsa.components.component_attrs.items()})

override_component_attrs['Generator'].loc['p_nom_opt', 'type'] = 'series'
override_component_attrs['Line'].loc['s_nom_opt', 'type'] = 'series'
override_component_attrs['Link'].loc['p_nom_opt', 'type'] = 'series'


class Network(pypsa.Network):
    def __init__(self,*args,**kwargs):
        kwargs["override_components"] = override_components
        kwargs["override_component_attrs"] = override_component_attrs
        super().__init__(*args,**kwargs)
