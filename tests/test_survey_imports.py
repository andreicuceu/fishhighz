"""Core paths stay usable when optional imports are actively blocked."""

import subprocess
import sys


def test_blocked_optional_core_and_adapter():
    script = r"""
import sys
import importlib.abc
class BlockOptional(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'astropy','scipy'}:
            raise ImportError('blocked optional support for test')
sys.meta_path.insert(0, BlockOptional())
import numpy as np
import fishhighz
from fishhighz.geometry import prepare_geometry, prepare_astropy_geometry
from fishhighz.response import velocity_response
from fishhighz.models.p1d import default_p1d
from fishhighz.models.external import evaluate_p1d, BoundParameters
from fishhighz.parameters import Parameter, ParameterRegistry
assert not {'astropy','scipy'} & set(sys.modules)
g = prepare_geometry(2,3,z_eval=2.4,area_deg2=1,h_fid=.7,z_order=4,
    hubble=lambda z:z*100,transverse_distance=lambda z:z*1000)
assert g.volume > 0
assert velocity_response([0],pixel_width_velocity=0,gaussian_sigma_velocity=0)[0] == 1
assert default_p1d([],3,[0])[0] > 0
r=ParameterRegistry([Parameter('a',1,'target')])
assert evaluate_p1d(lambda t,z,k: t[0]+k, BoundParameters(r,['A'],{'A':'a'}),[1],3,[0])[0] == 1
assert not {'astropy','scipy'} & set(sys.modules)
try:
    prepare_astropy_geometry(None,2,3,z_eval=2.4,area_deg2=1,h_fid=.7,z_order=4)
except ImportError as e:
    assert 'fishhighz[cosmology]' in str(e)
else:
    raise AssertionError('adapter must require optional dependencies')
"""
    completed = subprocess.run(
        [sys.executable, "-I", "-c", script], capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr
    assert not completed.stdout
