"""Actively reject optional imports during normalized multi-bin orchestration."""

import os
import subprocess
import sys

# Also extracted as a standalone installed-wheel probe (no pytest dependency).
BASE_PROGRAM = r"""
import importlib.abc
import sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in {'scipy','astropy'}:
            raise ImportError('optional libraries deliberately blocked')
sys.meta_path.insert(0,Block())
import numpy as np
from fishhighz.fields import ObservedField,PairSelection
from fishhighz.geometry import prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters,P3DProvider,PreparedP3D
from fishhighz.parameters import Parameter,ParameterRegistry
from fishhighz.response import InstrumentResponse
from fishhighz.survey import BinSpec,ForestInput
from fishhighz.forecast import prepare_bin,run_forecast
from fishhighz.adapters.legacy_inputs import DensityReader,SNRReader
registry=ParameterRegistry([Parameter('A',2,'target',step=.01),Parameter('b',0,'nuisance',step=.01)])
field=ObservedField('g','galaxy','g')
selection=PairSelection([field])
grid=gauss_legendre_grid([.02,.2],k_order=2,mu_order=3,h_fid=.7)
binding=BoundParameters(registry,('a','b'),{'a':'A','b':'b'})
bins=[]
for i,s in enumerate([1,-1]):
    geometry=prepare_geometry(2+i,3+i,z_eval=2.5+i,area_deg2=10,h_fid=.7,z_order=4,
        hubble=lambda z:np.full_like(z,200),transverse_distance=lambda z:1000*(1+z))
    model=lambda t,z,k,mu,p,s=s:np.full((len(k),len(p)),t[0]+s*t[1])
    p3d=PreparedP3D(registry,selection,[P3DProvider('scalar',model,binding,[(0,0)])])
    bins.append(prepare_bin(BinSpec(str(i),geometry,grid,p3d,{'g':InstrumentResponse(0,0)},
        galaxies={'g':1},independent_sampling=True)))
prior=np.diag([0,4.])
result=run_forecast(bins,prior_fisher=prior,batch_size=5)
expected=np.zeros((2,2))
for b,s,r in zip(bins,[1,-1],result.bins):
    f=b.modes.sum()/18*np.array([[1,s],[s,1]])
    np.testing.assert_allclose(r.result.data_fisher,f,rtol=5e-13,atol=0)
    assert r.result.diagnostics.rank==1
    expected+=f
np.testing.assert_allclose(result.combined.data_fisher,expected,rtol=5e-13,atol=0)
np.testing.assert_array_equal(result.combined.prior_fisher,prior)
assert result.combined.diagnostics.rank==2
# A normalized forest goes through independent external P1D and fixed weights.
f=ObservedField('f','forest','lya','qso')
selection=PairSelection([f])
p3d=PreparedP3D(registry,selection,[P3DProvider('forest',model,binding,[(0,0)])])
source=ForestInput(dict(z_source=5,magnitudes=[20,21],quadrature=[1,1],rho=[.01,.02],variance=[1,2],
    length_velocity=10000,method='supplied',weights=[1,1]),lambda t,z,k:np.ones_like(k),BoundParameters(registry,(),{}),registry.fiducials)
b=prepare_bin(BinSpec('forest',geometry,grid,p3d,{'f':InstrumentResponse(30,10)},forests={'f':source},independent_sampling=True))
assert run_forecast([b]).combined.diagnostics.rank==1
for reader,args in [(DensityReader,dict(path='unused',semantics='cell_count_per_deg2',target_density=None,z_norm_min=None)),
                    (SNRReader,dict(paths=['unused'],smoothing='legacy'))]:
    try: reader(**args)
    except ImportError as error: assert 'fishhighz[survey]' in str(error)
    else: raise AssertionError('missing survey extra must fail at setup')
assert not {'scipy','astropy','vega','lyaforecast','camb'} & set(sys.modules)
print('NumPy-only two-bin/scalar and normalized forest oracles passed; both missing-extra errors passed')
"""


def test_blocked_optional_multi_bin():
    result = subprocess.run(
        [sys.executable, "-I", "-c", BASE_PROGRAM],
        capture_output=True,
        text=True,
        env=dict(
            os.environ,
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
        ),
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_raw_example():
    import runpy
    from pathlib import Path

    report = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "examples/survey_forecast.py")
    )["run"]()
    assert report["bin_ids"] == ("low", "high")
    assert report["weights"][0]["fq"] != report["weights"][0]["fl"]
