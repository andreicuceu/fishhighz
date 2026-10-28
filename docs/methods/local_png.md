# Leading local-PNG forecast integration

This sibling is based on the unchanged full-EFT integration branch
d68a591c50317857dbc78a66f86659233f9594ee. Its optional companion is
eftkernels branch codex/lya-galaxy-fnl, based on
f5f8427b49814dd0e910a7a5aa60ed53e74cff09.

The companion pure-array local_png module adds the leading scale-dependent
bias to a Gaussian tree or Gaussian one-loop mean. It is not a full PNG
one-loop EFT model. Use the generic P3DProvider interface to expose fNL and
the applicable nuisance parameters. Keep signed forest/galaxy kernels and
conditional response coefficients explicit. Known observational noise stays
in the covariance and must not be added to the intrinsic mean twice.

For partial forest overlap, add nonoverlap QSO-auto information and the joint
overlap FF/FQ/QQ Fisher in a common unmarginalized registry. The QSO overlap
must be counted once. Select spectrum covariance before inversion, preserve
null directions, and marginalize all active linear/EFT biases. Historical
DESI forecast source scripts and dependency pins are archived in the private
companion repository; survey inputs and posterior archives are not public
code assets.
