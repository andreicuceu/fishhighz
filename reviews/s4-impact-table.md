# S4 joint BAO effects from saved direct controls

Each range spans bins 1–6; components are radial and transverse. Forward changes introduce one accuracy setting into fixed-compatibility. Reverse changes restore one legacy setting in accuracy; their signs are **not reversed**. Percentages are relative to that direction’s own endpoint. They are conditional effects, not additive shares of the endpoint difference.

| Setting (legacy → accuracy) | Forward radial (%) | Forward transverse (%) | Reverse radial (%) | Reverse transverse (%) |
| --- | ---: | ---: | ---: | ---: |
| 107-node full-endpoint rectangular → composite Gauss–Legendre | -0.0137684 to +0.000425783 | -0.0156424 to +0.000430085 | -0.0107707 to +0.147708 | -0.0181169 to +0.206639 |
| Signed interpolated density → negative-density floor 1e−20 | -0.116379 to -0.0425688 | -0.204943 to -0.0611562 | +0.108805 to +0.200282 | +0.185123 to +0.288514 |
| Legacy rectangular k/midpoint mu → Gauss–Legendre | -0.865845 to -0.630887 | -0.311236 to -0.223242 | +0.630581 to +0.900455 | +0.230276 to +0.336356 |
| Centre-volume approximation → integrated volume | +0.00640434 to +0.0103511 | +0.00640434 to +0.0103511 | -0.01035 to -0.00640393 | -0.01035 to -0.00640393 |
| Arithmetic mean/damping redshift → geometric centre | +0.0140598 to +0.218064 | -0.00633529 to +0.173564 | -0.218257 to -0.0204079 | -0.177606 to +0.00302324 |
| Resolution sigma=c/R → FWHM conversion | -0.311348 to -0.0845352 | -0.0674066 to -0.017429 | +0.120673 to +0.35425 | +0.0242243 to +0.0811432 |
| EdS power scaling → CAMB sigma8 power scaling | -0.658974 to +0.172453 | -0.834158 to +0.242151 | -0.264013 to +0.553597 | -0.365225 to +0.700949 |
| CAMB interpolated linear spectrum → supplied Vega template | -0.067954 to -0.0380399 | -0.0969604 to -0.0650112 | -0.00249627 to -0.000850602 | -0.00204499 to -0.000501736 |
| Undamped covariance power → damped template wiggle | +0.0215074 to +0.0462278 | +0.00812482 to +0.0242962 | -0.0491547 to -0.0229366 | -0.0264892 to -0.00928686 |
| Mixed forest–galaxy forest width → mean auto-width squares | -4.10543 to -2.37603 | -2.38308 to -0.954074 | +2.46025 to +4.30386 | +0.981388 to +2.47954 |
| Observed-power polynomial residual → template PK−PKSB wiggle | -0.411502 to -0.271058 | -0.426712 to -0.269178 | +0.101309 to +0.23927 | +0.0698021 to +0.193221 |
| Backward k derivative → inverse AP derivative with fixed observed response | +1.01977 to +1.17409 | +1.12319 to +1.30083 | -1.17585 to -1.11237 | -1.44869 to -1.29297 |
| Pair response/forest roundtrip → per-field response/exact angle | +4.91993e-09 to +5.89701e-08 | +2.10369e-09 to +2.77487e-08 | -5.79374e-08 to -5.74867e-09 | -2.87388e-08 to -2.43718e-09 |
| c=299800 → 299792.458 in velocity/noise/response, volume fixed | -0.00109667 to -0.000802554 | -0.00144256 to -0.00110598 | +0.000778393 to +0.00106994 | +0.00108805 to +0.00142113 |
| Weight convergence rtol 1e−4 → 1e−5 | +6.79456e-12 to +1.65477e-09 | +7.06102e-12 to +1.8354e-09 | -1.89284e-09 to -6.02851e-12 | -2.16367e-09 to -6.33937e-12 |

Per-bin values, individual-spectrum effects, exact tested interactions, rank, correlation and ellipse-area changes are retained in [attribution-summary.json](../.validation/s4/attribution-summary.json). Source definitions and unchanged assumptions are in [s4-profile-inventory.md](s4-profile-inventory.md).

Maximum endpoint error closure residual: 1.28664e-06 fraction. Maximum tested refinement change across joint and individual errors: 0.0107825%. Maximum joint Fisher Frobenius relative refinement: 0.000214949. These refinements cover the Gauss–Legendre endpoint and single-switch controls (including the accuracy-side rectangular-425 check), not all two-switch interactions or the literal rectangular Fourier estimator.

The interaction statistic is ln[sigma_AB sigma_0/(sigma_A sigma_B)], computed separately in each direction, bin, spectrum and component. The maximum absolute value across retained controls is 0.00413745. A nonzero value measures failure of multiplicative single-switch attribution; it is not a statistical significance.

| Bin | Accuracy/fixed radial (%) | Accuracy/fixed transverse (%) | Correlation change | Ellipse-area change (%) | Rect425 forward radial (%) | Rect425 forward transverse (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | -3.860897 | -1.32995 | +0.01737468 | -4.285552 | +0.0006355858 | +0.0006428048 |
| 2 | -3.838112 | -1.354774 | +0.01702293 | -4.306164 | -0.01887624 | -0.02144259 |
| 3 | -3.685243 | -1.331909 | +0.01666438 | -4.159099 | -0.01114055 | -0.01594434 |
| 4 | -2.963541 | -0.8795054 | +0.01602495 | -3.04646 | -0.002760111 | -0.002975794 |
| 5 | -2.669346 | -0.6664862 | +0.01595763 | -2.548718 | -0.001404837 | -0.001259183 |
| 6 | -3.822668 | -1.608909 | +0.01689114 | -4.548846 | -0.005822513 | -0.005091345 |

Any joint or individual rank change among single-switch controls: False.
Maximum absolute refined joint changes (radial, transverse), percent: [0.01048912404106872, 0.010782498426287201]. Corresponding individual maxima: [0.01004836387699104, 0.010217452297134333].

| Interaction | Tested bins | Joint radial log range | Joint transverse log range | Largest individual absolute log |
| --- | --- | ---: | ---: | ---: |
| forward-fourier+operator | 1, 2, 3, 4, 5, 6 | -5.72156e-05 to +0.000829163 | +0.000149193 to +0.00197042 | 0.00158946 |
| forward-growth+operator | 1 | +3.01697e-06 to +3.01697e-06 | -8.84354e-06 to -8.84354e-06 | 2.02325e-05 |
| forward-growth+peak | 1 | +0.00207986 to +0.00207986 | +0.00207796 to +0.00207796 | 0.00207836 |
| forward-magnitude+negative | 2, 3, 4, 5, 6 | -0.00148219 to -0.000487093 | -0.00205199 to -0.000656572 | 0.00413745 |
| forward-magnitude+operator | 6 | -8.17891e-08 to -8.17891e-08 | +6.60054e-08 to +6.60054e-08 | 2.28259e-06 |
| forward-magnitude+peak | 2, 3, 4, 5 | +1.53112e-08 to +2.97043e-07 | +1.63346e-08 to +3.54361e-07 | 9.66888e-06 |
| forward-negative+operator | 6 | -7.69197e-07 to -7.69197e-07 | +2.51289e-06 to +2.51289e-06 | 6.91136e-06 |
| forward-negative+peak | 2, 3, 4, 5 | +6.53527e-07 to +1.94877e-06 | +8.65993e-07 to +3.17514e-06 | 1.8659e-05 |
| forward-operator+peak | 1, 2, 3, 4, 5, 6 | +0.000526775 to +0.00137943 | +0.00151828 to +0.00281285 | 0.00272494 |
| reverse-fourier+operator | 1, 2, 3, 4, 5, 6 | -3.13337e-05 to +7.70985e-06 | -0.000146168 to -3.94787e-05 | 0.00182599 |
| reverse-growth+operator | 1 | +3.18425e-06 to +3.18425e-06 | -1.35318e-05 to -1.35318e-05 | 3.22613e-05 |
| reverse-growth+peak | 1 | -7.53231e-06 to -7.53231e-06 | -8.95528e-06 to -8.95528e-06 | 1.54499e-05 |
| reverse-magnitude+negative | 2, 3, 4, 5, 6 | -0.00142145 to -0.000478841 | -0.00201636 to -0.000647614 | 0.0041049 |
| reverse-magnitude+operator | 6 | +2.18422e-07 to +2.18422e-07 | +8.39151e-06 to +8.39151e-06 | 1.3604e-05 |
| reverse-magnitude+peak | 2, 3, 4, 5 | +1.46913e-06 to +2.39442e-06 | +1.56691e-06 to +2.45582e-06 | 1.77254e-05 |
| reverse-negative+operator | 6 | +5.50029e-07 to +5.50029e-07 | +1.20658e-05 to +1.20658e-05 | 1.87631e-05 |
| reverse-negative+peak | 2, 3, 4, 5 | +2.25407e-06 to +3.30522e-06 | +2.22582e-06 to +3.64123e-06 | 1.77954e-05 |
| reverse-operator+peak | 1, 2, 3, 4, 5, 6 | -0.000229246 to -5.99695e-05 | +0.000156309 to +0.000498731 | 0.00281239 |

The rectangular-425 reverse check restores a finer rectangular measure on accuracy: bin 1: -0.000605575%, -0.00116482%; bin 2: -0.00748741%, -0.00944257%; bin 3: -0.00573855%, -0.00874716%; bin 4: -0.00406226%, -0.00582926%; bin 5: -0.00581668%, -0.00773619%; bin 6: -0.0062395%, -0.00862527%.


Shared prescriptions have no independent parameter change here: early-lyaforecast recurrence; representative angular/velocity mode; density normalization and support; raw first-spacing density-width convention; SNR bright clamp, sentinel and floor; P1D prescription; independent-sampling noise; forecast selection; covariance-derivative exclusion. Pair versus field instrument response agrees for the equal forest pixel widths of these inputs. The field-response switch additionally removes the legacy forest-auto mu=k_parallel/(k+1e−10) roundtrip and therefore can retain a tiny nonzero effect.

The redshift switch changes the derivative mean and damping evaluation; covariance, source queries and geometry already use geometric redshift. The direct reconstruction uses the newly prepared CAMB redshift samples for galaxy f in both hybrid endpoints. The captured legacy forecast interpolated f between its arithmetic-bin samples. That small interpolation-versus-exact difference is not separately switched by redshift: it is included in the reported endpoint reconstruction residual, with no independently measured BAO contribution assigned. The covariance-damping hybrid adds the template-wiggle damping correction to the selected linear power. The peak switch adopts the template wiggle even when the covariance power source is legacy. These hybrid definitions must be retained when interpreting forward/reverse asymmetry. The constants switch changes remaining velocity, noise and response inputs with volume held fixed. The separate volume switch replaces saved centre volume with the physical integral and therefore already includes the c convention within volume. Refined rectangular integration, floor/width changes, fixed-comoving representative modes, McDonald and W12 are not new adopted profiles.
