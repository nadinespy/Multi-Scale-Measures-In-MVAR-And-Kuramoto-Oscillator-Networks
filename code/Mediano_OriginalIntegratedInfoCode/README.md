Measuring integrated information
================================

This zip file contains code necessary to reproduce the key results in

Mediano, P.A.; Seth, A.K.; Barrett, A.B. *Measuring Integrated Information:
Comparison of Candidate Measures in Theory and Simulation*. Entropy 2019, 21,
17.

The two provided Matlab functions, `TwoNodeExample` and `PlotNetworkSummary`
reproduce figures 1 and 5 in the paper. The rest of the figures can be easily
computed with additional parameter sweeps using the same base code. See the
docstrings and code of both functions for more information and usage of the
estimators.


## Minimal example

In addition to computing the (analytical) results in the paper, this bundle
contains code to compute integration measures in data. This is a minimal
example to compute whole-minus-sum integrated information in a time series of
two random Gaussian variables:

```octave
X = randn(100, 2);
phiCalc = infodynamics.measures.continuous.gaussian.IntegratedInformationCalculatorGaussian();
phiCalc.compute(X)
```

The integration calculators used here (except for PhiG) are part of a private
fork of Lizier's JIDT toolbox, available in `infodynamics.jar`. Further
information on how to use JIDT calculators can be found in JIDT's website and
in the two functions provided.


## Relevant links

* **Original paper**:
    https://www.mdpi.com/1099-4300/21/1/17

* **JIDT toolbox**:
    https://www.github.com/jlizier/jidt

* **Phi toolbox** (for PhiG only):
    https://figshare.com/articles/code/phi_toolbox_zip/3203326/5


Pedro Mediano, 2020

