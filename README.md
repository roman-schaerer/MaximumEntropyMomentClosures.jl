# MaximumEntropyMomentClosures

Methods to evaluate maximum entropy distributions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://roman-schaerer.github.io/MaximumEntropyMomentClosures.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://roman-schaerer.github.io/MaximumEntropyMomentClosures.jl/dev/)
[![Build Status](https://github.com/roman-schaerer/MaximumEntropyMomentClosures.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/roman-schaerer/MaximumEntropyMomentClosures.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/roman-schaerer/MaximumEntropyMomentClosures.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/roman-schaerer/MaximumEntropyMomentClosures.jl)


**Warning**: Currently, this package is at an early stage and considered experimental. The interfaces are not yet stable and are expected to change.

For now this package provides routines for the evaluation of one-dimensional maximum entropy distributions over an interval.

Roadmap:
- Stabilize interface
- Improve robustness of Lagrange parameter evaluation
- Improve run-time efficiency by using blocked evaluation of moments
- Extend to multidimensional distributions
