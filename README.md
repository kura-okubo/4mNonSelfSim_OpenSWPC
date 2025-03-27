# Extended version of OpenSWPC for laboratory waveform propagation modeling

This repository includes an extended version of the FDM-based software [**OpenSWPC**](https://github.com/OpenSWPC/OpenSWPC/tree/master) for waveform modeling on laboratory rock specimens.

## Contents

1. Minimum working example for ball-drop source: [**example_balldrop**](example_balldrop).

2. Cross-verification of the extended software: [**cross-verification**](cross-verification).

>This extended software was developed for the analysis of [non-self-similar laboratory earthquakes](https://github.com/kura-okubo/4mNonSelfSim_Paper).


> **Reference:**
> Okubo, K., Yamashita, F., & Fukuyama, E. (2025). Dynamics of non-self-similar earthquakes illuminated by a controlled fault asperity, submitted.


## Installation
Installation instructions can be found [here](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/develop/example_balldrop/README.md) as well as in [the original documentation](https://openswpc.github.io/1._SetUp/0100_trial/).

## Newly Developed Modules

1. Velocity model for the rock specimen ([code](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/develop/src/swpc_3d/m_vmodel_balldropseg_sidecoord.F90)).
2. Ball-drop impact as a single-force source ([code](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/a4444600be3d6318acf16f3943c39b058d8b5268/src/shared/m_fdtool.F90#L542)).
3. Velocity model output (implemented for debugging; available only in single-core mode) ([code](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/a4444600be3d6318acf16f3943c39b058d8b5268/src/swpc_3d/m_medout.F90#L1)).
4. Extended free-surface detection algorithm:
   - [Detection implementation](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/a4444600be3d6318acf16f3943c39b058d8b5268/src/swpc_3d/m_medium.F90#L549)
   - [Reduction to second-order accuracy](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/a4444600be3d6318acf16f3943c39b058d8b5268/src/swpc_3d/m_kernel.F90#L170)
5. Updated makefile ([code](https://github.com/kura-okubo/4mNonSelfSim_OpenSWPC/blob/a4444600be3d6318acf16f3943c39b058d8b5268/src/swpc_3d/makefile#L44)).


## Software Version

This extended version was forked from `v5.1.0`. See the latest software updates in the [original repository](https://github.com/OpenSWPC/OpenSWPC).

## Original license
```
The MIT License (MIT)

Copyright 2013-2025 Takuto Maeda

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Accompanying Paper

Maeda, T., S. Takemura, and T. Furumura (2017),
OpenSWPC: An open-source integrated parallel simulation code for modeling seismic wave propagation in 3D heterogeneous viscoelastic media,
_Earth Planets Space_, 69, 102.
doi:[10.1186/s40623-017-0687-2](https://doi.org/10.1186/s40623-017-0687-2)