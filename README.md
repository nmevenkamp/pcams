### PCA-MS: Variational Multi-Phase Segmentation using High-Dimensional Local Features
___
A **variational multi-phase segmentation** framework based on the **Mumford-Shah** energy, combined with **PCA**-based dimension reduction is used to segment color or gray-value images into regions of different structure identified by **high-dimensional features**, such as **local spectral histograms (Texture)** or **localized Fourier transforms (Crystals)**.
___
**CONTENTS**: C++ **source code** reproducing texture and crystal segmentation results presented at the **IEEE Winter Conference on Applications of Computer Vision (WACV 2016)** - see reference below.

    |- quocmesh: source code
    |+ - quocmesh/finishedProjects/highDimFeatureSegmentation/applyTxtMerge.m: MATLAB script to apply TxtMerge post-processing
    |- quocGCC: compilation folder
    |+ - quocGCC/go.sh: bash script for CMake (see README.txt)
    |+ - quocGCC/finishedProjects/highDimFeatureSegmentation: WACV 2016 executables; created during compilation
    |- LICENSE.txt: Common Development and Distribution License
    |- README.txt: instructions for compilation (GCC; Linux or MacOSX) and execution
**INSTRUCTIONS**: Please see the accompanying **README.txt**
___
**LICENSE**: PCA-MS is distributed under the terms of the [Common Development and Distribution License](LICENSE.txt).
___
**REFERENCE**:
Mevenkamp, N., and Berkels, B. _Variational Multi-Phase Segmentation using High-Dimensional Local Features_. Applications of Computer Vision (WACV), 2016 IEEE Winter Conference on, 2016, (accepted)
___
**CONTACT**:<br>
Niklas Mevenkamp<br>
Aachen Institute for Advanced Study in Computational Engineering Science<br>
RWTH Aachen University<br>
mevenkamp@aices.rwth-aachen.de<br>
http://www.aices.rwth-aachen.de/people/mevenkamp
