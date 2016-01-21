### PCA-MS: Variational Multi-Phase Segmentation using High-Dimensional Local Features
A **variational multi-phase segmentation** framework based on the **Mumford-Shah** energy, combined with **PCA**-based dimension reduction is used to segment color or gray-value images into regions of different structure identified by **high-dimensional features**, such as **local spectral histograms (Texture)** or **localized Fourier transforms (Crystals)**.
___
**CONTENTS**: C++ **source code** reproducing texture and crystal segmentation results presented at the **IEEE Winter Conference on Applications of Computer Vision (WACV 2016)** - see reference below.

    |- quocmesh/                                        source code
    |+ - fini...jects/highDi...tation/applyTxtMerge.m   MATLAB script for TxtMerge post-processing
    |
    |- quocGCC/                                         compilation folder
    |+ - go.sh                                          bash script for CMake (see README.txt)
    |+ - finishedProjects/highDimFeatureSegmentation    executables (created during compilation)
    |
    |- LICENSE.txt                                      Common Development and Distribution License
    |- README.txt                                       intructions for compilation and execution
**INSTRUCTIONS**: Please refer to the accompanying **README.txt**
___
**DATA & RESULTS**: [http://nmevenkamp.github.io/pcams-data-wacv2016](http://nmevenkamp.github.io/pcams-data-wacv2016)
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
