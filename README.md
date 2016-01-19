### PCA-MS: Variational Multi-Phase Segmentation using High-Dimensional Local Features
___
A **variational multi-phase segmentation** framework based on the **Mumford-Shah** energy, combined with **PCA**-based dimension reduction is used to segment color or gray-value images into regions of different structure identified by **high-dimensional features**, such as **local spectral histograms (Texture)** or **localized Fourier transforms (Crystals)**.
___
**CONTENTS**: This package contains the C++ **source code** used to produce texture and crystal segmentation results published at WACV 2016 (see reference below).

    |- quocmesh: source code
    |+ - quocmesh/finishedProjects/highDimFeatureSegmentation: WACV 2016 executables
    |- quocGCC: compilation folder
    |+ - quocGCC/go.sh: bash script for CMake (see README.txt)
    |- LICENSE.txt: Common Development and Distribution License
    |- README.txt: instructions for compilation (GCC, Linux or MacOSX) and execution
    |- data
    |+ - data/Crystals: ground truth crystal images (Figure 2)
    |+ - data/ICPR2014: ground truth segments and texture mosaics from the ICPR 2014 contest (Table 2)
    |+ - data/Outex: ground truth segments and texture mosaics from the Outex_US_00000 test suite converted from .ras to .png (Table 1)
    |- results: same structure as data; contains segmentations produced by this code
    |+ - results/ICPR2014/Raw: segmentations before TxtMerge post-processing
    |+ - results/ICPR2014/TxtMerge: segmentations after TxtMerge post-processing (using quocGCC/finishedProjects/highDimFeatureSegmentation/applyTxtMerge.m)
___
**INSTRUCTIONS**: Please see the accompanying README.txt
___
**LICENSE**: PCA-MS is distributed under the terms of the [Common Development and Distribution License](LICENSE.txt).
___
**CONTACT**: mevenkamp@aices.rwth-aachen.de
___
**REFERENCE**:
Mevenkamp, N., and Berkels, B. _Variational Multi-Phase Segmentation using High-Dimensional Local Features_. Applications of Computer Vision (WACV), 2016 IEEE Winter Conference on, 2016, (accepted)
