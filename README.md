### PCA-MS: Variational Multi-Phase Segmentation using High-Dimensional Local Features
___
A **variational multi-phase segmentation** framework based on the **Mumford-Shah** energy, combined with **PCA**-based dimension reduction is used to segment color or gray-value images into regions of different structure identified by **high-dimensional features**, such as **local spectral histograms (Texture)** or **localized Fourier transforms (Crystals)**.
___
**CONTENTS**: This package contains the source code used to produce texture and crystal segmentation results published at WACV 2016 (see reference below).
|- quocmesh (source code)
|+ - quocmesh/finishedProjects/highDimFeatureSegmentation (WACV 2016 executables)
|- quocGCC (compilation folder)
|+ - quocGCC/go.sh (adjust compiler path and call to execute CMake with required options; afterwards compile with make)
|- LICENSE.txt (Common Development and Distribution License)
|- README.md (this README)
|- data
|+ - data/Crystals (ground truth crystal images in Figure 2)
|+ - data/ICPR2014 (ground truth segments and texture mosaics from the ICPR 2014 contest; used for the results in Table 2)
|+ - data/Outex (ground truth segments and texture mosaics from the Outex_US_00000 test suite converted from .ras to .png; used for the results in Table 1)
|- results (same structure as data; contains the segmentations produced by this code)
|+ - results/ICPR2014/Raw (contains the segmentations before TxtMerge post-processing)
|+ - results/ICPR2014/TxtMerge (contains the segmentations after TxtMerge post-processing (use quocGCC/finishedProjects/highDimFeatureSegmentation/applyTxtMerge.m)

___
**LICENSE**: PCA-MS is distributed under the terms of the [Common Development and Distribution License](LICENSE.txt).
___
In case you want to refer to PCA-MS, please use the following citation:

Mevenkamp, N., and Berkels, B. **Variational Multi-Phase Segmentation using High-Dimensional Local Features**. Applications of Computer Vision (WACV), 2016 IEEE Winter Conference on, 2016, (accepted)
