PCA-MS: Variational Multi-Phase Segmentation using High-Dimensional Local Features

A variational multi-phase segmentation framework based on the Mumford-Shah energy, combined with PCA-based dimension reduction is used to segment color or gray-value images into regions of different structure identified by high-dimensional features, such as local spectral histograms (Texture) or localized Fourier transforms (Crystals).

LICENSE: PCA-MS is distributed under the terms of the Common Development and Distribution License (see LICENSE.txt).

COMPILE INSTRUCTIONS (for GCC under Linux or MacOSX)
1) go to the quocGCC folder and open the go.sh with a text editor.
2) replace the paths in "export CC=..." and "export CXX=..." with the according paths to gcc and g++ on your machine
3) open a terminal and cd to quocGCC
4) execute "sh go.sh"
5) execute "make"

Afterwards, inside the terminal, you can execute either of the following

Run_Clustering_Outex
Run_PCAMS_Crystals
Run_PCAMS_ICPR2014
Run_PCAMS_Outex

inside "quocGCC/finishedProjects/highDimFeatureSegmentation".
The usage is: Run_... <sourceDirectory> <outputDirectory>

In order to perform the TxtMerge post-processing used in Table 2, use
"quocmesh/finishedProjects/applyTxtMerge.m".
Please regard the required setup procedure detailed in the accompanying README file.

The data and results used in the WACV 2016 paper can be downloaded here:


If you have any further questions regarding the code or the algorithm,
please feel free to contact me: mevenkamp@aices.rwth-aachen.de


In case you want to refer to PCA-MS, please use the following citation:

Mevenkamp, N., and Berkels, B. Variational Multi-Phase Segmentation using High-Dimensional Local Features. Applications of Computer Vision (WACV), 2016 IEEE Winter Conference on, 2016, (accepted)