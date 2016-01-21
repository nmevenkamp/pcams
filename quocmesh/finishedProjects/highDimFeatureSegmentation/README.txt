------------------ Instructions for applyTxtMerge.m ------------------

usage: 			applyTxtMerge ( imgDir, segDir )
parameters: 	imgDir = path to directory containing the ICPR 2014 contest texture mosaics (tm*_*_*.png)
				segDir = path to directory containing raw segmentations (seg*_*_*.png)

Calling will create a sub-folder /TxtMerge inside segDir and store the post-processed segmentations there.


IMPORTANT NOTE: applyTxtMerge.m DEPENDS on the ICPR2014 version of the FSEG MATLAB code,
				which can be downloaded here: http://web.ornl.gov/~jiy/FSEG_contest.zip
				
				Make sure you extract all required files (see list below) from the FSEG_contest.zip
				to ../finishedProjects/highDimFeatureSegmentation
				AND mex the files AdjMxCrt0.c, RmSmRg.c, SHcomp.c, SHedge_1s.c (see ReadMe.txt inside FSEG_contest.zip)
				before calling the applyTxtMerge function!
				
List of dependencies:
	AdjMxCrt0.c
	gabor_fn.m
	RmSmRg.c
	SHcomp.c
	SHedge_1s.c
	subImg.m
	TxtMerge.m
	
	
	
------------------ Instructions for Run_Evaluation_Outex.cpp ------------------

usage:			Run_Evaluation_Outex <segmentationDirectory> <pathToGroundTruthSegments>
parameters:		segmentationDirectory = path to the directory containing all Outex segmentations (seg000.png - seg099.png)
				pathToGroundTruthSegments = path to the image containing the ground truth labels (gt.png)

Calling will create a .txt file inside the <segmentationDirectory> containing the results of the quantitative evaluation.


IMPORTANT NOTE: Run_Evaluation_Outex requires the Hungarian external library in order to match segments
				in the ground truth to those detected by the segmentation algorithm.
				
				To connect the source code with the Hungarian external library,
				add "-DBUILD_AND_USE_HUNGARIAN=1" to the CMake command.
				This should automatically download and link the library.
				
				In case this does not work, try downloading the Hungarian library manually from
				http://robotics.usc.edu/~lantao/codes/hungarian-v2.0.zip
				and place all of its .h and .cpp files inside "quocmesh/external/hungarian-source".
				Then, re-run CMake.