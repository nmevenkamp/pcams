---------------- Instructions for using the applyTxtMerge MATLAB function ------------------

usage: 			applyTxtMerge ( imgDir, segDir )
parameters: 	imgDir = path to directory containing the ICPR 2014 contest texture mosaics (tm*_*_*.png)
				segDir = path to directory containing raw segmentations (seg*_*_*.png)

Calling will create a sub-folder /TxtMerge inside segDir and store the post-processed segmentations there.


IMPORTANT NOTE: the applyTxtMerge.m MATLAB function provided here
				DEPENDS on the ICPR2014 version of the FSEG MATLAB code,
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