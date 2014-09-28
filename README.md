SegWorm - Matlab Classes
========================

SegWorm is Matlab code from Dr. Eviatar Yemini built as part of the [WormBehavior database](http://wormbehavior.mrc-lmb.cam.ac.uk/)

This version of the code is an attempt to rewrite/clean up the code by using classes and packages to better encapsulate functionality.

As of Sep 2014, the idea is for this repo to be completely ported over to [movement_validation](https://github.com/openworm/movement_validation), in Python, soon.

Project Status
==============
- video to normalized worm - unfinished, current work at [`seg_worm.parseWormFromVideo`](https://github.com/JimHokanson/SegwormMatlabClasses/blob/master/%2Bseg_worm/parseWormFromVideo.m)
- normalized worm to features - finished, see example of usage in [`seg_worm.testing.features.t_001__mrcCodeVsNewCode`](https://github.com/JimHokanson/SegwormMatlabClasses/blob/master/%2Bseg_worm/%2Btesting/%2Bfeatures/t_001__mrcCodeVsNewCode.m)
- features to stats - finished - see example usage in [`seg_worm.testing.stats.t001_oldVsNewStats`](https://github.com/JimHokanson/SegwormMatlabClasses/blob/master/%2Bseg_worm/%2Btesting/%2Bstats/t001_oldVsNewStats.m)

Test files
==========
Files for features, somewhere in the [Google Drive folder example_data](https://drive.google.com/folderview?id=0B7to9gBdZEyGNWtWUElWVzVxc0E&usp=sharing)

[Jim & Michael need to work this out ... (i.e. issue 80 of the movement_validation repo)](https://github.com/openworm/movement_validation/issues/80)

Files for stats are all on the MRC ftp, see [@slarson's discussion of the how to connect to the FTP](https://github.com/openworm/OpenWorm/issues/82)

Setting up code
===============
1. Get the files and change any constants (these should be obvious in the files)
2. Add the repository folder to the Matlab path
3. If comparing to old code, add the folders oldFeatures and oldStats to the path


