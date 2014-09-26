SegWorm - Matlab Classes
========================

SegWorm is Matlab code from Dr. Eviatar Yemini built as part of the WormBehavior database (http://wormbehavior.mrc-lmb.cam.ac.uk/)

This version of the code is an attempt to rewrite/clean up the code by using classes and packages to better encapsulate functionality.

Project Status
==============
video to normalized worm - unfinished, current work at seg_worm.parseWormFromVideo
normalized worm to features - finished, see example of usage in seg_worm.testing.features.t_001__mrcCodeVsNewCode
features to stats - finished - see example usage in seg_worm.testing.stats.t001_oldVsNewStats

Test files
==========
Files for features, somewhere in:
https://drive.google.com/folderview?id=0B7to9gBdZEyGNWtWUElWVzVxc0E&usp=sharing

Jim & Michael need to work this out ...

Files for stats are all on the MRC ftp, see ...

Setting up code
===============
1. Get the files and change any constants (these should be obvious in the files)
2. Add the repository folder to the Matlab path
3. If comparing to old code, add the folders oldFeatures and oldStats to the path


