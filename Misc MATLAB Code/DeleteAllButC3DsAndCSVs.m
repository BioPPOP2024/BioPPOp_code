


%mainfolder = 'I:\NeilGroins181017\DatasetsMTPc3d' ;
mainfolder = 'F:\02 LOWER LIMB PROCESSED DATA REPOSITORY\NFR' ;
filestodelete = dirPlus(mainfolder, 'FileFilter', '\.(history|mkr|mp|system|vsk|x1d|x2d|xcp|xml)$');
%filestodelete = dirPlus(mainfolder, 'FileFilter', '\.(history|mkr|mp|system|vsk|x1d|x2d|xcp|xml|avi|c3d|enf)$');
%filestodelete = dirPlus(mainfolder, 'FileFilter', '\.(Cut|HH|DL|DJ|Hop|Squat|CMJ1|CMJ2|CMJ3|CMJ4)$');
%filestodelete = dirPlus(mainfolder, 'FileFilter', '\.avi');
%filestodelete = dirPlus(mainfolder, 'FileFilter', '\.csv');
%filestodelete = dirPlus(mainfolder, 'FileFilter', 'SLDL');
%filestodelete = dirPlus(mainfolder, 'FileFilter', '(ROM|Static|DJ|Rebound|Cut|Squat)');
%filestodelete = dirPlus(mainfolder, 'FileFilter', 'Cut|HH|DL|DJ|Hop|Squat|CMJ1|CMJ2|CMJ3|CMJ4');

delete(filestodelete{:,:});