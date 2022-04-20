clear all;
close all;
clc;

cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts';



% add that folder plus all subfolders to the path.
addpath(genpath(pwd));

gpterm = 'x11';
sim_model; % can run it, if called directly, can't call it within this script...

mdim = 3;

preoptimise_dd;
gpterm = 'x11';
optimise_dd;



