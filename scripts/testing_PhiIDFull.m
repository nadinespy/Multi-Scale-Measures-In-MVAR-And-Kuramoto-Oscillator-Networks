cd '/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab';
addpath('/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts');
javaaddpath('/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/infodynamics.jar');
            
eval('pkg load statistics');

redundancy_function = 'mmi';
tau = 1;
data = rand(2, 1000);

phiid = PhiIDFull(data, tau, redundancy_function)

pkg install -forge statistics