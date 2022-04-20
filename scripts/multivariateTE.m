% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/media/nadinespy/NewVolume/my_stuff/work/PhD/my_projects/EmergenceComplexityMeasuresComparison/EmergenceComplexityMeasuresComparison_Matlab/scripts/infodynamics.jar');


% infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate: for continuous multivariate variables
% TransferEntropyCalculatorMultiVariateGaussian: for Gaussian variables
% TransferEntropyCalculatorMultiVariateKraskov: for non-Gaussian variables
% --> calculates I(Y1_t, Y2_t, ... , Yn_t ; X1_t-1, X2_t-1, ... , Xn_t-1 | Y1_t-1, Y2_t-1, ... , Yn_t-1). i.e., 
% pairwise TE between the whole multivariates

% Change location of jar to match yours:
javaaddpath('../../infodynamics.jar');

% Generate some random binary data.
% Note that we need the *1 to make this a number not a Boolean,
%  otherwise this will not work (as it cannot match the method signature)
% numObservations = 100;
sourceArray = phase(1:end-1,:);

% sourceArray2 = (rand(numObservations,2)>0.5)*1;
sourceArray_random = sourceArray(randperm(numel(sourceArray)));
sourceArray_random_matrix = reshape(sourceArray_random, [1999,256]);

% Destination variable takes a copy of the first bit of the source in bit 1,
%  and an XOR of the two bits of the source in bit 2:
%destArray = [0, 0; sourceArray(1:numObservations-1, 1), xor(sourceArray(1:numObservations-1, 1), sourceArray(1:numObservations-1, 2))];
destArray = phase(2:end,:);

sourceDim = size(sourceArray, 2);
destDim = size(destArray, 2);

% Class TransferEntropyCalculatorMultiVariateKraskov: http://lizier.me/joseph/software/jidt/javadocs/v1.5/
% "Specifically, this class implements the pairwise or apparent transfer entropy"
teCalc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov');
teCalc.setProperty('k', '4');
teCalc.initialise(1, sourceDim, destDim);

% teCalc.initialise(1);	% Joe: "The initialise part is important, so the the estimator knows how many dimensions 
				% you will have for source and target. In the first example, it will assume 1 dimension only, 
				% so the calculation it will do is basically univariate (that's why the answers are different). 
								
teCalc.setObservations(octaveToJavaDoubleMatrix(sourceArray), ... % The octaveToJavaDoubleMatrix conversion is less 
			octaveToJavaDoubleMatrix(destArray));		% important - if you're in matlab that actually 
											% does nothing, it's only important if you're 
											% running octave instead.
result1_nats = teCalc.computeAverageLocalOfObservations();
result1_bits = result1_nats/(1/log(2));

% Class TransferEntropyCalculatorMultiVariateGaussian: http://lizier.me/joseph/software/jidt/javadocs/v1.5/:
teCalc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorMultiVariateGaussian');
teCalc.setProperty('k_HISTORY', '3');
teCalc.setProperty('k_TAU', '1');
teCalc.initialise();
teCalc.setObservations(sourceArray, destArray);
result3_nats = teCalc.computeAverageLocalOfObservations(); % result is given in nats! 1 nat/(1/ln(2)) = 1 nat/1.4426950408889 = 1 bit
result3_bits = result3_nats/(1/log(2));

% Class TransferEntropyCalculatorDiscrete: http://lizier.me/joseph/software/jidt/javadocs/v1.5/
% Joe Lizier: "For discrete, the key is that we convert any of the Y_t , X_t-1 or Y_t-1 which are multi-
% variate into a univariate time series, by basically converting into a larger alphabet space, e.g. 
% 2x binary variables convert into a base-4 univariate. This is done simply by using the utility 
% computeCombinedValues method, and setting the base or alphabet appropriately in the constructor, 
% and then you just use the calculator as is."

% "Specifically, this class implements the pairwise or apparent transfer entropy; i.e. we compute the 
% transfer that appears to come from a single source variable, without examining any other potential sources"

numObservations = 100;
sourceArray=(rand(numObservations,2)>0.5)*1;
sourceArray2=(rand(numObservations,2)>0.5)*1;
% Destination variable takes a copy of the first bit of the source in bit 1,
%  and an XOR of the two bits of the source in bit 2:
destArray = [0, 0; sourceArray(1:numObservations-1, 1), xor(sourceArray(1:numObservations-1, 1), sourceArray(1:numObservations-1, 2))];
% Create a TE calculator and run it:
teCalc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 4, 1); % 4: base - number of symbols for each variable. 
															   % E.g. binary variables are in base-2.
															   % 1: destHistoryEmbedLength - embedded history 
															   % length of the destination to condition on
teCalc.initialise();
% We need to construct the joint values of the dest and source before we pass them in,
% and need to use the matrix conversion routine when calling from Matlab/Octave:
mUtils = javaObject('infodynamics.utils.MatrixUtils');
teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(sourceArray), 2), ...
		mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(destArray), 2));
fprintf('For source which the 2 bits are determined from, result should be close to 2 bits : ');
result = teCalc.computeAverageLocalOfObservations()
teCalc.initialise();
teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(sourceArray2), 2), ...
		mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(destArray), 2));
fprintf('For random source, result should be close to 0 bits in theory: ');
result2 = teCalc.computeAverageLocalOfObservations()
fprintf('\nThe result for random source is inflated towards 0.3 due to finite observation length (%d).', ...
	'One can verify that the answer is consistent with that from a random source by checking:', ...
	'teCalc.computeSignificance(1000); ans.pValue\n', teCalc.getNumObservations());