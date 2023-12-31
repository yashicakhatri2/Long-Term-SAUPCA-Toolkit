function [OUT1, OUT2, OUT3] = savedResults(case_flag, type)
% This function contains saves results for some test cases. These can be
% retrieved from the main function as needed.
%
% =======================================================================
% INPUTS = 
% case_flag: Test number from define_cases function 
% type: 1 - Monte Carlo results of probability of collision, 0 - Semi-analytical results of probability of collision
% 
% VARIABLES = 
% MC_Pc: Monte Carlo probability of collision
% points: Number of points used in computation of MC_Pc
% GMMSTT_Pc: Semi-analytical probability of collision
% JMAX: Number of GMM components used in semi-analytical computation 
% MCRunTime: Runtime of the Monte Carlo calculations
% GMMSTTRunTime: Runtime of the semi-analytial method calculations
% 
% OUTPUTS = 
% OUT1: Probability of collision
% OUT2: Number of points or number of GMM components used
% OUT3: Time taken for calculations
% =======================================================================

if case_flag == 1
    MC_Pc = 0.0001325; 
    points = 1E7;
    GMMSTT_Pc = [0.00013101284684046];
    JMAX = [101];
    MCRunTime = NaN; % Use tic toc output
    GMMSTTRunTime = 3033;
end

if type == 1
    OUT1 = MC_Pc;
    OUT2 = points;
    OUT3 = MCRunTime;
else
    OUT1 = GMMSTT_Pc;
    OUT2 = JMAX;
    OUT3 = GMMSTTRunTime;
end

end