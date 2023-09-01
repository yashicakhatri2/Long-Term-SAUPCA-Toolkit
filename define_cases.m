function [COE2, P_prop] = define_cases(case_flag)
% Defined test cases that have pre-saved results to confirm functionality.
% =======================================================================
% INPUTS = 
% case_flag: Test number from define_cases function 
% 
% OUTPUTS = 
% COE2: Classical orbital element set for object 2 in conjunction
% P_prop: Nominal Time of Closest Approach for the two object in conjunction
% =======================================================================

if case_flag == 1
    COE2 = [42096.4926672584        0.0976248366420683         0.784134429378206     -0.000121979064338973       -0.0174246915204144         0.053500575481622]';
    P_prop = 1.5 * 24 * 3600;
end

end