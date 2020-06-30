function osim2bio
% Gets moment arms and lenghts from the Opensim DAS3 model and generates polynomials
% that replace the ones in the .bio file
%
% The original bio file das3.bio (in the main das3implicit folder) is not
% changed; the new .bio file is saved in tools/<mydir>

osimfilename = '../das3.osim';
mydir = 'opensim_momarms';

% get joint and muscle information from osim file
[joints, muscles, error] = opensim_get_parameters(osimfilename); 
if error==-1, return; end

if ~exist(mydir), mkdir(mydir); end

% build .mat files that contain all moment arms of the DAS3 model
error = get_momentarms(joints,muscles,osimfilename,mydir);
if error==-1, return; end

% build .mat files that contain all muscle lengths of the DAS3 model
error = get_lengths(joints,muscles,osimfilename,mydir);
if error==-1, return; end

% compute best fitting polynomials for muscle-tendon lengths as a function
% of kinematic degrees of freedom q
error = musclepath_poly(muscles,mydir);
if error==-1, return; end

% compares polynomials to original Opensim lengths/moment arms
check_momarms(muscles,mydir);

% update the bio file with new moment arm coefficients
change_momarms_bio(muscles,mydir);