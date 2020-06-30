function [RMSlen1,RMSmom1,RMSlen2,RMSmom2] = compare_momarms
% compares three different sets of polynomials to Opensim
%
% dof (nodof): makes (does not make) sure muscles cross all degrees 
%              of freedom they are supposed to cross
% length (nolength): fits (does not fit) lengths as well as moment arms

osimfilename = '../das3.osim';
mydir = 'opensim_momarms';

% get joint and muscle information from osim file
[~, muscles, error] = opensim_get_parameters(osimfilename); 
if error==-1, return; end

% compares polynomials to original Opensim lengths/moment arms:
polydir = 'nolength_nodof';        
[RMSlen1,RMSmom1] = check_momarms(muscles,mydir,polydir);
% polydir = 'nolength_dof';
% check_momarms(muscles,mydir,polydir);
polydir = 'length_dof';
[RMSlen2,RMSmom2] = check_momarms(muscles,mydir,polydir);

for imus=1:length(muscles)
    fprintf('%s   length: %f m /%f m, moment arm: %f m /%f m\n',muscles{1,imus}.name,...
        RMSlen1(imus),RMSlen2(imus),RMSmom1(imus),RMSmom2(imus));
end

figure; subplot(2,1,1);
plot(RMSlen1-RMSlen2,'o');
subplot(2,1,2);
plot(RMSmom1-RMSmom2,'o');