function length = opensim_get_length(Mod, angles)
% function length = opensim_get_length(Mod, angles)
%
% This function calculates the length of all muscles at a given angle 
% matrix "angles" of opensim model "Mod"
%
% "Angles" can be a vector (one hand position) or a matrix (one hand
% position per row)
%
% Adapted from opensim_muscle_analysis.m by Glen Lichtwark
% 
% Dimitra Blana, February 2012
%
% 28/3/2012: Use setValue (with active constraints) instead of set state

import org.opensim.modeling.*
% initialize the system to get the initial state
state = Mod.initSystem;

Mod.equilibrateMuscles(state);

% get the coordinates structure
CoordSet = Mod.getCoordinateSet();
DofNames = ArrayStr();
CoordSet.getNames(DofNames);
nDofs = DofNames.getSize;

% get the muscles
MuscleSet = Mod.getMuscles();
MuscleNames = ArrayStr();
MuscleSet.getNames(MuscleNames);
nMus = MuscleNames.getSize;

% angles matrix: one position per row
[nrows,ncols] = size(angles);

if ncols~=nDofs
    if nrows~=nDofs
        errordlg('Angle matrix not the right size','Input Error');
        length = [];
        return;
    else
        angles = angles';
    end
end

% initialise vector for length output
length = zeros(size(angles,1),nMus);

for istep = 1:size(angles,1)    
    for idof = 1:nDofs
        currentDof = CoordSet.get(idof-1);    
        currentDof.setValue(state,angles(istep,idof),1);
    end

    for imus = 1:nMus
        length(istep,imus) = MuscleSet.get(imus-1).getLength(state);
    end
end
