function minusdLdq = opensim_get_dLdq(Mod, angles, Mus, Dofs)
% function minusdLdq = opensim_get_dLdq(Mod, angles, Mus, Dofs)
%
% This function calculates -dL/dq of muscle "Mus" about dof set 
% "Dofs" at a given angle matrix "angles" of opensim model "Mod"
%
% "Angles" can be a vector (one hand position) or a matrix (one hand
% position per row)
% "Dofs" can a string (one dof) or a cell array (multiple dofs)
%
% Adapted from opensim_get_momarm.m 
% Dimitra Blana, March 2012
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
for idof = 1:nDofs
    DofNames_str{idof} = char(DofNames.getitem(idof-1));
end

% How many moment arms do we want?
if iscell(Dofs)
    num_request_dofs = length(Dofs);
else
    num_request_dofs = 1;
end

% find indeces of dofs for moment arms
if num_request_dofs==1
        dofindex = find(strcmp(Dofs,DofNames_str));
else
    for idof = 1:num_request_dofs
        dofindex(idof) = find(strcmp(Dofs{idof},DofNames_str));
    end
end

% get the muscles
MuscleSet = Mod.getMuscles();
MuscleNames = ArrayStr();
MuscleSet.getNames(MuscleNames);
nMus = MuscleNames.getSize;
for imus = 1:nMus
    MuscleNames_str{imus} = char(MuscleNames.getitem(imus-1));
end
musindex = find(strcmp(Mus,MuscleNames_str));

% angles matrix: one position per row
[nrows,ncols] = size(angles);

if ncols~=nDofs
    if nrows~=nDofs
        errordlg('Angle matrix not the right size','Input Error');
        minusdLdq = [];
        return;
    else
        angles = angles';
    end
end

% initialise matrix for -dL/dq output
minusdLdq = zeros(size(angles,1),num_request_dofs);

for istep = 1:size(angles,1)
    if ~mod(istep,50)
        disp(['Muscle ',char(MuscleSet.get(musindex-1).getName()), ' - step ',...
            num2str(istep),' of ', num2str(size(angles,1))]);
    end
    
    % set dof values for this step
    for idof = 1:nDofs
        currentDof = CoordSet.get(idof-1);    
        currentDof.setValue(state,angles(istep,idof),1);
    end
    
    for idof = 1:num_request_dofs
        currentDof = CoordSet.get(dofindex(idof)-1);    
        
        % change the value of the dof by +-0.0001 and find length:
        currentDof.setValue(state,angles(istep,dofindex(idof))-0.0001,1);
        L1 = MuscleSet.get(musindex-1).getLength(state);
                        
        currentDof.setValue(state,angles(istep,dofindex(idof))+0.0001,1);                
        L2 = MuscleSet.get(musindex-1).getLength(state);
        
        minusdLdq(istep,idof) = -(L2-L1)/0.0002;
        
        % set dof to original value
        currentDof.setValue(state,angles(istep,dofindex(idof)),1);                
    end
end

