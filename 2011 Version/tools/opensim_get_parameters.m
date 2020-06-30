function [joints, muscles, error] = opensim_get_parameters(osimfilename) 
% returns structure with joint parameters from osim file
% input: osimfilename
% outputs: joints, muscles and error=0 if successful, -1 if not
% fields of structure joints:
%   name
%   proximal_body
%   distal_body
%   limits(2)
%
% fields of structure muscles:
%   name
%   fmax
%   lceopt
%   lslack
%   pennopt
%   width
%   timeconstants(2)
%   origin
%   insertion
%   dofnames(num_dofs)
%   dofnums(num_dofs)
%   finger

% Dimitra Blana, June 2012
% (to replace get_info_osim that used Load_OSIM20)

import org.opensim.modeling.*
Mod = Model(osimfilename);
% initialize the system to get the initial state
state = Mod.initSystem;

Mod.equilibrateMuscles(state);

% define a context object which will allow us to access the state
context = OpenSimContext(state,Mod);

% determine the number of states in the model
nStates = Mod.getNumStates();

% get the coordinates structure
CoordSet = Mod.getCoordinateSet();
DofNames = ArrayStr();
CoordSet.getNames(DofNames);
nDofs = DofNames.getSize;

% loop through dofs
index=1;
for idof = 1:nDofs
    currentDof = CoordSet.get(idof-1);    
    dofname = char(currentDof.getName());
    joints{index}.limits(1) = currentDof.getRangeMin*180/pi;
    joints{index}.limits(2) = currentDof.getRangeMax*180/pi;
    
    dofs{index}=dofname;    
    currentJnt = currentDof.getJoint();    
    joints{index}.name = dofname;
    
    % proximal and distal body
    prox_body{index} = char(currentJnt.getParentBody());
    dist_body{index} = char(currentJnt.getBody());
    
    joints{index}.proximal_body = prox_body{index};
    joints{index}.distal_body = dist_body{index};
    
    index=index+1;    
end

% Get the set of bodies
BodySet = Mod.getBodySet();
%Count the bodies
nBodies = BodySet.getSize();

% loop through bodies 
for ibody = 2:nBodies
    currentBody = BodySet.get(ibody-1);    
    segments{ibody-1}=(char(currentBody.getName()));
end

% get the muscles
MuscleSet = Mod.getMuscles();
MuscleNames = ArrayStr();
MuscleSet.getNames(MuscleNames);
nMus = MuscleNames.getSize-12;

%%%%%%%%%%%%%%%%%%%%

% From Lisa Schutte's dissertation, appendices 2 and 3.
% "Using musculoskeletal models to explore strategies for improving
% performance in electrical stimulation-induced leg cycle ergometry"
% LM Schutte - 1992 - Stanford University
% k1 = activation1
% k2 = activation2
% t_act = 1/(k1+k2)
% t_deact = 1/k2

% From default Schutte1993Muscle in the original osim model:
% <!--Parameter used in time constant of ramping up of muscle force-->
% <activation1>       7.66700000 </activation1>
% <!--Parameter used in time constant of ramping up and ramping down of
%     muscle force-->
% <activation2>       1.45985400 </activation2>

k1 = 7.66700000;
k2 = 1.45985400; 
t_act = 1/(k1+k2);
t_deact = 1/k2;

%%%%%%%%%%%%%%%%%%%%%%

for imus = 1:nMus

    currentMuscle = MuscleSet.get(imus+11);
    
    muscles{imus}.name = char(currentMuscle.getName());
    muscles{imus}.fmax = currentMuscle.getMaxIsometricForce();
    muscles{imus}.lceopt = currentMuscle.getOptimalFiberLength(); 
    muscles{imus}.lslack = currentMuscle.getTendonSlackLength();
    muscles{imus}.pennopt = currentMuscle.getPennationAngleAtOptimalFiberLength()*180/pi;
    muscles{imus}.width = 0.56;
    muscles{imus}.timeconstants = [t_act   t_deact];

    muspath = currentMuscle.getGeometryPath();
    PtSet = muspath.getPathPointSet();
    PtNames = ArrayStr();
    PtSet.getNames(PtNames);
    nPts = PtNames.getSize;
    originPt = PtSet.get(0);
    muscles{imus}.origin = char(originPt.getBody());
    insertionPt = PtSet.get(nPts-1);
    muscles{imus}.insertion = char(insertionPt.getBody());

    proximal_segment = find(strcmp(muscles{imus}.origin,prox_body),1,'first');
    distal_segment = find(strcmp(muscles{imus}.insertion,dist_body),1,'last');
    if isempty(distal_segment)
        distal_segment = 0; 
    end
    if isempty(proximal_segment)
        proximal_segment = nDofs+1; 
    end
    
    if proximal_segment>distal_segment,
        helpseg = proximal_segment;
        proximal_segment = distal_segment+1;
        distal_segment = helpseg-1;
    end
    
    candidate_dofs = dofs(proximal_segment:distal_segment);

    muscles{imus}.dofnames = [];
    muscles{imus}.dofnums =[];
        
    for idof=1:length(candidate_dofs)
        muscles{imus}.dofnames = [muscles{imus}.dofnames candidate_dofs(idof)];
        muscles{imus}.dofnums =  [muscles{imus}.dofnums find(strcmp(candidate_dofs(idof),dofs))];
    end
        
end

error=0;