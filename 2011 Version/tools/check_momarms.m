function [RMSlen,RMSmom] = check_momarms(muscles,mydir,polydir)
% compares polynomials to original Opensim lengths/moment arms

if nargin>2
    momarm_coef = load([mydir '/' polydir '/poly_results']);
else
    momarm_coef = load([mydir '/poly_results']);
end

RMSlen = zeros(length(muscles),1);
RMSmom = zeros(length(muscles),1);
nDofs = zeros(length(muscles),1);
nDofs_poly = zeros(length(muscles),1);

for imus = 1:length(muscles)
    
    onemus = muscles{1,imus};
    fprintf('(%i) %s...\n',imus,onemus.name);

    % get moment arm out of .mat file created by get_momentarms.m
    musfilename = [mydir,'\momarms_',muscles{imus}.name,'.mat'];
    ma = load(musfilename);
    % get lengths out of .mat file created by get_lengths.m
    musfilename = [mydir,'\lengths_',muscles{imus}.name,'.mat'];
    lg = load(musfilename);
    
    numframes = size(lg.alljnts,1);
        
    num_musdofs = length(onemus.dofnums);
    nDofs(imus) = num_musdofs;
    musdofs = onemus.dofnums;    
    mus = momarm_coef.mus_model{1,imus};
    
    nDofs_poly(imus) = sum(sum(mus.lparams)~=0);
    
    allframes = 1:1:numframes;
    poly_len=zeros(length(allframes),1);
    poly_momarms=zeros(length(allframes),num_musdofs);
    osim_len=zeros(length(allframes),1);
    osim_momarms=zeros(length(allframes),num_musdofs);
    index=1;
    for iframe=allframes
        q = (ma.alljnts(iframe,musdofs) + 1e-6);	% protect against angle = 0.0
        poly_length = 0;
        poly_ma = zeros(1,num_musdofs);
        for iparam = 1:mus.num_lparams
            term = mus.lcoef(iparam);
            for idof=1:num_musdofs
                term = term*q(idof)^mus.lparams(iparam,idof);
            end
            poly_length = poly_length+term;        
            for idof=1:num_musdofs
                poly_ma(idof) = poly_ma(idof)-mus.lparams(iparam,idof)*term/q(idof);
            end
        end
        poly_len(index)=poly_length;
        poly_momarms(index,:)=poly_ma;
        osim_len(index)= lg.alllengths(iframe);
        osim_momarms(index,:)= ma.allmomarms(iframe,:);
        index=index+1;
    end
    RMSlen(imus) = (sqrt(mean((poly_len-osim_len).^2)));       
%    fprintf('Newpoly-Opensim length RMS error: %f m\n\n',RMSlen(imus));

    for i=1:num_musdofs
        RMSfull = (sqrt(mean((poly_momarms(:,i)-osim_momarms(:,i)).^2)));       
%        fprintf('Newpoly-Opensim %s moment arm RMS error:   %f m\n\n',onemus.dofnames{i},RMSfull);
    end

    RMSmom(imus) = (sqrt(mean(mean((poly_momarms-osim_momarms).^2))));       
%    fprintf('Newpoly-Opensim mean moment arm RMS error:   %f m\n',RMSmom(imus));
    
end

figure;
subplot('Position',[0.1 0.4 0.8 0.5]); 
plot(RMSlen,'bo'); hold on; plot(RMSmom,'r^'); title('Polynomial-Opensim RMS error (m)') 
legend('length','moment arm'); 
subplot('Position',[0.1 0.1 0.8 0.25]);
plot(nDofs,'b*'); hold on; plot(nDofs_poly,'r*'); legend('#dofs Opensim','#dofs polynomial');
xlabel('Muscle number');

