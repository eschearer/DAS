function act  = computeActivations(M,tau,alpha0)


% solve Malpha=tau for alpha

    x=alpha0;
    Bk=eye(9);
    func=@penaltyFunc;
    [f,g]=func(x,M,tau);
    
    it=0;
    while norm(g)>1e-3 && it<120
        p=-Bk*g;
        alpha2=armijo(func,x,f,g,p,M,tau);
        xold=x;
        x=x+alpha2*p;
        sk=x-xold;
        gold=g;
        [f,g]=func(x,M,tau);
        yk=g-gold;
        rho=1/(yk'*sk);
        Bk=(eye(9)-rho*sk*yk')*Bk*(eye(9)-rho*yk*sk')+rho*sk*sk';
        it=it+1;
    end
    for j=1:length(x)
        if x(j)<0
            x(j)=0;
        end
    end
    if max(x)>1
        x=x/max(x);
    end
    
    act=x;
end

function alpha1 = armijo(func,xk,fk,gradk,pk,M,tau)
% Originally from Dr. Schearer. Modified by Derek Wolf
alpha1   = 1;
c       = 0;
rho     = 0.5;

n=0;
while fk + c*alpha1*gradk'*pk < func(xk+alpha1*pk,M,tau) % sufficient descent condition is not met
    alpha1 = rho*alpha1;
    n=n+1;
end
end

function [f,g]=penaltyFunc(x,M,tau)

c=1;
c2=5000;

% Define function
K=0;
for i=1:length(x)
    if x(i)<0
        K=K+x(i)^2;
    elseif x(i)>1
        %K=K+(x(i)-1)^2;
        %else add 0 to K if 0<=x(i)<=1
    end
end

f=norm(x)^2+c*norm(M*x-tau)^2+c2*K;


if nargout>1
    % Define gradient
    K=zeros(9,1);
    for i=1:length(x)
        if x(i)<0
            K(i)=2*x(i);
        elseif x(i)>1
            %K(i)=2*(x(i)-1);
        end
    end
    
    g=2*x+2*c*M'*(M*x-tau)+c2*K;
end

end



