% calculating hand position hessians
function [hessian ] = das_hessian( x,model )
% x need only include the 11 angles [x(1:11)]
for i=1;
if nargin<2
    load model_struct
end
p1=model.joints{1, 1}.location';    
p2=model.joints{1, 2}.location';    
p3=model.joints{1, 3}.location';    
p4=model.joints{1, 4}.location';    
p5=model.joints{1, 5}.location';    
p6=model.joints{1, 6}.location';    
p7=model.joints{1, 7}.location';    
p8=model.joints{1, 8}.location';    
p9=model.joints{1, 9}.location';    
p10=model.joints{1, 10}.location';  
p11=model.joints{1, 11}.location';  
p12=model.joints{1, 12}.location';  
p13=model.joints{1, 13}.location'; 

theta=zeros(4,1);
r11=model.joints{1, 11}.r1_axis;
r12=model.joints{1, 12}.r1_axis;
theta(1)=atan2(r11(2),r11(3));
theta(2)=acos(r11(1));
theta(3)=atan2(r12(2),r12(3));
theta(4)=acos(r12(1));
end
% y z x y z x y z y 'x' 'y'
hessian=zeros(3,11,11);
for m=1:11;
for i=1:11;
    n=zeros(1,11);j=ones(10);
    j(1:i)=0;j(1:m)=0;n(i)=n(i)+1;n(m)=n(m)+1;
  hessian(:,i,m)=((p1+p2)*j(1)+r3(x(1),'y',n(1))*...
    (p3*j(2)+r3(x(2),'z',n(2))*...
    (p4*j(3)+r3(x(3),'x',n(3))*...
    (p5*j(4)+r3(x(4),'y',n(4))*...
    (p6*j(5)+r3(x(5),'z',n(5))*...
    (p7*j(6)+r3(x(6),'x',n(6))*...
    (p8*j(7)+r3(x(7),'y',n(7))*...
    (p9*j(8)+r3(x(8),'z',n(8))*...
    (p10*j(9)+r3(x(9),'y',n(9))*...
    (p11*j(10)+r3(-theta(1),'x')*r3(-theta(2),'y')*r3(x(10),'x',n(10))*r3(theta(2),'y')*r3(theta(1),'x')*...
    (p12*j(11)+r3(-theta(3),'x')*r3(-theta(4),'y')*r3(x(11),'x',n(11))*r3(theta(4),'y')*r3(theta(3),'x')*p13)))))))))));
end
end
end