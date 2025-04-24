% F=Fxi +Fyj + Fzk
% M=Mxi +Myj+ Mzk
% let F1,F2,F3,F4,F5,F6 be  the forces acting on the 6 links which are also
% in form of vectors 
%directions of Forces are

%F1= 156.6i - 35.31j -355.99k   0.40i -0.09j -0.91k left down front
%F2= -245.45i - 8.37j -330.1k     -0.59i -0.02j -0.80k left down back
%F3= 157.11i -51.85j -273.34k     0.49i -0.16j -0.85k left up front
%F4= -199.16i -64.96j -287.3k      -0.56i -0.18j -0.80k left up back
%F5= 325i - 36.33j                     0.99i-0.11 j     tie rod
%F6= 259.79i -270.53j               0.69i-0.72j  push rod
% wheel centre= 9.91i +2.50j
%assume
Fx=3572.183;
Fy=3087.04;
Fz=2534.42;
Mx=-66.46;
My=0;
Mz=-196.51;

syms F1 F2 F3 F4 F5 F6
eq1= F1*0.40 + F2*(-0.59) +F3*0.49 +F4*(-0.56) + F5*0.99 + 0.69*F6 ==Fx;
eq2= -(F1*0.09) -(F2*0.02) -(F3*0.16)-(F4*0.18) +F5*(-0.11)+ 0.72*F6== Fy;
eq3=  F1*(-0.91) -(F2*0.80) +F3*(-0.85) -(F4*0.80)==Fz;

eq4=  F1*(-0.91)*(0.124) + F2*(-0.80)*0.124 +F3*(-0.85)*0.244 +F4*(-0.80)*0.244 ==Mx;
eq5=  F1*(0.40)*(0.807) + F2*(-0.59)*(0.807) + F3*(-0.16)*(0.788)+ F4*(-0.22)*(0.788) + F4*0.71*(-0.007) +F5*0.99*0.15+F6*0.69*0==My;
eq6=  F1*0.29*0.577 +F2*(-0.33)*0.577 + F3*0.36*0.564 +F4*(-0.33)*0.564 +F5*(-0.11)*0.975+F6*(-0.72)*0.099==Mz;


sol=solve([eq1,eq2,eq3,eq4,eq5,eq6],[F1,F2,F3,F4,F5,F6]);

disp(sol.F1)
disp(sol.F2)