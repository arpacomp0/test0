function [ mismatch ] = ffeq(x_arg, ...
                   AA,TE_whole_0,TJ_whole_0, ...
                   TE_whole_1,TJ_whole_1,genlist, ...
                   loadlist,snet,droop_gain)
%
%
%
gen_set=x_arg(1:3);
Vm_0=x_arg(5:13);
theta_0=x_arg(14:21);
Vbus_0=Vm_0.*exp(j*[0;theta_0]);
%
%
Jm_0=x_arg(22:39);
zeta_0=x_arg(40:57);
Jbr_0=Jm_0.*exp(j*zeta_0);
%
Vm_1=x_arg(59:67);
theta_1=x_arg(68:75);
Jm_1=x_arg(76:93);
zeta_1=x_arg(94:111);
Vbus_1=Vm_1.*exp(j*[0;theta_1]);
Jbr_1=Jm_1.*exp(j*zeta_1);
P_0=real(Vbus_0.*conj(AA*Jbr_0));
gentarget=P_0(genlist);

Q_0=imag(Vbus_0.*conj(AA*Jbr_0));
P_1=real(Vbus_1.*conj(AA*Jbr_1));
Q_1=imag(Vbus_1.*conj(AA*Jbr_1));
%

mismatch=[ P_0(genlist)-x_arg(1:3);
      P_1(genlist)-x_arg(1:3)+x_arg(58)*droop_gain;
      real(TJ_whole_0*Jbr_0+TE_whole_0*AA'*Vbus_0);
      imag(TJ_whole_0*Jbr_0+TE_whole_0*AA'*Vbus_0); 
      real(TJ_whole_1*Jbr_1+TE_whole_1*AA'*Vbus_1); ...
      imag(TJ_whole_1*Jbr_1+TE_whole_1*AA'*Vbus_1); ...
      P_0(loadlist)+real(snet(loadlist)); ...
      Q_0(loadlist)+imag(snet(loadlist));
      P_1(loadlist)+real(snet(loadlist)); ...
      Q_1(loadlist)+imag(snet(loadlist));
      x_arg(5)-x_arg(59);  %hold gen voltage setpoint from base case
      x_arg(8)-x_arg(62);
      x_arg(11)-x_arg(65)];


end

