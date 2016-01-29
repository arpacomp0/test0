function [ deriv ] = dffeq_dx(x_arg, ...
                   AA,TE_whole_0,TJ_whole_0, ...
                   TE_whole_1,TJ_whole_1,genlist, ...
                   loadlist,snet,droop_gain)
%
%
%
gen_set=x_arg(1:3); % Generator active power setpoints are variables 1:3
                    % System frequency in base case is variable 4(not used)
Vm_0=x_arg(5:13);   % Base case bus voltage magnitudes are variables 5:13
theta_0=x_arg(14:21); % Base case bus voltage phase angles are variables 14:21
                      % (note that the angle at bus 1 is the reference angle for 
                      % for all phasor quantities in the bases case solution.
                      % Hence it is identically zero, and is NOT among the
                      % variables; i.e., only angles at buses 2 through 9 are
                      % variables.
Vbus_0=Vm_0.*exp(j*[0;theta_0]); % construct complex bus voltage phasor
                                 % as an intermediate variable (base case).
%
%
Jm_0=x_arg(22:39);  % Base case port current magnitudes are variables 22:39
zeta_0=x_arg(40:57); % Base case port current phase angles are variables 40:57
                      % Recall that the VOLTAGE angle at bus 1 is the reference angle for 
                      % for all phasor quantities in the bases case solution.
                      % Hence ALL 18 port current phase angles ARE variables,
                      % with their angles defined relative to the reference
                      % of bus voltage 1's phase.
Jbr_0=Jm_0.*exp(j*zeta_0);  % construct complex port current phasor
                            % as an intermediate variable (base case).
%
                      % System frequency in the contingency case is
                      % variable 58
%
Vm_1=x_arg(59:67); % Contingency case bus voltage magnitudes are variables 59:67
theta_1=x_arg(68:75); % Contingency case bus voltage phaseangles are variables 68:75.
                      % Observe that the contingency case is a new
                      % operating point, distinct from the base case, and
                      % hence has anew reference angle.  This new reference
                      % is the CONTINGENCY case bus 1 phase angle.  So as
                      % above, bus angle 1 is identically zero in the
                      % contingency case, and is not in the variable set,
                      % leaving only 8 phase angle unknowns.
Jm_1=x_arg(76:93);    % Contingency case port current magnitudes are variables 76:93
zeta_1=x_arg(94:111); % Contingency case port current phase angles are variables 94:11
%
% Below is construction of additional complex quantities that
% serve as convenient intermediate variables in subsequent
% derivative computations
%
Vbus_1=Vm_1.*exp(j*[0;theta_1]);
Jbr_1=Jm_1.*exp(j*zeta_1);
P_0=real(Vbus_0.*conj(AA*Jbr_0));
Q_0=imag(Vbus_0.*conj(AA*Jbr_0));
P_1=real(Vbus_1.*conj(AA*Jbr_1));
Q_1=imag(Vbus_1.*conj(AA*Jbr_1));
%
%
%  Below are computations of sub-matrices that will appear as
%  blocks in the overall Jacobian matrix of partial derivatives
%  being computed
dP0_dVm0=real(diag(conj(AA*Jbr_0))*diag(exp(j*[0;theta_0])));
dP0_dtheta0=real(j*diag(conj(AA*Jbr_0))*[zeros(1,8);diag(Vbus_0(2:9))]);
dP0_dJm0=real(diag(Vbus_0)*AA*diag(exp(-j*zeta_0)));
dP0_dzeta0=real([-j*diag(Vbus_0)*AA*diag(conj(Jbr_0))]);
%
dQ0_dVm0=imag(diag(conj(AA*Jbr_0))*diag(exp(j*[0;theta_0])));
dQ0_dtheta0=imag(j*diag(conj(AA*Jbr_0))*[zeros(1,8);diag(Vbus_0(2:9))]);
dQ0_dJm0=imag(diag(Vbus_0)*AA*diag(exp(-j*zeta_0)));
dQ0_dzeta0=imag([-j*diag(Vbus_0)*AA*diag(conj(Jbr_0))]);
%
%
dP1_dVm1=real(diag(conj(AA*Jbr_1))*diag(exp(j*[0;theta_1])));
dP1_dtheta1=real(j*diag(conj(AA*Jbr_1))*[zeros(1,8);diag(Vbus_1(2:9))]);
dP1_dJm1=real(diag(Vbus_1)*AA*diag(exp(-j*zeta_1)));
dP1_dzeta1=real([-j*diag(Vbus_1)*AA*diag(conj(Jbr_1))]);
%
%
dQ1_dVm1=imag(diag(conj(AA*Jbr_1))*diag(exp(j*[0;theta_1])));
dQ1_dtheta1=imag(j*diag(conj(AA*Jbr_1))*[zeros(1,8);diag(Vbus_1(2:9))]);
dQ1_dJm1=imag(diag(Vbus_1)*AA*diag(exp(-j*zeta_1)));
dQ1_dzeta1=imag([-j*diag(Vbus_1)*AA*diag(conj(Jbr_1))]);
%
dportr_dVm0=real(TE_whole_0*AA'*diag(exp(j*[0;theta_0])));
dportr_dtheta0=real(j*TE_whole_0*AA'*[zeros(1,8);diag(Vbus_0(2:9))]);
dportr_dJm0=real(TJ_whole_0*diag(exp(j*zeta_0)));
dportr_dzeta0=real(j*TJ_whole_0*diag(Jbr_0));
%
%
dporti_dVm0=imag(TE_whole_0*AA'*diag(exp(j*[0;theta_0])));
dporti_dtheta0=imag(j*TE_whole_0*AA'*[zeros(1,8);diag(Vbus_0(2:9))]);
dporti_dJm0=imag(TJ_whole_0*diag(exp(j*zeta_0)));
dporti_dzeta0=imag(j*TJ_whole_0*diag(Jbr_0));
%
dportr_dVm1=real(TE_whole_1*AA'*diag(exp(j*[0;theta_1])));
dportr_dtheta1=real(j*TE_whole_1*AA'*[zeros(1,8);diag(Vbus_1(2:9))]);
dportr_dJm1=real(TJ_whole_1*diag(exp(j*zeta_1)));
dportr_dzeta1=real(j*TJ_whole_1*diag(Jbr_1));
%
dporti_dVm1=imag(TE_whole_1*AA'*diag(exp(j*[0;theta_1])));
dporti_dtheta1=imag(j*TE_whole_1*AA'*[zeros(1,8);diag(Vbus_1(2:9))]);
dporti_dJm1=imag(TJ_whole_1*diag(exp(j*zeta_1)));
dporti_dzeta1=imag(j*TJ_whole_1*diag(Jbr_1));
%
deriv=[-eye(3) zeros(3,1) dP0_dVm0(genlist,:) dP0_dtheta0(genlist,:) dP0_dJm0(genlist,:) dP0_dzeta0(genlist,:) zeros(3,54); % derivatives of base case generator active power constraints
       -eye(3) zeros(3,54) droop_gain dP1_dVm1(genlist,:) dP1_dtheta1(genlist,:)  dP1_dJm1(genlist,:) dP1_dzeta1(genlist,:); % derivatives of contingency case generator active power constraints
       zeros(18,4) dportr_dVm0 dportr_dtheta0 dportr_dJm0 dportr_dzeta0 zeros(18,54); % derivatives of base case two-port element constraints, real part
       zeros(18,4) dporti_dVm0 dporti_dtheta0 dporti_dJm0 dporti_dzeta0 zeros(18,54); % derivatives of base case two-port element constraints, imaginary part
       zeros(18,58) dportr_dVm1 dportr_dtheta1 dportr_dJm1 dportr_dzeta1; % derivatives of contingency case two-port element constraints, real part
       zeros(18,58) dporti_dVm1 dporti_dtheta1 dporti_dJm1 dporti_dzeta1; % derivatives of contingency case two-port element constraints, imaginary part
       zeros(6,4) dP0_dVm0(loadlist,:) dP0_dtheta0(loadlist,:) dP0_dJm0(loadlist,:) dP0_dzeta0(loadlist,:) zeros(6,54); % derivatives of load buses active power constraints, base case
       zeros(6,4) dQ0_dVm0(loadlist,:) dQ0_dtheta0(loadlist,:) dQ0_dJm0(loadlist,:) dQ0_dzeta0(loadlist,:) zeros(6,54); % derivatives of load buses reactive power constraints, base case
       zeros(6,58) dP1_dVm1(loadlist,:) dP1_dtheta1(loadlist,:)  dP1_dJm1(loadlist,:) dP1_dzeta1(loadlist,:); % derivatives of load buses active power constraints, contingency case
       zeros(6,58) dQ1_dVm1(loadlist,:) dQ1_dtheta1(loadlist,:)  dQ1_dJm1(loadlist,:) dQ1_dzeta1(loadlist,:); % derivatives of load buses reactive power constraints, contingency case
       zeros(1,4) 1 zeros(1,53) -1 zeros(1,52); % derivative of constraint holding generator voltage magnitude setpoint in continegncy case equal to base case
       zeros(1,7) 1 zeros(1,53) -1 zeros(1,49); % derivative of constraint holding generator voltage magnitude setpoint in continegncy case equal to base case
       zeros(1,10) 1 zeros(1,53) -1 zeros(1,46)];% % derivative of constraint holding generator voltage magnitude setpoint in continegncy case equal to base case


end