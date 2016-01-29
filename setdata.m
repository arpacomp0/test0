nbus=9;
nelement=9;
nbranch=2*nelement;
genlist=[1 4 7];
loadlist=[2 3 5 6 8 9];
slacklist=[1]; % this is used only in a powerflow
               % that sets the initial condition;
               % the actual OPF does not use a slack bus
%
% Parameters for quadratic generator cost curves
%
% Cost will be expressed in terms of generator power.
% Common convention writes these cost functions with
% the argument in units of MWs. However, our constraint
% functions will all use generator active powers in per unit
% nomralized form, with a base of 1000 MVA.  Hence, the cost
% function evaluation will incluse a factor of 1000 multiplier
% on the per unit powers.
%
% $/hour operating cost = alpha + beta*(gen MW output) 
%                         + gamma*{(gen MW output)^2}

alph=[400;500;450];
beta=[25;30;28]; 
gamma=[.01;0.012;0.013];
%
% Parameters for generator min and max active power output
Pmin=[200;250;175];
Pmax=[600;650;900];
%
% Parameters for generator "droop"
% Actual generator output, "P_actual," in per unit, will be:
% P_actual=P_setpoint + (1/R)*(frequency error, per unit)
%
% The R values below are specified according to standard
% power systems engineering convention, as "percentage"
% values, but with a renomalization from the generator's
% power base, to the system's power base.  To understand this,
% first recognize that per unit frequency error may 
% be viewed as percentage frequency error relative to rated frequency.
% Hence, the R value of 0.04 used below, or 4%,  indicates that
% a 4 % frequency error induces a change in generator power output,
% away from setpoint, equal to 100% of the generator's rated power.
% Note that this requires a multiplier of 
% [Generator Rated Power]/[System Base Power].
% System base power here is 1000 MVA; generator rated power
% is assumed equal to Pmax above.
%
droop_gain=(1/0.04)*(1/1000)*Pmax;
%
AA=zeros(nbus,nbranch);
Zseries=zeros(nelement,1);
Yshunt=zeros(nelement,1);
Current_lim=zeros(nbranch,1);
TE_element=zeros(2,2,nelement);
TJ_element=zeros(2,2,nelement);
% element 1 - assumed to be a transformer - leakage reactance only
AA(1,1)=1;
AA(2,2)=1;
Zseries(1)=j*0.12;
Yshunt(1)=0;
% element 2 - assumed to be a transformer - leakage reactance only
AA(2,3)=1;
AA(3,4)=1;
Zseries(2)=j*0.11;
Yshunt(2)=0;
% element 3 - transmission line modeled with a series impedance
%             and a shunt admittance
AA(3,5)=1;
AA(5,6)=1;
Zseries(3)=0.5318+1i*2.486;
Yshunt(3)=j*6.128E-3;
% element 4 - assumed to be a transformer - leakage reactance only
AA(4,7)=1;
AA(5,8)=1;
Zseries(4)=j*0.3;
Yshunt(4)=0;
% element 5 - transmission line modeled with a series impedance
%             and a shunt admittance
AA(5,9)=1;
AA(6,10)=1;
Zseries(5)=0.5648+j*2.1073;  % old value 0.6648+j*3.1073
Yshunt(5)=j*7.66E-3;
% element 6 - assumed to be a transformer - leakage reactance only
AA(7,11)=1;
AA(6,12)=1;
Zseries(6)=1i*0.25;
Yshunt(6)=0;
% element 7 - transmission line modeled with a series impedance
%             and a shunt admittance
AA(6,13)=1;
AA(8,14)=1;
Zseries(7)=0.4432+j*2.0715;
Yshunt(7)=j*5.1E-3;
% element 8  - assumed to be a transformer - leakage reactance only
AA(8,15)=1;
AA(9,16)=1;
Zseries(8)=j*0.11;
Yshunt(8)=0;
% element 9 - transmission line modeled with a series impedance
%             and a shunt admittance
AA(2,17)=1;
AA(9,18)=1;
Zseries(9)=5.923E-2+j*0.6259;
Yshunt(9)=j*10.46E-2;
for kk=1:nelement
    A=1+0.5*Zseries(kk)*Yshunt(kk);
    B=Zseries(kk);
    C=Yshunt(kk)*(1+0.25*Zseries(kk)*Yshunt(kk));
    D=A;
    TE_element(1,1,kk)=1;
    TE_element(1,2,kk)=-A;
    TE_element(2,1,kk)=0;
    TE_element(2,2,kk)=-C;
    TJ_element(1,1,kk)=0;
    TJ_element(1,2,kk)=B;
    TJ_element(2,1,kk)=1;
    TJ_element(2,2,kk)=D;
    index1=2*kk-1;
    index2=2*kk;
    Current_lim(index1)=0.5*abs(inv(Zseries(kk)));
    Current_lim(index2)=Current_lim(index1);
%Manually set tighter bound on select line current
      pick_line=3;
      Current_lim(2*pick_line-1)=.1815;
      Current_lim(index2)=Current_lim(index1);
      pick_line=7;
      Current_lim(2*pick_line-1)=.1100;
      Current_lim(index2)=Current_lim(index1);
%
end
TE_whole_0=blkdiag(TE_element(:,:,1),TE_element(:,:,2), ...
                   TE_element(:,:,3),TE_element(:,:,4), ...
                   TE_element(:,:,5),TE_element(:,:,6), ...
                   TE_element(:,:,7),TE_element(:,:,8), ...
                   TE_element(:,:,9));
TJ_whole_0=blkdiag(TJ_element(:,:,1),TJ_element(:,:,2), ...
                   TJ_element(:,:,3),TJ_element(:,:,4), ...
                   TJ_element(:,:,5),TJ_element(:,:,6), ...
                   TJ_element(:,:,7),TJ_element(:,:,8), ...
                   TJ_element(:,:,9));

%
%  Contingency index 1 case created by removing element line_out
%  This forces a condition for which the two port currents
%  associated withe element line_out must be zero.  In the implicit
%  two-port representation, that implies TE=zero_matrix;
%  TJ=identity.
%
line_out=5;
% element line_out _ if lins below temporarily commented out, no contingency
    TE_element(1,1,line_out)=0;
    TE_element(1,2,line_out)=0;
    TE_element(2,1,line_out)=0;
    TE_element(2,2,line_out)=0;
    TJ_element(1,1,line_out)=1;
    TJ_element(1,2,line_out)=0;
    TJ_element(2,1,line_out)=0;
    TJ_element(2,2,line_out)=1;
%
%
%
TE_whole_1=blkdiag(TE_element(:,:,1),TE_element(:,:,2), ...
                   TE_element(:,:,3),TE_element(:,:,4), ...
                   TE_element(:,:,5),TE_element(:,:,6), ...
                   TE_element(:,:,7),TE_element(:,:,8), ...
                   TE_element(:,:,9));
TJ_whole_1=blkdiag(TJ_element(:,:,1),TJ_element(:,:,2), ...
                   TJ_element(:,:,3),TJ_element(:,:,4), ...
                   TJ_element(:,:,5),TJ_element(:,:,6), ...
                   TJ_element(:,:,7),TJ_element(:,:,8), ...
                   TJ_element(:,:,9));
%
%
TV=[TE_whole_0*AA' TE_whole_1*AA'];
TJ=[TJ_whole_0 TJ_whole_1];
% Matrices above will be used to form linearequality constraints
% associated with the two-port behavior of each element, over
% base case and contingency, in terms of bus volages and
% branch currents; that is:
%  [TV]*[v_solution_0;v_solution_1]
%      + [TJ]*[branch_currents_0;branch_currents_1] = 0
%
%
% This initial example does not include any element behavior
% that is incompatible with an admittance, "Ybus" description.
% As a debugging check against more traditional methods, the
% Ybus representation for the base case (index 0) and the
% contingency case (index 1) are constructed.  But recognize
% this step would NOT be possible if any of the elements
% above were of a type that does not allow an admittance
% representation; for example, a circuit breaker element.
%
yprim_0=-inv(full(TJ_whole_0))*full(TE_whole_0);
Ybus_0=sparse(AA*yprim_0*AA');
%
yprim_1=-inv(full(TJ_whole_1))*full(TE_whole_1);
Ybus_1=sparse(AA*yprim_1*AA');
%
% Observe that the incidence matrix, AA, that describes the
% network interconnection topology does NOT change from
% base case to any contingency case.  ONLY the element
% consituative equations, summarized in TE and TJ, are
% allowed to change; topology stays the same accross
% the base case and all contingency cases.
%
%
%
% Remaining data describes loads as a mix of constant
% complex impedance, constant complex current, and
% constant complex power.  These will be given in
% a complex matrix, "ZIP," of dimension busx3.
% In this initial text implementation, values will
% be used only at load buses (i.e., no load will be
% present at generator buses), and values will be specified
% that are the same across base case and all contingencies.
% Later improvements are anticipated to remove these restrictions.
% 
% A "demand" sign convention is assumed, wherein a postive
% quantity for current or power indicates a wthdrawal from
% the network, being consumed by the load element.
%
ZIP=zeros(nbus,3);
%
% Bus 2
ZIP(2,:)=[0+j*0, 0+j*0, 0.2+j*0.02];
ZIP(3,:)=[0+j*0, 0+j*0, 0.2+j*0.02];
ZIP(5,:)=[0+j*0, 0+j*0, 0.2+j*0.02];
ZIP(6,:)=[0+j*0, 0+j*0, 0.2+j*0.02];
ZIP(8,:)=[0+j*0, 0+j*0, 0.2+j*0.02];
ZIP(9,:)=[0+j*0, 0+j*0, 0.2+j*0.02];
%
%
% Solve "standard" power flow problem for base case and contingency:
%  Constant P-Q loads set by third column of "ZIP"
%  Generators' voltage manitudes set to 1.04 pu
%  Generator at bus 1 treated at slack;
%  Generator active power setpoints for gen #4 is midpoint of Pmin, Pmax;
%  Generator active power setpoints for gen #4 is midpoint of Pmin, Pmax;
%
vbus=ones(nbus,1);
vbus(1)=1.04;
vbus(4)=1.04;
vbus(7)=1.04;
snet=ZIP(:,3);
snet(4)=-((1/1000)*Pmin(2)+(0.5/1000)*(Pmax(2)-Pmin(2)));
snet(7)=-((1/1000)*Pmin(3)+(0.5/1000)*(Pmax(3)-Pmin(3)));
%
% Base case solution given by bus voltages "v_solved_0" below:
v_solved_0=pf_solve(Ybus_0,vbus,snet,slacklist,[4 7],loadlist);
%
% Find corresponding complex power delivered to network
% at each bus, consistent with base case solution above:
%
S_inject_0=v_solved_0.*conj(Ybus_0*v_solved_0);
%
branch_currents_0=yprim_0*AA'*v_solved_0;
%
% Contingency case solution (with line #5 out of service)
% given by bus voltages "v_solved_1" below:
v_solved_1=pf_solve(Ybus_1,vbus,snet,slacklist,[4 7],loadlist);
%
branch_currents_1=yprim_1*AA'*v_solved_1;
%
% The vector of decision variables is defined in partitions, as 
% described below
%
% First ngen components  - generator active power setpoints, p.u.
%
% 1 component - base case frequency error, p.u.
%
% Next nbus components - bus voltge magnitudes, p.u., base case network
% 
% Next (nbus-1) components - bus voltge angles, radians, base case network
%                            bus #1 is taken as angle reference

% 
% Next nbranch components - branch current magnitudes, base case network
% 
% Next nbranch components - branch current angles, base case network
%                           (note - bus #1 serves as angle reference
%                            for ALL angular quantities - hence no
%                            reduction in dimension here).
%
% The pattern of patitiooned variables below repeats for however
% many contingency configurations are being considered in the
% the constraint set.  For this example, only one contingency
% configuration is being considered
%
% 1 component - contingency case frequency error, p.u.
%
% Next nbus components - bus voltge magnitudes, p.u., contingency case 
% 
% Next (nbus-1) components - bus voltge angles, radians, contingency case
%                            bus #1 is taken as angle reference
% 
% Next nbranch components - branch current magnitudes, contingency case 
% 
% Next nbranch components - branch current angles, contingency case 
%
% So,  this example of 9 buses, 9 two-port elements (hence 18 currents
% when accounting for sending end and receiving end of the two-port),
% 3 generators, and the base case plus one contingency, yields a
% decision variable vector of form and dimension given by:
%
%   3 generator setpoints - these are same for base case and contingency; 
%     component indices 1:3
%   1 base case frequency error (note - held to bounds very near zero); 
%     component index 4
%   9 base case bus voltage magnitudes; 
%     component indices 5:13
%   8 base case bus voltage angles;
%     component indices 14:21
%  18 base case branch current magnitudes; 
%     component indices 22:39
%  18 base case branch current angles;
%     component indices 40:57
%   1 contingency case freq error; 
%     component index 58
%   9 contingency case bus voltage magnitudes; 
%     component indices 59:67
%   8 contingency case bus voltage angles;
%     component indices 68:75
%  18 contingency case branch current magnitudes;
%     component indices 76:93
%  18 contingency case branch current angles;
%     component indices 94:111
%
% 111 total variables
%
% COMMENT FOR FUTURE WORK: For initial development and ease of debugging,
% index locations above are "hard-wired" for this examle.  A key step
% in generalizing for arbtrary dimensioned networks will be the creation
% of a flexible indexing scheme for the decision variable vector.
%
% Set initial guess values for decision variables
%
xx0=zeros(111,1);
xx0(1)=real(S_inject_0(1)); % Gen @ bus 1 power setpoint
xx0(2)=real(S_inject_0(4)); % Gen @ bus 4 power setpoint
xx0(3)=real(S_inject_0(7)); % Gen @ bus 7 power setpoint
xx0(4)=0;  % freq error base case
xx0(5:13)=abs(v_solved_0);
xx0(14:21)=angle(v_solved_0(2:9));
xx0(22:39)=abs(branch_currents_0);
xx0(40:57)=angle(branch_currents_0);
xx0(58)=0;
xx0(59:67)=abs(v_solved_1);
xx0(68:75)=angle(v_solved_1(2:9));
xx0(76:93)=abs(branch_currents_1);
xx0(94:111)=angle(branch_currents_1);
%
LB=-1000*ones(111,1);
UB=1000*ones(111,1);
LB(1:3)=(1/1000)*Pmin;
UB(1:3)=(1/1000)*Pmax;
LB(4)=-0.0001; % 0.0001 pu or 0.01% frequency error allowed base case
UB(4)=0.0001;
LB(5:13)=0.94*ones(9,1);
UB(5:13)=1.06*ones(9,1);
LB(58)=-0.02; % 0.02 pu or 5% frequency error allowed in contingency
UB(58)=0.02;
LB(59:67)=0.92*ones(9,1); % Voltage magnitude min of 0.92 pu
UB(59:67)=1.08*ones(9,1); % Voltage magnitude max of 1.08 pu
LB(22:39)=0;
UB(22:39)=Current_lim;
LB(76:93)=0;
UB(76:93)=Current_lim;













