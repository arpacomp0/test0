vbus0=x(5:13).*exp(j*[0;x(14:21)]);
Jbranch0=x(22:39).*exp(j*x(40:57));
S_injections_0=vbus0.*conj(AA*Jbranch0);
vbus1=x(59:67).*exp(j*[0;x(68:75)]);
Jbranch1=x(76:93).*exp(j*x(94:111));
S_injections_1=vbus1.*conj(AA*Jbranch1);
disp('The base case bus voltage solution for the study network is:')
disp('    Bus#      V p.u.  Angle-degrees')
b_indices=(1:nbus)';
disp([b_indices abs(vbus0)  (180/pi)*angle(vbus0)])
%disp('With associated per unit active and reactive power injections of:')
%disp('(Note that the system power base is 1000 MVA, so p.u. quantities')
%disp('below multiplied by 1000 yield actual MW or MVar)')


%disp('    Bus#      P pu      Q pu')
%b_indices=(1:nbus)';
%disp([b_indices real(S_injections_0)  imag(S_injections_0)])
%
op_cost = Gencost(x,alph,beta,gamma);
disp('The operating cost (evaluated for base case disptach) is:'), ...
    op_cost, disp('in units of $/hour') 
%
base_losses=1000*sum(real(vbus0.*conj(AA*Jbranch0)));
disp('The system losses for the base case are:'), ...
base_losses, disp('in MWs')
%
disp('The contingency case bus voltage solution is:')
disp('    Bus#      V p.u.  Angle-degrees')
b_indices=(1:nbus)';
disp([b_indices abs(vbus1)  (180/pi)*angle(vbus1)])
%disp('With associated per unit active and reactive power injections of:')
%disp('(Note that the system power base is 1000 MVA, so p.u. quantities')
%disp('below multiplied by 1000 yield actual MW or MVar)')
%
%
%disp('    Bus#      P pu      Q pu')
%b_indices=(1:nbus)';
%disp([b_indices real(S_injections_1)  imag(S_injections_1)])
%
contingency_case_losses=1000*sum(real(vbus1.*conj(AA*Jbranch1)));
disp('The system losses for the contingency case are:'), ...
contingency_case_losses, disp('in MWs')
freq_error_contingency=60*x(58);
generator_droop_response=1000*x(58)*droop_gain;
disp('The system frequency error for the contingency case is:'), ...
freq_error_contingency, disp('in Hz')
%
disp_indices=[1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9]';
disp('The magnitudes of line flow currents are shown below')
disp('for base case and contingency, followed by the limit value.')
disp('Each line has two rows: first is sending-end current magnitude,')
disp('                        second is receiving-end current magnitude')
disp('    line #  base-case  contingency limit')
disp([ disp_indices x(22:39) x(76:93) Current_lim])
