function outcome = funcminf(Vm, Vmx, Smx)
outcome = 1./(1+exp(-(Vm-Vmx)./Smx));