function contSig = controlSignal(t,amp,time,step,npulse,offset)
contSig =  amp*pulstran(t,(offset:time(end)/npulse:time(end)),...
    rectpuls(time,time(end)*step/npulse));