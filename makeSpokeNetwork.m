function NT = makeSpokeNetwork(Rin,Rout,Nspk)
% spoke network with inner radius Rin, outer Rout, and Nspk number of
% spokes

dth=2*pi/Nspk; % spoke separation

% define outer circle
th=(0:dth:2*pi);
th=th(1:end-1);

% outer points
xout=Rout*cos(th)';
yout=Rout*sin(th)';

% inner points
xin=Rin*cos(th)';
yin=Rin*sin(th)';

NT=NetworkObj;
% set positions
NT.nodepos=[xin yin ; xout yout];

envec=reshape([round(1:.5:length(xin)) 1],2,[])';
envec2=[(1:length(xin))' (1:length(xin))'+length(xin)];
% set neighbors
NT.edgenodes=[envec ; envec+length(xin) ; envec2];
NT.setupNetwork(true);

end
