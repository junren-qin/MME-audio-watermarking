%%%%%%%%%%%%%%%%%%%������������%%%%%%%%%%%%%%%%%%%%
function [dBer] = BER( wm,wr )
[M,N]=size(wm);
nDiff= find(wm~=wr);
dBer = length(nDiff)/(M*N);
end