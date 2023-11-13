function m=count_SWRdB(O,W)
    o=reshape(O,1,[]);
    w=reshape(W,1,[]);
    d=w-o;
    ps=mean(o.^2);
    pw=mean(d.^2);
    m=10*log10(ps/pw);
end