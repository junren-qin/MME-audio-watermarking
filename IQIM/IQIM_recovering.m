function m=IQIM_recovering(A,step,N,R)
    A1=reshape(A,1,[])*type;
    len=floor(size(A1,2)/N);
    A1=reshape(A1(1:len*N),N,[]);
    m=A1;
    len=size(A1,2);
    for i=1:len
        p=A1(:,i);
        t=r_IQIM(step,p,R);
        m(:,i)=t;
    end
end

function m=r_IQIM(step,p,R)
    k=floor(p./step);
    a=p-step*k;
    m=step*R*floor(k./R)+R*a;
end