function m=IQIM_embedding(A,w,step,N,R,type)
    [M1,N1]=size(A);
    A1=reshape(A,1,[])*type;
    len=floor(size(A1,2)/N);
    A1=reshape(A1(1:len*N),N,[]);
    w1=reshape(w,1,[])*1;
    len=floor(size(w1,2)/N);
    w1=reshape(w1(1:len*N),N,[]);
    m=A1;
%     len=size(w1,2);
    
    ss=1;
    for i=1:len
%         p=A1(:,i);
        p=A1(:,(i-1)*ss+1);
        s=w1(:,i);
        t=e_IQIM(step,p,s,R);
        m(:,(i-1)*ss+1)=t;
%         m(:,i)=t;
    end
    m=reshape(m,M1,N1)/type;
end
function m=e_IQIM(step,p,s,R)
    N=size(p,1);
    k=floor(p./(R*step));
    a=p-R*step*k;
    m=R*step*k+step*s+(1/R)*a;   
end