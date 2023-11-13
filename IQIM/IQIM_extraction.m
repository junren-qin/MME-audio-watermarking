function m=IQIM_extraction(A,w,step,N,R,type)
    A1=reshape(A,1,[])*type;
    len=floor(size(A1,2)/N);
    A1=reshape(A1(1:len*N),N,[]);
%     m=A1;
    len=size(A1,2);
    out=reshape(w,1,[])*1;
    len=floor(size(out,2)/N);
    m=reshape(out(1:len*N),N,[]);

    ss=1;
    for i=1:len
%         p=A1(:,i);
        p=A1(:,(i-1)*ss+1);
        t=ex_IQIM(step,p,R);
        m(:,i)=t;
    end
    m=reshape(m,1,[]);
end

function w=ex_IQIM(step,p,R)
    k=floor(p./step);
    w=mod(k,R);
end