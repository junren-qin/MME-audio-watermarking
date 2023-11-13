function out=Nishimura_embedding(A,w,alpha,type)
    [M1,N1]=size(A);
    A1=reshape(A,1,[])*type;
    w1=reshape(w,1,[])*1;
    k=1;
    for i=1:length(w1)
        h=A1(:,i);
        k1=changeable(alpha,h);
        if k1==1
            s=w1(:,k);
            t=Ni_reQIM(alpha,h,s);
            A1(:,i)=t;
            k=k+1;
            if k>length(w1)
                break;
            end
        end
    end
    out=reshape(A1/type,M1,N1);
end

%--------------------------------------------------------------------------
%Nishimura's method
function m=E(a,h)
    m=round(a*h);
end
function m=changeable(alpha,h)
    if (E(alpha,h-1)~=E(alpha,h)-1&h>0)|(E(alpha,h+1)~=E(alpha,h)+1&h<0)
        m=1;
    else
        m=0;
    end
end
function t=Ni_reQIM(alpha,h,s)%Nishimura's method:embedding for 1 sample
    if E(alpha,h-1)~=E(alpha,h)-1&h>0
        t=E(alpha,h)-s;
    elseif E(alpha,h+1)~=E(alpha,h)+1&h<0
        t=E(alpha,h)+s;
    else
        t=E(alpha,h);
    end
end
