function out=embedding_process(A,w,B_name,a,R,type,CoreNum)
    [B,rp,rc,G]=lattice_information(B_name);
    [M1,N1]=size(A);
    N=size(B,1);
    A1=reshape(A,1,[])*type;
    len=floor(size(A1,2)/N);
    A1=reshape(A1(1:len*N),N,[]);
    w1=reshape(w,1,[])*1;
    len=floor(size(w1,2)/N);
    w1=reshape(w1(1:len*N),N,[]);
    out=A;
    
    ss=1;
    if nargin<7
        for i=1:len
%             p=A1(:,i);
            p=A1(:,(i-1)*ss+1);
            s=w1(:,i);
            t=sample_embedding(p,s,B,a,R);
%             A1(:,i)=t;
            A1(:,(i-1)*ss+1)=t;
        end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
%         for i=1:len
%             A2(:,i)=A1(:,(i-1)*300+1);
%         end
        A2=A1;
        parfor i=1:len
            p=A2(:,i);
            s=w1(:,i);
            t=sample_embedding(p,s,B,a,R);
            A2(:,i)=t;
        end
        A1=A2;
%         for i=1:len
%            A1(:,(i-1)*300+1)=A2(:,i);
%         end
%         delete(par);
    end
    
    A1=reshape(A1,1,[])/type;
    for i=1:length(A1)
        out(i)=A1(i);
    end
end
function m=sample_embedding(p,s,B,a,R)
    k=SDCVP(B*s-p,R*B);
    xi=B*s-R*B*k;
    m=p+a*(xi-p);
end
% function m=sample_embedding(p,s,B,a,R)
%     k=SDCVP(B*s-a*p,R*B);
%     xi=B*s-R*B*k;
%     m=xi+(1-a)*p;
% end