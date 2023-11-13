function out=MME_embedding_process(A,w,B_name,a,R,type,CoreNum)
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
    S=300;
    if nargin<7
        for i=1:len
            temp=A1(:,(i-1)*S+1:i*S);
            kmean=mean(transpose(temp))';
            kerror=temp-kmean;
            ksign=sign(kerror);
            p=sum(transpose(abs(kerror)))';
            s=w1(:,i);
            t=sample_embedding(p,s,B,a,R);
            wtemp=t-p;
            square_signal=ksign.*(wtemp/S);
            A1(:,(i-1)*S+1:i*S)=temp+square_signal;
        end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
%         for i=1:len
%             A2(:,i)=A1(:,(i-1)*300+1);
%         end
        for i=1:len
            temp=A1(:,(i-1)*S+1:i*S);
            kmean=mean(transpose(temp))';
            kerror=temp-kmean;
            ksign(:,(i-1)*S+1:i*S)=sign(kerror);
            A2(:,i)=sum(transpose(abs(kerror)))';
        end
        parfor i=1:len
            p=A2(:,i);
            s=w1(:,i);
            t=sample_embedding(p,s,B,a,R);
            A2(:,i)=t-p;
        end
        for i=1:len
            wtemp=A2(:,i);
            square_signal=ksign(:,(i-1)*S+1:i*S).*(wtemp/S);
            A1(:,(i-1)*S+1:i*S)=A1(:,(i-1)*S+1:i*S)+square_signal;
        end
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