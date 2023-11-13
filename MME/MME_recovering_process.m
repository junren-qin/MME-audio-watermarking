function out=MME_recovering_process(A,w,B_name,a,type,CoreNum)
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
    
    if nargin<6
        for i=1:len
            temp=A1(:,(i-1)*S+1:i*S);
            kmean=mean(transpose(temp))';
            kerror=temp-kmean;
            ksign=sign(kerror);
            p=sum(transpose(abs(kerror)))';
            rtt=sample_recovering(p,B,alpha);
            rt=rtt-p;
            inverse_square_signal=ksign.*(rt/S);
            A1(:,(i-1)*S+1:i*S)=(kerror+inverse_square_signal+kmean);
        end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
            for i=1:len
                temp=A1(:,(i-1)*S+1:i*S);
                kmean=mean(transpose(temp))';
                kerror=temp-kmean;
                ksign(:,(i-1)*S+1:i*S)=sign(kerror);
                A2(:,i)=sum(transpose(abs(kerror)))';
            end
            parfor i=1:len
                p=A2(:,i);
                t=sample_recovering(p,B,a);
                A2(:,i)=t;
            end
            for i=1:len
                wtemp=A2(:,i);
                inverse_square_signal=ksign(:,(i-1)*S+1:i*S).*(wtemp/S);
                A1(:,(i-1)*S+1:i*S)=A1(:,(i-1)*S+1:i*S)+inverse_square_signal;
            end
    end
    
    A1=reshape(A1,1,[])/type;
    for i=1:length(A1)
        out(i)=A1(i);
    end
end
function m=sample_recovering(p,B,a)
    k=SDCVP(p,B);
    xi=B*k;
    m=(1/(1-a))*p-(a/(1-a))*xi;
end
% function m=sample_recovering(p,B,a)
%     k=SDCVP(a*p,B);
%     xi=B*k;
%     m=(1/(1-a))*p-(1/(1-a))*xi;
% end