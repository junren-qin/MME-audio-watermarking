function out=MME_extraction_process(A,w,B_name,R,type,CoreNum)
    [B,rp,rc,G]=lattice_information(B_name);
    [M1,N1]=size(w);
    N=size(B,1);
    A1=reshape(A,1,[])*type;
    len=floor(size(A1,2)/N);
    A1=reshape(A1(1:len*N),N,[]);
    out=reshape(w,1,[])*1;
    len=floor(size(out,2)/N);
    out=reshape(out(1:len*N),N,[]);
    
%     ss=floor(300*length(A)/300000);
    ss=1;
    S=300;
    if nargin<7
        for i=1:len
            temp=A1(:,(i-1)*S+1:i*S);
            kmean=mean(transpose(temp))';
            kerror=abs(temp-kmean);
            p=sum(transpose(kerror))';
            out(:,i)=sample_extraction(p,B,R);
        end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
            for i=1:len
                temp=A1(:,(i-1)*S+1:i*S);
                kmean=mean(transpose(temp))';
                kerror=abs(temp-kmean);
                A2(:,i)=sum(transpose(kerror))';
            end
            parfor i=1:len
                p=A2(:,i);
                t=sample_extraction(p,B,R);
                out(:,i)=t;
            end
            delete(par);
    end   
    out=reshape(out,M1,N1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
end
% function m=sample_extraction(p,B,R,a)
%     k=SDCVP(a*p,B);
%     m=mod(k,R);
% end
function m=sample_extraction(p,B,R)
    k=SDCVP(p,B);
    m=mod(k,R);
end