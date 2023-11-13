function out=dct_embedding_process(A,w,B_name,a,R,type,CoreNum)
    [B,rp,rc,G]=lattice_information(B_name);
    [ht,lt]=lwt(A,'db4');
    [M1,N1]=size(ht);
    A1=reshape(ht,1,[])*type;
%     [M1,N1]=size(A);
    N=size(B,1);
%     A1=reshape(A,1,[])*type;
    Alen=floor(size(A1,2)/N);
    A1=reshape(A1(1:Alen*N),N,[]);
    w1=reshape(w,1,[])*1;
    len=floor(size(w1,2)/N);
    w1=reshape(w1(1:len*N),N,[]);
%     A1=dct(A1);
    out=ht;
    ss=300;
    if nargin<7
        for i=1:len
            p=A1(:,i);
%             p=A1(:,(i-1)*ss+1);
%             p=A1(:,floor(Alen/2)+i-1);
            s=w1(:,i);
            t=sample_embedding(p,s,B,a,R);
            
            A1(:,i)=t;
%             A1(:,(i-1)*ss+1)=t;
%             A1(:,floor(Alen/2)+i-1)=t;
        end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
%         for i=1:len
%             A2(:,i)=A1(:,(i-1)*ss+1);
%         end
%         A2=A1(:,floor(Alen/2):Alen);
        A2=A1;
        parfor i=1:len
            p=A2(:,i);
            s=w1(:,i);
            t=sample_embedding(p,s,B,a,R);
            A2(:,i)=t;
        end
        A1=A2;
%         A1(:,floor(Alen/2):Alen)=A2;
%         for i=1:len
%            A1(:,(i-1)*ss+1)=A2(:,i);
%         end
%         delete(par);
    end
%     A1=idct(A1);
    A1=reshape(A1,1,[])/type;
    for i=1:length(A1)
        out(i)=A1(i);
    end
    out=ilwt(out,lt,'db4');
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
