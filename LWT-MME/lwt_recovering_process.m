function out=recovering_process(A,w,B_name,a,type,CoreNum)
    [B,rp,rc,G]=lattice_information(B_name);
    [ht,lt]=lwt(A,'db2');
    [M1,N1]=size(ht);
%     [M1,N1]=size(A);
    N=size(B,1);
    A1=reshape(lt,1,[])*type;
%     A1=reshape(A,1,[])*type;
    Alen=floor(size(A1,2)/N);
    A1=reshape(A1(1:Alen*N),N,[]);
    w1=reshape(w,1,[])*1;
    len=floor(size(w1,2)/N);
    w1=reshape(w1(1:len*N),N,[]);
%     A1=dct(A1);
%     [ht,lt]=lwt(A1,'db2');
%     out=A;
    out=ht;
    ss=300;
    if nargin<6
%         for i=1:len
%             A2(:,i)=A1(:,(i-1)*ss+1);
%         end
        A2=A1;
        for i=1:len
            p=A2(:,i);
%             p=A1(:,floor(Alen/2)+i-1);
            t=sample_recovering(p,B,a);
            A2(:,i)=t;
%             A1(:,floor(Alen/2)+i-1)=t;
        end
        A1=A2;
%         for i=1:len
%             A1(:,(i-1)*ss+1)=A2(:,i);
%         end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
%             for i=1:len
%                 A2(:,i)=A1(:,(i-1)*ss+1);
%             end
            A2=A1;
%             A2=A1(:,floor(Alen/2):Alen);
            parfor i=1:len
                p=A2(:,i);
                t=sample_recovering(p,B,a);
                A2(:,i)=t;
            end
            A1=A2;
%             A1(:,floor(Alen/2):Alen)=A2;
%             for i=1:len
%                 A1(:,(i-1)*ss+1)=A2(:,i);
%             end
%             delete(par);
    end
%     A1=idct(A1);
    A1=reshape(A1,1,[])/type;
    for i=1:length(A1)
        out(i)=A1(i);
    end
    out=ilwt(out,lt,'db2');
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
