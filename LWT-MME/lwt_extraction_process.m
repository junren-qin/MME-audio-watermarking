function out=extraction_process(A,w,B_name,R,type,CoreNum)
    [B,rp,rc,G]=lattice_information(B_name);
    [M1,N1]=size(w);
    N=size(B,1);
    [ht,lt]=lwt(A,'db4');
    A1=reshape(ht,1,[])*type;
%     A1=reshape(A,1,[])*type;
    Alen=floor(size(A1,2)/N);
    A1=reshape(A1(1:Alen*N),N,[]);
    out=reshape(w,1,[])*1;
    len=floor(size(out,2)/N);
%     A1=dct(A1);
    out=reshape(out(1:len*N),N,[]);
    
%     ss=floor(300*length(A)/300000);
    ss=300;
    if nargin<6
%         A2=A1(:,floor(Alen/2):Alen);
        A2=A1;
%         for i=1:len
%             A2(:,i)=A1(:,(i-1)*ss+1);
%         end
        for i=1:len
%             p=A1(:,i);
            p=A2(:,i);
            t=sample_extraction(p,B,R);
            out(:,i)=t;

        end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
%             A2=A1(:,floor(Alen/2):Alen);
            A2=A1;
%             for i=1:len
%                 A2(:,i)=A1(:,(i-1)*ss+1);
%             end
            parfor i=1:len
                p=A2(:,i);
                t=sample_extraction(p,B,R);
                out(:,i)=t;
            end
%             delete(par);
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
