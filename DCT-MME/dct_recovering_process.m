function out=recovering_process(A,w,B_name,a,type,win,CoreNum)
    [B,rp,rc,G]=lattice_information(B_name);
    [M1,N1]=size(A);
    N=size(B,1);
    A1=reshape(A,1,[])*type;
    Alen=floor(size(A1,2)/N);
    A1=reshape(A1(1:Alen*N),N,[]);
    w1=reshape(w,1,[])*1;
    len=floor(size(w1,2)/N);
    w1=reshape(w1(1:len*N),N,[]);
%     A1=dct(A1);
    out=A;
    ss=300;
    if nargin<7
        for i=1:len
            [ptemp,pb,pl]=make_window(A1,win,i);
            ptemp=dct(ptemp);
            A2(:,i)=ptemp(:,ceil(size(ptemp,2)/2));
        end
%         for i=1:len
%             A2(:,i)=A1(:,(i-1)*ss+1);
%         end
        for i=1:len
            p=A2(:,i);
%             p=A1(:,floor(Alen/2)+i-1);
            t=sample_recovering(p,B,a);
            A2(:,i)=t;
%             A1(:,floor(Alen/2)+i-1)=t;
        end
        for i=1:len
            [ptemp,pb,pl]=make_window(A1,win,i);
            ptemp=dct(ptemp);
            ptemp(:,ceil(size(ptemp,2)/2))=A2(:,i);
            ptemp=idct(ptemp);
            A1(:,pb:pl)=ptemp;
        end

%         for i=1:len
%             A1(:,(i-1)*ss+1)=A2(:,i);
%         end
    else
            if isempty(gcp('nocreate'))
                par=parpool(CoreNum);
            end
            for i=1:len
                [ptemp,pb,pl]=make_window(A1,win,i);
                ptemp=dct(ptemp);
                A2(:,i)=ptemp(:,ceil(size(ptemp,2)/2));
            end
%             for i=1:len
%                 A2(:,i)=A1(:,(i-1)*ss+1);
%             end
%             A2=A1;
%             A2=A1(:,floor(Alen/2):Alen);
            parfor i=1:len
                p=A2(:,i);
                t=sample_recovering(p,B,a);
                A2(:,i)=t;
            end
            for i=1:len
                [ptemp,pb,pl]=make_window(A1,win,i);
                ptemp=dct(ptemp);
                ptemp(:,ceil(size(ptemp,2)/2))=A2(:,i);
                ptemp=idct(ptemp);
                A1(:,pb:pl)=ptemp;
            end

%             A1=A2;
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
function [m,begin,last]=make_window(data,window_size,i)
    len=size(data,2);
    if i*window_size<=len
        m=data(:,(i-1)*window_size+1:i*window_size);
        begin=(i-1)*window_size+1;
        last=i*window_size;
    else
        m=data(:,(i-1)*window_size+1:len);
        begin=(i-1)*window_size+1;
        last=len;
    end
end