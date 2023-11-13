function m=extracting_data(A,B_name,R,method,w,...%1-5
                            rate,...%6
                            S,T,G,NN,TN)%7-11
    name=["GLMR","Xing","IQIM","DCT_MME","LWT_MME"];
    m=struct(...
    "GLMR",...
    struct("isuse",0,"watermark",zeros(1,1)),...
    "Xing",...
    struct("isuse",0,"watermark",zeros(1,1)),...
    "IQIM",...
    struct("isuse",0,"watermark",zeros(1,1)),...
    "DCT_MME",...
    struct("isuse",0,"watermark",zeros(1,1)),...
    "LWT_MME",...
    struct("isuse",0,"watermark",zeros(1,1))...
    );
    
    if nargin<11
       TN=2;
    end
    if nargin<10
       NN=1;
    end
    if nargin<7
        S=300;
        T=100;
        G=7300;
    end
    if nargin<5
        rate=16;
    end
    
    [B,rp,rc,G1]=lattice_information(B_name,NN);
    N=size(B,1);
    [M1,N1]=size(A);
    len=floor(M1*N1/N);
    w1=reshape(w,1,[])*1;
    len=floor(size(w1,2));
%     w1=reshape(w1(1:len*N),N,[]);

    type=2^(rate-1);
    step=get_step;
    
    for i=1:length(method)
        j=method(i);
        m.(name{get_name(j)}).isuse=1;
        temp_A=A.(name{get_name(j)}).watermarked_signal;
        switch j
            case 1
                m.(name{get_name(j)}).watermark=extraction_process(temp_A,w,B_name,R,type,50);             
%                 m.(name{j}).watermark=MME_extraction_process(temp_A,w,B_name,R,type);
            case 2
                
                A1=reshape(temp_A,1,[])*type;
                wr=zeros(4,len);
                [wr(1,:),wr(2,:),wr(3,:),wr(4,:),Ew,dw]=Detecting_A(A1',S,T,G,len,TN);
                ber=zeros(1,4);
                for q=1:4
                    ber(q)=BER(w,wr(q,:));
                end
                [wber wNO]=min(ber);
                m.(name{get_name(j)}).watermark=wr(wNO,:);
            case 3
                m.(name{get_name(j)}).watermark=IQIM_extraction(temp_A,w,step,N,R,type);
            case 6
                m.(name{get_name(j)}).watermark=dct_extraction_process(temp_A,w,B_name,R,type,S,50);
            case 7
                m.(name{get_name(j)}).watermark=lwt_extraction_process(temp_A,w,B_name,R,type,50);
        end
        m.(name{get_name(j)}).watermark(w)=m.(name{get_name(j)}).watermark;
    end
end

function m=get_name(j)
    if j>=6
        m=j-2;
    else
        m=j;
    end
end