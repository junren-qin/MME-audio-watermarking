function m=embedding_data(A,B_name,R,method,w,...%1-5
                            rate,a,...%6-7
                            S,T,G,NN,TN,delta)%8-11
    name=["GLMR","Xing","IQIM","DE","Nishimura","DCT_MME","LWT_MME"];
    m=struct(...
    name{1},...
    struct("isuse",0,"watermarked_signal",zeros(1,1)),...
    name{2},...
    struct("isuse",0,"watermarked_signal",zeros(1,1)),...
    name{3},...
    struct("isuse",0,"watermarked_signal",zeros(1,1)),...
    name{4},...
    struct("isuse",0,"watermarked_signal",zeros(1,1)),...
    name{5},...
    struct("isuse",0,"watermarked_signal",zeros(1,1)),...
    name{6},...
    struct("isuse",0,"watermarked_signal",zeros(1,1)),...
    name{7},...
    struct("isuse",0,"watermarked_signal",zeros(1,1))...
    );
    if nargin<13
        delta=1;
    end
    if nargin<12
        TN=2;
    end
    if nargin<11
        NN=1;
    end
    if nargin<8
        S=300;
        T=100;
        G=7300;
    end
    if nargin<7
        a=Wang_a(B_name,R);
    end
    if nargin<6
        rate=16;
    end
    
    [B,rp,rc,G1]=lattice_information(B_name,NN);
    N=size(B,1);
    [M1,N1]=size(A);
    len=floor(M1*N1/N);
    type=2^(rate-1);
    step=get_step;
    
    for i=1:length(method)
        j=method(i);
        m.(name{j}).isuse=1;
        switch j
            case 1
                
                m.(name{j}).watermarked_signal=embedding_process(A,w,B_name,a,R,type,20);%%7 means set for corenum; 
%                 m.(name{j}).watermarked_signal=MME_embedding_process(A,w,B_name,a,R,type);%%7 means set for corenum; 
            case 2
                A1=reshape(A,1,[])*type;
                watermarked_signal=Embeding(A1',w,S,T,G,len,TN)';
                m.(name{j}).watermarked_signal=watermarked_signal/type;
            case 3
                m.(name{j}).watermarked_signal=IQIM_embedding(A,w,step,N,R,type);
            case 4%DE:2:1
                m.(name{j}).watermarked_signal=DE_embedding(A,w,type);
            case 5%Nishimura's method:n:1(order by the number of changeable host signal)
                m.(name{j}).watermarked_signal=Nishimura_embedding(A,w,step,type);
            case 6%DCT-MME
                m.(name{j}).watermarked_signal=dct_embedding_process(A,w,B_name,a,R,type,S,50);
            case 7%LWT-MME
                m.(name{j}).watermarked_signal=lwt_embedding_process(A,w,B_name,a,R,type,50);
        end
    end
end



