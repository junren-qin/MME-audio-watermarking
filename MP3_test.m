clc;clear all;
addpath('Lattice','MME','IQIM','Xing','Nishimura','DCT-MME','LWT-MME');
mname=["GLMR","Xing","IQIM","DCT_MME","LWT_MME"];
cover=load("cover_data.mat");
color={[37 151 213]/255,[4 149 174]/255,[3 107 168]/255,[22 60 126]/255,[0 191 255]/255};

name=cover.name;
cname=1:6;
tlen=1:6;
method=[1:3,6,7];
B_name="Z";
S=300;
T=7500;
G=20000;
rate=16;
R=2;
a_range=proposed_a(B_name,R)+.1;

ber=zeros(length(method),length(tlen));
NN=2;
TN=2;
bitrate_range=[64,96,128,192,256,320];
[B,rp,rc,Gi]=lattice_information(B_name,NN);
result=zeros(length(tlen),length(method),length(bitrate_range));
rSWR=result;
r_noise_mean=result;
r_noise_mean2=result;
r_ber=zeros(length(method),length(bitrate_range),4);
r_swr=r_ber;
for q=1:length(bitrate_range)
    bitrate=bitrate_range(q);
for a=1:length(a_range)
    aa=a_range(a);
    tic;
    for ii=1:length(tlen) 
        i=tlen(ii);
        name2="wav"+num2str(name(cname(i)));
        [watermarked_signal,w,perwl]=watermarked_signal_create(cover,name2,B_name,R,method,rate,aa,S,T,G,NN,TN);
        attacked_signal=mp3_compress(watermarked_signal,rate,method,mname,bitrate,char(name2));
        extracted_watermark=extracting_data(attacked_signal,B_name,R,method,perwl,rate,S,T,G,NN,TN);
        for tt=1:length(method)
            j=method(tt);
            result(ii,tt,q)=BER(extracted_watermark.(mname{get_name(j)}).watermark,w);
            A=cover.(name2).data';
            A=A(80001:380000);
            sita1(i)=NN*mean(reshape(A,1,[]).^2);
            rSWR(ii,tt,q)=count_SWRdB(A,watermarked_signal.(mname{get_name(j)}).watermarked_signal);
            [r_noise_mean(ii,tt,q),r_noise_mean2(ii,tt,q)]=noise_analysis(watermarked_signal.(mname{get_name(j)}).watermarked_signal,attacked_signal.(mname{get_name(j)}).watermarked_signal);
        end
    end
    time=toc;
end

    N=size(B,1);
    disp("time="+num2str(time)+"(s)"+"  bitrate: "+num2str(bitrate));
    temp_ber=result(:,:,q);
    temp_swr=rSWR(:,:,q);
    temp_noise1=r_noise_mean(:,:,q);
    temp_noise2=r_noise_mean2(:,:,q);
    k_ber(:,1)=transpose(min(temp_ber));
    k_ber(:,2)=transpose(mean(temp_ber));
    k_ber(:,3)=transpose(max(temp_ber));
    k_ber(:,4)=transpose(std(temp_ber));
    k_swr(:,1)=transpose(min(temp_swr));
    k_swr(:,2)=transpose(mean(temp_swr));
    k_swr(:,3)=transpose(max(temp_swr));
    k_swr(:,4)=transpose(std(temp_swr));
    noise1=transpose(mean(temp_noise1));
    noise2=transpose(mean(temp_noise2));
    disp("ber:"+num2str(k_ber(:,2)'));
    disp("swr:"+num2str(k_swr(:,2)'));
    disp("the mean of the noise caused by mp3 compression:"+num2str(noise1'));
    disp("the second order moment of the noise caused by mp3 compression:"+num2str(noise2'));
    for p=1:4
        r_ber(:,q,p)=k_ber(:,p);
        r_swr(:,q,p)=k_swr(:,p);
    end
end

kmove=0;%0.8;
legend_name=["MME","Liang","IQIM","DCT-MME","LWT-MME"];
figtype={'-d','-.o','--s','-.x','-^'};
figure
hold on
rtemp=r_ber(:,:,3)-r_ber(:,:,1);
for i=1:length(method)
    j=method(i);
    kk=[r_ber(i,:,1); rtemp(i,:)];
    va=plot(bitrate_range+(i-2)*kmove,r_ber(i,:,2),figtype{get_name(j)},'LineWidth',2);
    va.Color=color{get_name(j)};
    llegend(i,1)=va;
end
hold off
llegend1=reshape(llegend,1,[]);
legend(llegend1,legend_name{[1,2,3,4,5]});
xlabel('MP3 compression bit rate');
ylabel("BER");
set(gca,'Yscale','log');


x=(1:length(method));
figure
hold on
up=r_swr(:,1,3);
low=r_swr(:,1,1);
ktemp=[low,abs(up-low)]';
upb=barh(x,ktemp',.5,'stacked','FaceColor',color{1},'FaceAlpha',.5,'EdgeAlpha',.0);
upb(1,1).FaceAlpha=.0;
err=errorbar(r_swr(:,1,2),x,r_swr(:,1,4),'horizontal','.k','LineWidth',3);
h=scatter(r_swr(:,1,2),x,100,'sk','MarkerFaceColor','k','MarkerFaceAlpha',.9,'MarkerEdgeAlpha',.9,'LineWidth',3);
text(r_swr(:,1,2)-8,x+.16,num2str(r_swr(:,1,2))); 
hold off
set(gca,'YTick',1:5,'YTickLabel',["MME","Liang","IQIM","DCT-MME","LWT-MME"]);
legend([upb(1,2),err,h],["CI","SD","Average Value"],'NumColumns',3,'Location','southwest');
xlabel("SWR(dW)");

function [m,w,perwl]=watermarked_signal_create(cover,name,B_name,R,method,rate,a,S,T,G,NN,TN)
    A=cover.(name).data';
    A=A(80001:380000);
    fs=cover.(name).fs;
    [M1,N1]=size(A);
    rand('seed',0);
    w=randi([0 R-1],1,floor(M1*N1/S));
    wlen=1:size(w,2);
    perwl=randperm(size(w,2));
    w1=w(perwl);
    m=embedding_data(A,B_name,R,method,w1,rate,a,S,T,G,NN,TN);
end

function m=mp3_compress(embedding_result,rate,method,mname,bitrate,name)
    type=2^(rate-1);
    m=embedding_result;
    for i=1:length(method)
        j=method(i);
        signal=embedding_result.(mname{get_name(j)}).watermarked_signal;
        m.(mname{get_name(j)}).watermarked_signal=mp3compression(signal,44100,name,bitrate)';
    end
end

function m=get_name(j)
    if j>=6
        m=j-2;
    else
        m=j;
    end
end
function [m1,m2]=noise_analysis(O,W)
    temp=reshape(W-O,1,[]);
    m1=mean(temp);
    m2=10*log10(mean(reshape(O,1,[]).^2)/mean(temp.^2));
end