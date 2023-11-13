function [m,fs]=mp3compression(A,fs,name,bitrate)
    name_mp3=[name, '.mp3'];
    name_mp3_1=[name,'_stego','.mp3'];
    name_wav=[name, '.wav'];
    
%     A=[1,zeros(1,100),A,zeros(1,100),1];
    len=length(A);
%     A=[A(1:begin-1),zeros(1,20000),1,zeros(1,20000),A(begin:last),zeros(1,20000),1,zeros(1,20000),A(last+1:len)];
%     A=[1,zeros(1,20000-1),A,zeros(1,20000-1),1];
%     cd audio_mp3
%     audiowrite(name_wav,A,fs);
%     cd ..
%     file=audioload_1('audio_mp3',name_wav);
    file=audioload_mp3(name_wav,A,fs);
    audiosave(A,file,'.mp3',bitrate);
%     figure,plot(A);
%     mp3write(A,fs,name,[])

    file_mp3=audioload_1('audio_out',name_mp3_1);
    temp=file_mp3.data';
    nfs=file_mp3.fs;
    
    
%     cd audio_out
%     temp=audioread(name_mp3_1);
%     info=audioinfo(name_mp3_1);
%     cd ..
%     fs=info.SampleRate;
%     
    cd audio_mp3
    audiowrite(name_wav,temp,fs);
    m=audioread(name_wav);
    cd ..
% %     figure,plot(m);
%     len1=floor(len/44100*fs);
%     kk=find(m>.15);
%     m=m(kk(1)+floor(20000/44100*fs):end,1);
%     m=m(1:len1,1);
% %     figure,plot(m);
%     figure,plot(m);
end

function file=audioload_mp3(name,A,fs)
    PathName='E:\audio watermarking\audio_GLMR_robust\audio_mp3\';
    FileName=name;
    [file.path,file.name,file.ext] = fileparts([PathName,FileName]);
%     [file.data,file.fs,file.bitrate] = mp3read([PathName,FileName]);
    file.data=A';
    file.fs=fs;
    file.nbit=16;
%     file.bitrate=bitrate;
    [file.len, file.ch] = size(A');
end
