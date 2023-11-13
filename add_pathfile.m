clc;clear all;
cd test_audio;
main_num=1:2;
num=cell(1,2);
num{1}=1:5;
num{2}=1:70;
k2=size(num{1},2)+size(num{2},2);
original=cell(3,k2);
for i=1:size(main_num,2)
    k=size(num{i},2);
    for j=1:k
        name1=i*1000+j;
        name2="test_audio"+num2str(i)+" ("+num2str(j)+").wav";
        k1=(i-1)*size(num{1},2)+j;
        original{1,k1}=name1;
        [original{3,k1},original{2,k1}]=audioread(name2);
    end
end
cd ..
name=original{1,1};
for i=2:k2
    name3=original{1,i};
    name=[name,name3];
end
wav="wav";
t=cell(1,k2);
for i=1:k2
    t{i}=struct(wav+num2str(name(i)),{original{2,i},original{3,i}});
end
save("cover_data.mat","name");
for i=1:k2
    sname=wav+num2str(name(i));
    sdata=struct("fs",original{2,i},"data",original{3,i});
    eval(sname+"=sdata;");
    save("cover_data.mat",sname,"-append");
end

system("mkdir audio_mp3")
