clear;
%命令字符串部分
ontStr='blastpgp -j 3 -e 1e-3 -i F:/RPI_SAN/PSSM/input.fa -d F:/RPI_SAN/PSSM/swissprot -Q F:/RPI_SAN/PSSM/output1807/';
twoStr='.txt';  
%threeStr='.out';

%装载数据
[header,seq]=fastaread('F:/RPI_SAN/RPI_SAN/ncRNA-protein/RPI1807_protein_seq.fa');   %换文件仅需更改这里
myMap = containers.Map;
count = 0;
for i=1429:length(header)    %序列的个数
    %n=1:30;
    %m=n;
    CmdStr='';
    seqName=char(header(i)); %以序列名命名输出文件
    %out=num2str(i);
    CmdStr=[ontStr,seqName,twoStr];   %拼接生成命令字符串
    FileStr=char(seq(i));   %序列
    %保存序列到文件
    fid=fopen('F:/RPI_SAN/PSSM/input.fa','w'); %输出路径
    fprintf(fid,'%s',FileStr);
    fclose(fid);
    system(CmdStr);  %执行命令
    
    outFileName = ['F:/RPI_SAN/PSSM/output1807/', seqName];
    loadMat=importdata([outFileName,'.txt']);
    seqPssmMat=loadMat.data(1:length(loadMat.data)-5,1:20);
    pfm = zeros(length(seqPssmMat),20);
    a= zeros(length(seqPssmMat));
    bigrams=zeros(20,20);
    for i =1:length(seqPssmMat)
        for j =1:20
            pfm(i,j) =double(1/(1+exp(-seqPssmMat(i,j))));  %归一化
            a(i) = a(i)+pfm(i,j);   
        end
    end
    for i = 1:length(seqPssmMat)
        for j = 1:20
            pfm(i,j) = double(pfm(i,j)*(1/a(i)));   %使每一行的行向量相加为1
        end
    end
    sum = 0;
    for m = 1:20
       for n = 1:20
           for i = 1:length(seqPssmMat)-1   
               sum = sum+pfm(i,m)*pfm(i+1,n);     
           end
           bigrams(m,n)=sum;     %求出bi-gram矩阵一个行向量
        sum = 0;
        end
    end

    %[A_nm,zmlist,cidx,V_nm]=P_zernike(seqPssmMat,n,m);
  feature = zeros(400,1);
  count = 1;
  for i=1:20
      for j=1:20
          feature = bigrams(i,j);
          count = count + 1;
      end
  end
    myMap(seqName) = feature;
    %feature{i} = abs(A_nm);
    %cell(seqName)={seqPssmMat};
end
save('F:/RPI_SAN/PSSM/1807.mat','mymap');
%save('F:/result/PSSM488.mat','cell');