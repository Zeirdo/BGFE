clear;
%�����ַ�������
ontStr='blastpgp -j 3 -e 1e-3 -i F:/RPI_SAN/PSSM/input.fa -d F:/RPI_SAN/PSSM/swissprot -Q F:/RPI_SAN/PSSM/output369/';
twoStr='.txt';  
%threeStr='.out';

%װ������
[header,seq]=fastaread('F:/ipminer/ncRNA-protein/RPI369_protein.fa');   %���ļ������������
myMap = containers.Map;
count = 0;
for i=242:length(header)    %���еĸ���
    %n=1:30;
    %m=n;
    CmdStr='';
    seqName=char(header(i)); %����������������ļ�
    %out=num2str(i);
    CmdStr=[ontStr,seqName,twoStr];   %ƴ�����������ַ���
    FileStr=char(seq(i));   %����
    %�������е��ļ�
    fid=fopen('F:/RPI_SAN/PSSM/input.fa','w'); %���·��
    fprintf(fid,'%s',FileStr);
    fclose(fid);
    system(CmdStr);  %ִ������
    
    outFileName = ['F:/RPI_SAN/PSSM/output369/', seqName];
    loadMat=importdata([outFileName,'.txt']);
    seqPssmMat=loadMat.data(1:length(loadMat.data)-5,1:20);
%     pfm = zeros(length(seqPssmMat),20);
%     a= zeros(length(seqPssmMat));
%     for i =1:length(seqPssmMat)
%         for j =1:20
%             pfm(i,j) =double(1/(1+exp(-seqPssmMat(i,j))));  %��һ��
%             a(i) = a(i)+pfm(i,j);   
%         end
%     end
%     for i = 1:length(seqPssmMat)
%         for j = 1:20
%             pfm(i,j) = double(pfm(i,j)*(1/a(i)));   %ʹÿһ�е����������Ϊ1
%         end
%     end
%     sum = 0;
%     for m = 1:20
%        for n = 1:20
%            for i = 1:length(seqPssmMat)-1   
%                sum = sum+pfm(i,m)*pfm(i+1,n);     
%            end
%            bigram(m,n)=sum;     %���bi-gram����һ��������
%         sum = 0;
%         end
%     end

    %[A_nm,zmlist,cidx,V_nm]=P_zernike(seqPssmMat,n,m);
%   feature = zeros(400,1);
%   count = 1;
%   for i=1:20
%       for j=1:20
%           feature = bigram(i,j);
%           count = count + 1;
    myMap(seqName) = seqPssmMat;
    %feature{i} = abs(A_nm);
    %cell(seqName)={seqPssmMat};
end
save('F:/RPI_SAN/PSSM/369.mat','myMap');
%save('F:/result/PSSM488.mat','cell');



