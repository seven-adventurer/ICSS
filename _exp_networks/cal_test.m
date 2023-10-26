clc;
clear;

net_id=1;

% Net=f_read_Net(net_id);
A = load(strcat('./', num2str(net_id))).Net;
edge = sum(A(:));
n= size(A);
n = n(1);
DC = sum(A);

adj_idx=zeros(n+2,1);
adj_idx(1)=n;
for i=2:1:n+1
    adj_idx(i+1)=adj_idx(i)+DC(i-1);
end

adj=[];
for i =1:1:n
    index = find(A(i,:));    
    adj=[adj,index];
end
adj = adj';

[value,method_id]=sort(DC,'descend');

% save in binary
path='./networks_data/';
file1 = strcat(path,'net_adj_idx_',num2str(net_id),'.txt');
fid=fopen(file1, 'wb');
fwrite(fid,adj_idx,'int');
fclose(fid);

file2 = strcat(path,'net_adj_',num2str(net_id),'.txt');
fid=fopen(file2, 'wb');
fwrite(fid,adj,'int');
fclose(fid);

file3 = strcat(path,'method_',num2str(net_id),'.txt');
fid=fopen(file3, 'wb');
fwrite(fid,method_id','int');
fclose(fid);

% net_id=3;
% path='./networks_data/';
% file1 = strcat(path,'net_adj_idx_',num2str(net_id),'.txt');
% fid=fopen(file1, 'r');
% [meth, count]=fread(fid,10000,'int');
% fclose(fid);
