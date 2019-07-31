function [idx,central,number,cost]=samenumberkmeans(K,location)
%调试用
% location=100*rand(2,120);%输入的位置必须是2*x维
% clearvars -except location;
% clc
% K=6;
N=size(location,2);
R_up=ceil(N/K);%每组用户个数上限
R_down=fix(N/K);%每组用户个数下限
R=mod(N,K);%上限用户组的个数
stop=0;
A0=[];

% k-means++ 初始化中心
central = location(:,1+round(rand*(size(location,2)-1)));
idx = ones(1,size(location,2));
for i = 2:K
    D = location-central(:,idx);
    D = cumsum(sqrt(dot(D,D,1)));
    if D(end) == 0
        central(:,i:K) = location(:,ones(1,K-i+1));
        return;
    end
    central(:,i) = location(:,find(rand < D/D(end),1));
    [~,idx] = max(bsxfun(@minus,2*real(central'*location),dot(central,central,1).'));
end
D=[];

r=0;
cost=0;

while stop==0
    cost0=cost;
    central0=central;
    for n=1:N
        for k=1:K
            D(n,k)=norm(location(:,n)-central0(:,k));%计算点到中心距离
        end
    end
    [~,C]=sort(D,2);%按行排序
    
    %限制每组数量的初始化分组
    if r==0
        R0=0;%判断上限用户组数是否已满
        idx=zeros(N,1);
        C0=C;
        for n=1:N
            for k=1:K
                if idx(n)~=0
                    break
                end
                if C0(n,k)~=0
                    idx(n)=C0(n,k);
                else
                    continue
                end
            end
            if n==N
                break
            end
            if R0<R
                if length(find(idx==idx(n)))>=R_up%这组数量已满，则令后面用户属于最优中心的C为0
                    for i=(n+1):N
                        j=find(C0(i,:)==idx(n));
                        C0(i,j)=0;
                    end
                    R0=R0+1;
                    if R0==R%此时数量为R_down的中心也不需要再加入用户，之后的C也为0
                        for m=1:K
                            if length(find(idx==m))==R_down
                                for i=(n+1):N
                                    j=find(C0(i,:)==m);
                                    C0(i,j)=0;
                                end
                            end
                        end
                    end
                end
            else
                if length(find(idx==idx(n)))>=R_down%这组数量已满，则令后面用户属于最优中心的C为0
                    for i=(n+1):N
                        j=find(C0(i,:)==idx(n));
                        C0(i,j)=0;
                    end
                end
            end
        end
    end
    
    %置换操作and转移操作
    p=0;p0=-1;
    while p>p0
        p0=p;
        %置换操作
        for i=1:N
            a=zeros(1,N);
            for j=1:N
                if j~=i
                    d1=(D(i,idx(i)))^2+(D(j,idx(j)))^2;
                    d2=(D(i,idx(j)))^2+(D(j,idx(i)))^2;
                    if d2<d1
                        a(j)=d1/d2;%找到能和第i个用户换中心的用户
                    end
                end
            end
            if sum(a)>0
                [~,j]=max(a);%得到能换的最好的用户
                b=idx(i);
                idx(i)=idx(j);
                idx(j)=b;
                p=p+1;
            end
        end
        %转移操作(针对不能均分的情况）
        if mod(N,K)~=0
            for n=1:N
                if length(find(idx==idx(n)))==R_up
                    B=[];%少于最大数量的中心
                    for k=1:K
                        if length(find(idx==k))<R_up
                            B=cat(1,B,k);
                        end
                    end
                    if isempty(B)==1
                        break
                    end
                    [d,b]=min(D(n,B));%找到最好的中心
                    c=B(b);
                    if d<D(n,idx(n))
                        idx(n)=c;%赋值
                        p=p+1;
                    end
                else
                    continue;
                end
            end
        end
    end
    
    C=zeros(N,K);
    cost=0;
    for k=1:K
        c=find(idx==k);
        C(c,k)=ones;
        central(:,k)=mean(location(:,c),2);%重新计算中心
        for i=1:length(c)
            cost=cost+(norm(central(:,k)-location(:,c(i))))^2;%计算距离平方和
        end
    end
    if abs(cost-cost0)<1e-6
        stop=1;
    end
    r=r+1;
    A1=[r;cost];
    A0=cat(2,A0,A1);
end

%计算每组的数量
if mod(N,K)==0
    number=R_up*ones(1,K);
else
    for k=1:K
        number(1,k)=length(find(idx==k));
    end
end

% figure;
% plot(A0(1,:),A0(2,:));
% title('uniform kmeans cost');
% 
% color(1,:)=[0 0 0];color(2,:)=[1 0 0];color(3,:)=[0 1 0];color(4,:)=[0 0 1];color(5,:)=[0 1 1];color(6,:)=[1 0 1];color(7,:)=[1 1 0];
% figure;
% for i=1:K
%     j=i;
%     if j>7
%         j=j-7;
%     end
%     scatter(location(1,idx==i),location(2,idx==i),10,color(j,:),'filled');
%     hold on
% end
% plot(central(1,:),central(2,:),'ko',...
%     'MarkerSize',5,'LineWidth',3)
% title('uniform kmeans');
% hold off


