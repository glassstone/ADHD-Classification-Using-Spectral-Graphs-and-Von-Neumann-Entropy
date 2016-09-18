clear all
dir = 'C:\Users\Documents\ADHD_data\TimeSeriesData\avgTimeSeries';
[num, txt, raw] = xlsread('C:\Users\Documents\ADHD_data\ALL_phenotypic.xls') ;
k_min = 0
k_max = 350 ;
k = 4 ;
m=10 ;
mc = 3 ;
a =0.1 ;
Evalues= zeros(621, 14);
i = 1 ;
while i < 621
    i = i+ 1
   row = raw(i,1) ;
   mat = cell2mat(row) ;
   id = int2str(mat) ;
   if length(id) == 5
       id = strcat('00', id);
   end
   fname = strcat(dir,id) ;
   fname = strcat(fname,'.mat');
   if exist(fname, 'file') ~= 2
       continue
   end


   row = raw(i,2) ;
   site = cell2mat(row) ;

   row = raw(i,3) ;
   gen = cell2mat(row) ;

   row = raw(i,4) ;
   age = cell2mat(row) ;

   row = raw(i,5) ;
   hand = cell2mat(row) ;


   row = raw(i,6) ;
   dx = cell2mat(row) ;

   row = raw(i,8) ;
   adhdix = cell2mat(row) ;
   if isnumeric(adhdix) == 0
       adhdix = 0 ;
   end
   row = raw(i,9) ;
   inatt = cell2mat(row) ;
   if isnumeric(inatt) == 0
       inatt = 0 ;
   end

   row = raw(i,10) ;
   hyper = cell2mat(row) ;
   if isnumeric(hyper) == 0
       hyper = 0 ;
   end


   %%%%%%%%%%%
   load(fname) ;
   data = avgTimeSeriesFileMeanSubtract ;

   A = zeros(length(data(:,1)));

   for ii=1:1:length(data(:,1))
        idx=knnsearch(data,data(ii,:),'K',m,'Distance','cosine');
        for jj=1:1:m
            A(ii,idx(jj))=dot(data(ii,:),data(idx(jj),:))/(norm(data(ii,:),2)*
            norm(data(idx(jj),:),2));
            if((A(ii,idx(jj)) == 0) || (isnan(A(ii,idx(jj)))))
                A(ii,idx(jj)) = a;
            end
            A(idx(jj),ii)=A(ii,idx(jj));
        end
   end


   D=zeros(length(data(:,1)));

   for ii=1:1:length(data(:,1))
      D(ii,ii)=sum(A(ii,:));
   end

   %%%%%%%%%%%%


   L= D-A;
   L=D^(-1/2)*L*D^(-1/2);
%%%%%%%%%%%%% +++++++++++++++++++++++++
[U,V]= eig(L);

v=diag(V);

[vs, is] = sort(v,'ascend');

u=[];

for i1=1:1:k
    u = [u U(:,is(i1))];
end

% we re-normalize the matrix rows to unit norm

% for i1=1:1:length(u(:,1))
%     u(i1,:)=u(i1,:)/sqrt(sum(u(i,:).*u(i1,:)));
% end


% [membership, ctrs, sumd, pointd] = kmeans(input,...
%          clusters,'Replicates',num_replicates,...
%          'Distance','sqEuclidean');
%
% we are using k-means.  You can use any other clustering algorithm
% we also use the parameter replicates in order to obtain the "steady
% state" solution for k-means

[labels, C_spec, sumd] = kmeans(u,k,'Replicates',10);

% figure('Name', strcat('Figure: ', int2str(i)))
% scatter(u(:,3),u(:,4),50, labels,'filled')

 wc1 = data(labels==1,:)  ;
 wc2 = data(labels==2,:)  ;
 wc3 = data(labels==3,:)  ;
 wc4 = data(labels==4,:)  ;

mentwc1 = entropy(wc1) ;
mentwc2 = entropy(wc2) ;
mentwc3 = entropy(wc3) ;
mentwc4 = entropy(wc4) ;


%%%%**************************

ac1 = zeros(length(wc1(:,1)));

   for ii=1:1:length(wc1(:,1))
        idx=knnsearch(wc1,wc1(ii,:),'K',mc,'Distance','cosine');
        for jj=1:1:mc
            ac1(ii,idx(jj))=dot(wc1(ii,:),wc1(idx(jj),:))/(norm(wc1(ii,:),2)
            * norm(wc1(idx(jj),:),2));
            if((ac1(ii,idx(jj)) == 0) || (isnan(ac1(ii,idx(jj)))))
                ac1(ii,idx(jj)) = a;
            end
            ac1(idx(jj),ii)=ac1(ii,idx(jj));
        end
   end

   dc1=zeros(length(wc1(:,1)));

   for ii=1:1:length(wc1(:,1))
      dc1(ii,ii)=sum(ac1(ii,:));
   end

   lc1= dc1-ac1;
   lc1=dc1^(-1/2)*lc1*dc1^(-1/2);

   evc1 = eig(lc1) ;
   maxevc1 = max(evc1) ;
   av_degc1 = trace(dc1) ;
   entc1 = -evc1' * log2(evc1+(evc1==0)) ;

%%%%****************************

%%%%**************************

ac2 = zeros(length(wc2(:,1)));

   for ii=1:1:length(wc2(:,1))
        idx=knnsearch(wc2,wc2(ii,:),'K',mc,'Distance','cosine');
        for jj=1:1:mc
            ac2(ii,idx(jj))=dot(wc2(ii,:),wc2(idx(jj),:))/(norm(wc2(ii,:),2)
            * norm(wc2(idx(jj),:),2));
            if((ac2(ii,idx(jj)) == 0) || (isnan(ac2(ii,idx(jj)))))
                ac2(ii,idx(jj)) = a;
            end
            ac2(idx(jj),ii)=ac2(ii,idx(jj));
        end
   end

   dc2=zeros(length(wc2(:,1)));

   for ii=1:1:length(wc2(:,1))
      dc2(ii,ii)=sum(ac2(ii,:));
   end

   lc2= dc2-ac2;
   lc2=dc2^(-1/2)*lc2*dc2^(-1/2);

   evc2 = eig(lc2) ;
   maxevc2 = max(evc2) ;
   av_degc2 = trace(dc2) ;

   entc2 = real(-evc2' * log2(evc2+(evc2==0)) );

%%%%****************************

%%%%**************************

ac3 = zeros(length(wc3(:,1)));

   for ii=1:1:length(wc3(:,1))
        idx=knnsearch(wc3,wc3(ii,:),'K',mc,'Distance','cosine');
        for jj=1:1:mc
            ac3(ii,idx(jj))=dot(wc3(ii,:),wc3(idx(jj),:))/(norm(wc3(ii,:),2) *
            norm(wc3(idx(jj),:),2));
            if((ac3(ii,idx(jj)) == 0) || (isnan(ac3(ii,idx(jj)))))
                ac3(ii,idx(jj)) = a;
            end
            ac3(idx(jj),ii)=ac3(ii,idx(jj));
        end
   end

   dc3=zeros(length(wc3(:,1)));

   for ii=1:1:length(wc3(:,1))
      dc3(ii,ii)=sum(ac3(ii,:));
   end

   lc3= dc3-ac3;
   lc3=dc3^(-1/2)*lc3*dc3^(-1/2);

   evc3 = eig(lc3) ;
   maxevc3 = max(evc3) ;
   av_degc3 = trace(dc3) ;

   entc3 = real(-evc3' * log2(evc3+(evc3==0)) );

%%%%****************************

%%%%**************************

ac4 = zeros(length(wc4(:,1)));

   for ii=1:1:length(wc4(:,1))
        idx=knnsearch(wc4,wc4(ii,:),'K',mc,'Distance','cosine');
        for jj=1:1:mc
            ac4(ii,idx(jj))=dot(wc4(ii,:),wc4(idx(jj),:))/(norm(wc4(ii,:),2)
            * norm(wc4(idx(jj),:),2));
            if((ac4(ii,idx(jj)) == 0) || (isnan(ac4(ii,idx(jj)))))
                ac4(ii,idx(jj)) = a;
            end
            ac4(idx(jj),ii)=ac4(ii,idx(jj));
        end
   end

   dc4=zeros(length(wc4(:,1)));

   for ii=1:1:length(wc4(:,1))
      dc4(ii,ii)=sum(ac4(ii,:));
   end

   lc4= dc4-ac4;
   lc4=dc4^(-1/2)*lc4*dc4^(-1/2);

   evc4 = eig(lc4) ;
   maxevc4 = max(evc4) ;
   av_degc4 = trace(dc4) ;

   entc4 = real(-evc4' * log2(evc4+(evc4==0)) );

%%%%****************************
%%%%%%%%%%%%%% +++++++++++++++++++
%    ev = eig(L) ;
%    maxev = max(ev) ;
%    av_deg = trace(D) ;
%    ent = -ev' * log2(ev+(ev==0)) ;
%    ment = entropy(ev)  ; % entropy using native function
%    dataent = entropy(data) ;
   Evalues(i,1) = mat ;  %id
   Evalues(i,2)= site ; %site loc
   Evalues(i,3) = dx ;  %adhd versus td
   Evalues(i,4) = gen ; %gender
   Evalues(i,5) = age ;


   Evalues(i,6) = hand ; %handedness
   Evalues(i,7) = adhdix; % adhd index
   Evalues(i,8) = inatt ; % inattentive
   Evalues(i,9) = hyper ; % hyperactive

   Evalues(i,10) = entc1 ; %vn entropy
   Evalues(i,11) = entc2 ; % std entropy
   Evalues(i,12) = entc3 ; % max evalue
   Evalues(i,13) = entc4 ; % av degree
%    cluster1 = (labels == 1) ;
%    cluster2 = (labels == 2) ;
%    cluster3 = (labels == 3) ;
%    cluster4 = (labels == 4) ;
   figure('Name', strcat('Figure: ', int2str(i)))
   scatter(u(:,2),u(:,3),50, labels,'filled' )

   xlabel('Data Dimension  2') ;
   ylabel('Data Dimension  3') ;

   if dx > 0
       Evalues (i, 14) = 1 ;  % ADHD
        title('ADHD Subject') ;
   else
       Evalues(i,14) = 0 ;  % TD
        title('TD Subject') ;
   end



   if i > 20
       break ;
   end

end
xlswrite('C:\Users\Documents\ADHD_data\result.xls',Evalues)







%plot(vs,'.')



%Nsize = size(A,1);

%coords = [cos(2*pi*(1:Nsize)/Nsize); sin(2*pi*(1:Nsize)/Nsize)]';
%gplot(A, coords)
%text(coords(:,1) - 0.1, coords(:,2) + 0.1, num2str((1:Nsize)'), 'FontSize', 14)

%imagesc(A);
%colormap('gray');colorbar ;
