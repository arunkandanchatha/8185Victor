clear all;
clc;
A=xlsread('GDPDataV2',6,'E2:H259');
[n,k]=size(A);
w=1600; %the smoothing parameter for quarterly data%
quarter=(1947:1/4:2011.25)';

%-----------------------------------------------------------------
%Q1B
GDP=log((A(:,1)));

%HP filtering log of GDP%
[hptrend, hpcyc] = hpfilter([GDP A(:,2:4)]);
GDP_hp=hptrend(:,1);

figure
plot(quarter,[GDP GDP_hp]);
title('(Log of) Real GDP in billions of dollars (2005=100) and HP filtered series');
ylabel('');
xlim([1947  2011]);
legend('(Log of) GDP','HP filtered series', 'Location','SouthEast');

%now create cooley table
c = xcorr(hpcyc, 5, 'coeff');

%how is c setup now? As follows:
%  rows: correlation at each lag. Row 1 is x(-5), row 2 is x(-4), etc
%  columns: for a specific row (i.e. specific lag), the columns are
%           1: correlation of GDP with itself
%           2: correlation of GDP with investments
%           3: correlation of GDP with consumption
%           3: correlation of GDP with aggregate hours
%           4: correlation of investments with GDP
%           5: etc.
%  we want the correlation of GDP with four variables, for all lags
disp('   x(-5)      x(-4)     x(-3)     x(-2)     x(-1)      x         x(1)      x(2)       x(3)      x(4)    x(5)')
cooleytable = c(:,1:4)'

%-----------------------------------------------------------------
%Q1C, 1D
% we assume a linear realationship (a+tb) between gdp and time, where time
% is the number of periods from our start date to the data point in
% question.

%essentially, the time index
t=1:1:n;t=t';

%linear, so a constant and the time index
x=[ones(n,1) t];
y=GDP;

% run least squares regression
b=inv(x'*x)*x'*y;

% our trend line is basically the b(1)+b(2)*t
trend=b(1)+b(2)*t;

%now we find calculate the hp filter and the residuals
eGDP=GDP-trend;
[eGDP_hp, eGDP_res]=hpfilter(eGDP,w);

%let's plot gdp and each of the components
figure
quarter=(1947:1/4:2011.25);
plot(quarter,[GDP trend]);
title('GDP, Linear trend')
ylabel('')
xlim([1947  2011])
legend('RGDP', 'Trend', 'Location','SouthEast')

figure
quarter=(1947:1/4:2011.25);
plot(quarter,[eGDP_hp eGDP_res]);
title('GDP de-trended H-P Filter')
ylabel('')
xlim([1947  2011])
legend('H-P trend', 'Residuals', 'Location','SouthEast')

%growth rate is simply log(GDP)-log(GDP(-1))
GDP_growth=diff(GDP);

%let's plot growth versus residual
figure
quarter=(1947:1/4:2011.25);
plot(quarter(2:end),[GDP_growth eGDP_res(2:end)]);
title('Growth, HP Filtered residuals')
ylabel('')
xlim([1947  2011])
legend('(Log difference of) GDP', 'Residuals', 'Location','SouthEast')

% Q1E
% the original data has trends, so we need to detrend it and rescale
% so that everything is at the same scale
% but HP_cyc is exactly the detrended data using HP filtering
%logDataDiff = hpcyc;
%qend=size(quarter,2);
%       Aside: There seems to be questions on whether using an HP filter
%              is appropriate. Instead, I'm going to use a simple log
%              detrending
logDataDiff=diff(log(A));
qend=size(quarter,2)-1;

%impulse response is sensitive to order. So let's reorder from least
%exogenous to most exogenous
titleArray = {'Investment', 'Labour', 'Consumption', 'RGDP'};
logDataDiff=[logDataDiff(:,2) logDataDiff(:,4) logDataDiff(:,3) ...
    logDataDiff(:,1)];
%titleArray = {'RGDP', 'Investment', 'Consumption', 'Labour'};

%   The below code lets us see how the scaling is of these variables. We
%   can see that investment is an order of magnitude larger than everything
%   else, so we rescale
% figure
% subplot(4,1,1)
% plot(quarter(1:qend),logDataDiff(:,1),'r');
% title('GDP')
% hold('on')
% subplot(4,1,2);
% plot(quarter(1:qend),logDataDiff(:,2),'b');
% title('Inv')
% subplot(4,1,3);
% plot(quarter(1:qend),logDataDiff(:,3),'k');
% title('Consumption')
% subplot(4,1,4);
% plot(quarter(1:qend),logDataDiff(:,4),'k');
% title('Agg Hours')
% hold('off')

%everything looks okay. Don't rescale
factor = repmat([1 1 1 1],size(logDataDiff,1),1);
logDataDiff=factor.*logDataDiff;


% let's try a couple of models to see which fits best
%VAR5full = vgxset('n',4,'nAR', 5);
for i=1:5
%    dt=logical(triu(ones(4,4)));
%    dt=logical(triu(eye(4,4)));
%     VARarray(i) = vgxset('ARsolve',repmat({dt},i,1),'n',4);
    VARarray(i) = vgxset('n',4,'nAR', i);
end

%specify presample and estimation 
preInt=5;
estInt=ceil(0.95*size(logDataDiff,1));
dataPre=logDataDiff(1:preInt, :);
dataEst=logDataDiff(preInt+1:estInt,:);
dataFit=logDataDiff(estInt+1:end,:);

%define the model
for i=1:5
%   [EstSpec1,EstStdErrors1,LLF1,W1] = vgxvarx(VAR5full,dataEst,[],dataPre);
%    EstSpec1 = vgxset(EstSpec1,'Series', {'RGDP', 'Investment', 'Consumption', ...
%        'Labour'}); 
    temp=VARarray(i);
    [estspec(i),eststderror(i),llf(i),w] = vgxvarx(temp,...
                                            dataEst,[],dataPre);
    estspec(i) = vgxset(estspec(i),'Series', titleArray); 
end

[n1,n1p]=vgxcount(estspec(1));
[n2,n2p]=vgxcount(estspec(2));
[n3,n3p]=vgxcount(estspec(3));
[n4,n4p]=vgxcount(estspec(4));
[n5,n5p]=vgxcount(estspec(5));

whichOne = aicbic(llf, [n1p n2p n3p n4p n5p]);
[blah,minIndex]=min(whichOne);
EstSpec1 = estspec(minIndex);

%now we will plot the impulse response. We get the impulse from the
%Cholesky decomposition of Q in estspec
decomp = chol(EstSpec1.Q);

%lets plot response for next 20 periods (4 quarters times 5 years)
% w will represent impulses. 3 dimensional arrays are weird in matlab
% A(row,col,num) means 'num' arrays of dimension 'row'-by-'col'
%          w(:,:,1) is no impulse. 
%          w(:,:,2) is impulse on gdp
%          w(:,:,5) is impulse on investment
%          w(:,:,3) is impulse on consumption
%          w(:.:,4) is impulse on aggregate hours
numperiods=21;
w = zeros(numperiods,4,5);
impulseResp = zeros(numperiods,4,5);
realImpulse = zeros(numperiods,4,5);

% reldiff will be as follows
%          reldiff(:,:,1): impulse response to gdp
%          reldiff(:,:,2): impulse response to investment
%          etc
relDiff = zeros(numperiods,4,4);

%below loop but for the "no response" case
tempArray(:,:)=w(:,:,1);
impulseResp(:,:,1) = vgxproc(EstSpec1,tempArray,[],logDataDiff);
factor=factor(1:numperiods,:);
impulseResp(:,:,1)=factor.^-1.*impulseResp(:,:,1);
realImpulse(:,:,1) = exp(cumsum(impulseResp(:,:,1)));

%calculating the response to GDP impulse, investment, etc.
% for i=2, we find the decomposition of 1,1 which is the impulse of
% gdp. i=3 is investment, etc.
for i=2:5
    w(1,i-1,i) = decomp(i-1,i-1);
    
    %find the impulse response
    tempArray(:,:)=w(:,:,i);
    impulseResp(:,:,i) = vgxproc(EstSpec1,tempArray,[],logDataDiff);
    
    %undo scaling
    % recall, we divided investment by 10
    impulseResp(:,:,i)=factor.^-1.*impulseResp(:,:,i);
    realImpulse(:,:,i) = exp(cumsum(impulseResp(:,:,i)));
    
    %difference between impulse and no impulse
    temp1(:,:)=realImpulse(:,:,i);
    temp2(:,:)=realImpulse(:,:,1);
    temp3(:,:)=(temp1-temp2)./temp1;
    relDiff(:,:,i-1) = temp3;
end


figure
projection=(0:1/4:(numperiods-1)/4);
% for plotting, rows are as follows:
%     1: plot of each result to GDP shock
%     2: plot of each result to Investment shock
% etc
for i=1:4
    % i indexes the impulse. So 1: GDP, 2: Investment, etc.
    % j indexes the actual output.
    % so, (i,j)=(1,1) is GDP response to GDP shock
    %     (i,j)=(2,3) is Consumption response to Investment shock
    for j=1:4
        subplot(4,4,(i-1)*4+j)
        plot(projection,100*relDiff(:,j,i))
        title(strcat(titleArray(i),'-on-',titleArray(j)));
        xlabel('Years')
        ylabel('% change')
        xlim([0 (numperiods-1)/4]);
    end
end
