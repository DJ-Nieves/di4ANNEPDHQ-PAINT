%% Lr_vs_Lcross_w_Pearsons v1.0
% 
% This is used to calculate the L(r) vs Lcross(r) (according to Rossy et al., (2014) Histochem Cell Biol 141, 605-612) 
% from "marked" single molecule localisation data, e.g., here each point is
% associated with a value for membrane polarity. This code was implemented in 
% "Quantitative mapping of nanoenvironments through single-molecule imaging
% of solvatochromic probes", Nieves et al., (2020), bioarxiv, doi:
% 
% The L(r) vs Lcross(r) is calculated for the GP segregated channels, and
% also for points randomly assigned into the channels as a control. Finally
% the Pearson's correlation coefficient for the L(r) vs Lcross(r) data is
% calculated for each channel.
% 
% Author: Daniel J. Nieves, University of Birmingham
% email: d.j.nieves@bham.ac.uk
% 
% data to provide: 
% - a matrix of the form (x,y, GP value), called "GP_data".

%## Inputs and Outputs:

% Inputs:
% 
% - r :	         Radius of capture for Lr and Lcross (nm)
% 
% - GP_co:       GP value to use split the data. A vector using a series of cut-off values may also be used.
% 
% - pos:         Coordinates of the ROI from which to take [x,y,GP] data.
% 

% 
% Outputs:
% 
% - Lr_plus :       list of L(r) values for each point in the above GP cut-off
%                   data.
% - Lr_minus :      list of L(r) values for each point in the above GP cut-off
%                   data.
% - Lr_rand_plus :  list of L(r) values for each point assigned
%                   randomly in the above GP cut-off data.
% - Lr_rand_minus : list of L(r) values for each point assigned
%                   randomly in the below GP cut-off data.
% 
% - Lc_plus :       list of Lcross(r) values for each point in the above GP cut-off
%                   data.
% - Lc_minus :      list of Lcross(r) values for each point in the below GP cut-off
%                   data.
% - Lc_rand_plus :  list of Lcross(r) values for each point assigned
%                   randomly in the above GP cut-off data.
% - Lc_rand_minus : list of Lcross(r) values for each point assigned
%                   randomly in the below GP cut-off data.
% 
% - pLrLc_plus :            Pearson's coefficent of L(r) vs Lcross(r) for above GP cut-off data
% - pLrLc_minus :           Pearson's coefficent of L(r) vs Lcross(r) for below GP cut-off data
% - pLrLc_rand_plus :       Pearson's coefficent of L(r) vs Lcross(r) for randomly assigned above GP cut-off data
% - pLrLc_rand_minus :      Pearson's coefficent of L(r) vs Lcross(r) for randomly assigned below GP cut-off data
% 

%% GP_data must be loaded
% Data must be in the form of a 3 column matrix with (x, y, GP value), and
% called GP_data. x and y here are in nm

% user defined parameters for radius of search (r)
% GP cut off value (GP_co) and ROI (pos)

% radius of capture for Lr and Lcross (nm)
r = 100; 

% GP value to use split the data. 
%A vector using a series of cut-off values may also be used.
GP_co = 0.2;

%coordinates of the ROI to be examined (nm)
pos = [10000,15000; 10000,25000; 20000,25000; 20000,15000];
%% segregate data using GP cut offs
% pre-allocation
GP_seg_plus = cell(length(GP_co),1,1);
GP_seg_minus = cell(length(GP_co),1,1);
rand_seg_plus = cell(length(GP_co),1,1);
rand_seg_minus = cell(length(GP_co),1,1);

%fill the cell array with the GP segregations
for i = 1:length(GP_co)
    plus = GP_data(:,3) > GP_co(1,i);
    minus = GP_data(:,3) < GP_co(1,i);
    GP_seg_plus{i,1} = GP_data(plus==1,:);
    GP_seg_minus{i,1} = GP_data(minus==1,:);
end

%fill the cell array with the random control segregations
for i = 1:length(GP_seg_plus)
    p = randperm(length(GP_data));
    rand_seg_plus{i,1} = GP_data(p(1,1:length(GP_seg_plus{i,1})),1:3);
    rand_seg_minus{i,1} = GP_data(p(1,length(GP_seg_plus{i,1}):length(GP_seg_minus{i,1})+length(GP_seg_plus{i,1})),1:3);
end
%% Selection of x,y points within the region of interest
% pre-allocation
GPP_roi = cell(length(GP_seg_plus),1,1);
GPM_roi = cell(length(GP_seg_minus),1,1);
GPP_rand_roi = cell(length(GP_seg_plus),1,1);
GPM_rand_roi = cell(length(GP_seg_minus),1,1);

% now points from each segregated channel within the ROI are selected and
% stored
for jj = 1:length(GP_seg_plus)
kpin = inpolygon(GP_seg_plus{jj,1}(:,1),GP_seg_plus{jj,1}(:,2),pos(:,1), pos(:,2));
GPP_roi{jj,1} = GP_seg_plus{jj,1}(kpin==1,:);
end

for jj = 1:length(GP_seg_minus)
kmin = inpolygon(GP_seg_minus{jj,1}(:,1),GP_seg_minus{jj,1}(:,2),pos(:,1), pos(:,2));
GPM_roi{jj,1} = GP_seg_minus{jj,1}(kmin==1,:);
end

for jj = 1:length(rand_seg_plus)
kpin = inpolygon(rand_seg_plus{jj,1}(:,1),rand_seg_plus{jj,1}(:,2),pos(:,1), pos(:,2));
GPP_rand_roi{jj,1} = rand_seg_plus{jj,1}(kpin==1,:);
end

for jj = 1:length(rand_seg_minus)
kmin = inpolygon(rand_seg_minus{jj,1}(:,1),rand_seg_minus{jj,1}(:,2),pos(:,1), pos(:,2));
GPM_rand_roi{jj,1} = rand_seg_minus{jj,1}(kmin==1,:);
end
%% pre-allocation for L(r) and Lcross(r) above GP loops
Lr_plus = cell(length(GP_co),1,1);
for i = 1:length(GPP_roi)
    Lr_plus{i,1} = zeros(1,length(GPP_roi{i,1}));
end

Lc_plus = cell(length(GP_co),1,1);
for i = 1:length(GPP_roi)
    Lc_plus{i,1} = zeros(1,length(GPP_roi{i,1}));
end

Lr_rand_plus = cell(length(GP_co),1,1);
for i = 1:length(GPP_roi)
    Lr_rand_plus{i,1} = zeros(1,length(GPP_roi{i,1}));
end

Lc_rand_plus = cell(length(GP_co),1,1);
for i = 1:length(GPP_roi)
    Lc_rand_plus{i,1} = zeros(1,length(GPP_roi{i,1}));
end
%% L(r) and Lcross(r) loop, with Pearson's Correlation, for above GP cut-off data (GPP_roi)

for ii = 1:length(GPP_roi)
kp = GPP_roi{ii,1}(:,1:2);
km = GPM_roi{ii,1}(:,1:2);
L50p = zeros(1,length(GPP_roi{ii,1}));
p_counts = zeros(1,length(GPP_roi{ii,1}));
m_counts = zeros(1,length(GPM_roi{ii,1}));
LC50p = zeros(1,length(GPP_roi{ii,1}));
tic
    for i = 1:length(kp)
        for j = 1:length(kp)
            distances = sqrt((kp(j,1) - kp(i,1)).^2 + (kp(j,2) - kp(i,2)).^2);
            p_counts(j) = distances < r;
        end
        L50p(i) = sum(p_counts)-1;
    end
toc
Lr_plus{ii,1} = L50p;

tic
    for i = 1:length(kp)
        for j = 1:length(km)
            distances = sqrt((km(j,1) - kp(i,1)).^2 + (km(j,2) - kp(i,2)).^2);
            m_counts(j) = distances < r;
        end
        LC50p(i) = sum(m_counts);
    end
toc
Lc_plus{ii,1} = LC50p;
end

% Pearson's for above cut-off Lr vs Lcross
pLrLc_plus = zeros(1,length(Lr_plus));
for cc = 1:length(Lr_plus)
comb = [Lr_plus{cc,1}.' Lc_plus{cc,1}.'];
[mu,ic,ia] = unique(comb, 'rows', 'stable');
[N,edges] = histcounts(ia,max(ia));
Cplus = [mu N.'];
rCplus = corrcoef(Cplus(:,1),Cplus(:,2));
pLrLc_plus(cc) = rCplus(2,1);
end 

%loops for data assigned randomly to the above GP cut-off group (control)
%L(r) and Lcross(r)
for ii = 1:length(GPP_rand_roi)
kp = GPP_rand_roi{ii,1}(:,1:2);
km = GPM_rand_roi{ii,1}(:,1:2);
L50p = zeros(1,length(GPP_rand_roi{ii,1}));
p_counts = zeros(1,length(GPP_rand_roi{ii,1}));
m_counts = zeros(1,length(GPM_rand_roi{ii,1}));
LC50p = zeros(1,length(GPP_rand_roi{ii,1}));
tic
    for i = 1:length(kp)
        for j = 1:length(kp)
            distances = sqrt((kp(j,1) - kp(i,1)).^2 + (kp(j,2) - kp(i,2)).^2);
            p_counts(j) = distances < r;
        end
        L50p(i) = sum(p_counts)-1;
    end
toc
Lr_rand_plus{ii,1} = L50p;

tic
    for i = 1:length(kp)
        for j = 1:length(km)
            distances = sqrt((km(j,1) - kp(i,1)).^2 + (km(j,2) - kp(i,2)).^2);
            m_counts(j) = distances < r;
        end
        LC50p(i) = sum(m_counts);
    end
toc
Lc_rand_plus{ii,1} = LC50p;
end

% Pearson's for above cut-off, randomly assigned, Lr vs Lcross
pLrLc_rand_plus = zeros(1,length(Lr_plus));
for cc = 1:length(Lr_rand_plus)
comb = [Lr_rand_plus{cc,1}.' Lc_rand_plus{cc,1}.'];
[mu,ic,ia] = unique(comb, 'rows', 'stable');
[N,edges] = histcounts(ia,max(ia));
Cplus = [mu N.'];
rCplus = corrcoef(Cplus(:,1),Cplus(:,2));
pLrLc_rand_plus(cc) = rCplus(2,1);
end 

%% pre-allocation for L(r) and Lcross(r) for below cut-off GP data
Lr_minus = cell(length(GP_co),1,1);
for i = 1:length(GPM_roi)
    Lr_minus{i,1} = zeros(1,length(GPM_roi{i,1}));
end

Lc_minus = cell(length(GP_co),1,1);
for i = 1:length(GPM_roi)
    Lc_minus{i,1} = zeros(1,length(GPM_roi{i,1}));
end

Lr_rand_minus = cell(length(GP_co),1,1);
for i = 1:length(GPM_roi)
    Lr_rand_minus{i,1} = zeros(1,length(GPM_rand_roi{i,1}));
end

Lc_rand_minus = cell(length(GP_co),1,1);
for i = 1:length(GPM_roi)
    Lc_rand_minus{i,1} = zeros(1,length(GPM_rand_roi{i,1}));
end
%% L(r) and Lcross(r) loop, with Pearson's Pearson's Correlation, for below GP cutoff data (GPM_roi)

for ii = 1:length(GPM_roi)
kp = GPP_roi{ii,1}(:,1:2);
km = GPM_roi{ii,1}(:,1:2);
L50m = zeros(1,length(GPM_roi{ii,1}));
p_counts = zeros(1,length(GPP_roi{ii,1}));
m_counts = zeros(1,length(GPM_roi{ii,1}));
LC50m = zeros(1,length(GPM_roi{ii,1}));

tic
    for i = 1:length(km)
        for j = 1:length(km)
            distances = sqrt((km(j,1) - km(i,1)).^2 + (km(j,2) - km(i,2)).^2);
            m_counts(j) = distances < r;
        end
        L50m(i) = sum(m_counts)-1;
    end
toc
Lr_minus{ii,1} = L50m;

tic
    for i = 1:length(km)
        for j = 1:length(kp)
            distances = sqrt((kp(j,1) - km(i,1)).^2 + (kp(j,2) - km(i,2)).^2);
            p_counts(j) = distances < r;
        end
        LC50m(i) = sum(p_counts);
    end
toc
Lc_minus{ii,1} = LC50m;
end

% Pearson's for below cut off Lr vs Lcross
pLrLc_minus = zeros(1,length(Lr_minus));
for cc = 1:length(Lr_minus)
comb = [Lr_minus{cc,1}.' Lc_minus{cc,1}.'];
[mu,ic,ia] = unique(comb, 'rows', 'stable');
[N,edges] = histcounts(ia,max(ia));
Cminus = [mu N.'];
rCminus = corrcoef(Cminus(:,1),Cminus(:,2));
pLrLc_minus(cc) = rCminus(2,1);
end 

%loops for data assigned randomly to the below GP cut-off group (control)
%L(r) and Lcross(r)
for ii = 1:length(GPM_rand_roi)
kp = GPP_rand_roi{ii,1}(:,1:2);
km = GPM_rand_roi{ii,1}(:,1:2);
L50m = zeros(1,length(GPM_rand_roi{ii,1}));
p_counts = zeros(1,length(GPP_rand_roi{ii,1}));
m_counts = zeros(1,length(GPM_rand_roi{ii,1}));
LC50m = zeros(1,length(GPM_rand_roi{ii,1}));
tic
    for i = 1:length(km)
        for j = 1:length(km)
            distances = sqrt((km(j,1) - km(i,1)).^2 + (km(j,2) - km(i,2)).^2);
            m_counts(j) = distances < r;
        end
        L50m(i) = sum(m_counts)-1;
    end
toc
Lr_rand_minus{ii,1} = L50m;

tic
    for i = 1:length(km)
        for j = 1:length(kp)
            distances = sqrt((kp(j,1) - km(i,1)).^2 + (kp(j,2) - km(i,2)).^2);
            p_counts(j) = distances < r;
        end
        LC50m(i) = sum(p_counts);
    end
toc
Lc_rand_minus{ii,1} = LC50m;
end

% Pearson's for below cut off, randomly assigned, Lr vs Lcross
pLrLc_rand_minus = zeros(1,length(Lr_rand_minus));
for cc = 1:length(Lr_rand_minus)
comb = [Lr_rand_minus{cc,1}.' Lc_rand_minus{cc,1}.'];
[mu,ic,ia] = unique(comb, 'rows', 'stable');
[N,edges] = histcounts(ia,max(ia));
Cminus = [mu N.'];
rCminus = corrcoef(Cminus(:,1),Cminus(:,2));
pLrLc_rand_minus(cc) = rCminus(2,1);
end 