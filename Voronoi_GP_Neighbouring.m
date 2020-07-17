%% Voronoi_GP_Neighbouring v1.0

% This code is used to perform Voronoi tesselation analysis of "marked" single
% molecule localisation data, e.g., here each point/tile is associated with a
% value for membrane polarity. This code was implemented in "Quantitative
% mapping of nanoenvironments through single-molecule imaging of
% solvatochromic probes", Nieves et al., (2020),bioarxiv, doi:

% Here, the Voronoi tesselation of the marked data (in the form [x,y,GP
% value]) is calculated, and then neighbours withina certain GP range (rg)
% are found. If the tiles neighbour eachother, and are within the GP range,
% they are combined. This proceeds until the neighbouring converges and no
% new neighbours within the GP range can be made.

%This can also be compared to a random neighbouring situation by setting
%rand GP equal to 1.

% Author: Daniel J. Nieves, University of Birmingham
% email: d.j.nieves@bham.ac.uk

%% NOTE: this code requires the "uniquecell" function to run completely %%
% a copy will be provided at the repository with this code, or it can be
% downloaded from here: 
% https://uk.mathworks.com/matlabcentral/fileexchange/31718-unique-elements-in-cell-array
% make sure it is added to the path, or in the directory. 
%% Inputs and Outputs

% Data to provide: 
% - a matrix of the form (x,y, GP value), called "GP_data".

% Inputs (set in next section):

% - randGP:         Set to either 0 or 1. Default is 0. Choose 1 to run
%                   analysis again, with randomised neighbouring for
%                   comparison

% - rg:             This sets the +/- range of GP values to include in
%                   the neighbouring. E.g. a value of 0.1 would mean values
%                   0.1 either side of a GP value will be included in the
%                   neighbouring group.

% - area_co:        This sets the maximum area (nm^2) to allow into the Voronoi
%                   neighbouring, as larger tiles maybe from areas of low data 
%                   coverage.

% - pos:            Coordinates of the ROI from which to take [x,y,GP]
%                   data.

% Outputs:

% - GP_v:           Output from Voronoi tesselation of GP_data.
%                   *GP_v{:,1} - coordinates of tile vertices.
%                   *GP_v{:,2} - GP value.
%                   *GP_v{:,3} - Tile Area.
%                   *GP_v{:,4} - Indices of tiles within +/- rg of GP value.
%                   *GP_v{:,5} - Indices of the initial neighbouring result.

% - compGPv:        Final indices from complete neighbouring of Voronoi tiles by
%                   GP value.

% - compMeanGP:     Final mean GP value for combined tiles.

% - compGPareas:    Final combined tile areas.


%% GP_data must be loaded
% Data must be in the form of a 3 column matrix with (x, y, GP value), and
% called GP_data. x and y here are in nm

% user defined parameters for randomised neighbouring (randGP; equal to 0
% by default, but can be made equal to 1 for random nieghbouring
% comparison), range of GP values to neighbour (rg; will generate +/- the
% value) and the region to be analysed (pos; coordinates of region of
% interest in nm).

%randomised GPs?
randGP = 0;     % can be set to 1 for random neighbouring.

%Range for GP capture. 
rg = 0.1;       %This sets the +/- range of GP values to include in
                %the neighbouring. E.g. a value of 0.1 would mean values
                %0.1 either side of a GP value will be included in the
                %neighbouring group.

%Tile area cut off for size filtering                
area_co = 2500; %This sets the maximum area (nm^2) to allow into the Voronoi
                %neighbouring, as larger tiles maybe from areas of low data 
                %coverage.
                
%ROI selection. 
pos = [34000,25000; 34000,30000; 39000,30000; 39000,25000];
%Coordinates of the ROI from which to take [x,y,GP] data

%% Selection of [x,y,GP] data within the pos coordinates

% selects which data are within the pos ROI
in = inpolygon(GP_data(:,1),GP_data(:,2),pos(:,1), pos(:,2));
roi = GP_data(in==1,:);

%voronoi of ROI points
voronoi(roi(:,1),roi(:,2));
[v,c]= voronoin([roi(:,1) roi(:,2)], {'Qbb'});
looptimes = [];

%Removal of the infinite tiles on the edges of the Voronoi Tesselation.
tic
vor_coords = cell(length(roi),1,1);
infdex = cell(length(roi),1,1);
for i = 1:length(c)
vor_coords{i} = v(c{i},1:2);
infdex{i} = find(vor_coords{i}==inf);
end
infs = cellfun('isempty', infdex);
cleandex = find(infs==1);
inf_tiles = find(infs==0);

clean_coords = cell(length(cleandex),1,1);
centres = zeros(length(cleandex),2);
for i = 1:length(cleandex)
clean_coords{i} = vor_coords{cleandex(i),1};
centres(i,1:2) = [roi(cleandex(i),1) roi(cleandex(i),2)];
end

%The below code can be run uncommented to show the tiles remaining (blue in
%plot) after removal of the infinite tiles.

% for j = 1:length(clean_coords)
%         patch(clean_coords{j}(:,1),clean_coords{j}(:,2),'blue')
%         %text(centres(j,1),centres(j,2),num2str(j),'Color','black');
%         %text(mean(clean_coords{j}(:,1)),mean(clean_coords{j}(:,2)),num2str(j),'Color','black');
%         hold on
% end
% xlim([min(roi(:,1)) max(roi(:,1))]); ylim([min(roi(:,2)) max(roi(:,2))]);
looptimes(1) = toc
%% size filtering by area
tic
GP_area = zeros(length(c),1);
for i = 1:length(c)
GP_area(i) = polyarea(v(c{i},1),v(c{i},2));
end

GP_area_clean = GP_area(cleandex);

szfilt_coords = find(GP_area(cleandex)<area_co);

%The below code can be run uncommented to show the tiles remaining (red in
%plot) after removal of the tiles too large for the area cut off.

% for j = 1:length(szfilt_coords)
%         patch(clean_coords{szfilt_coords(j)}(:,1),clean_coords{szfilt_coords(j)}(:,2),'red')
%         hold on
% end
% xlim([min(roi(:,1)) max(roi(:,1))]); ylim([min(roi(:,2)) max(roi(:,2))]);
looptimes(2) = toc
%% size filtering and finding GP matching tiles
%selects the tiles that remain after size filtering
tic
GP_v = clean_coords(szfilt_coords);
rands = randperm(length(GP_v))';
for i = 1:length(szfilt_coords)
    GP_v{i,3} = GP_area_clean(szfilt_coords(i),1);
    if randGP == 1 %this will randomly assign tiles regardless of GP if 
                   %randGP is equal to 1
    GP_v{i,2} = roi(szfilt_coords(rands(i)),3);
    else
    GP_v{i,2} = roi(szfilt_coords(i),3);
    end
end

for i =  1:length(GP_v)
    % this selects for each tile, all the tiles in the ROI that are within
    % the rg of the GP value for that tile.
    GP_v{i,4} = find(cell2mat(GP_v(:,2))+2 > cell2mat(GP_v(i,2))+(2-rg) & cell2mat(GP_v(:,2))+2 < cell2mat(GP_v(i,2))+(2+rg));
end
looptimes(3) = toc

%The below code can be run uncommented to display the [x,y,GP] data from
%the ROI color-coded according to the GP value.

% GP_color = [GP_v{:,2}]';
% for j = 1:length(GP_v)
%         patch(GP_v{j}(:,1),GP_v{j}(:,2),GP_color(j,1))
%         hold on
% end
% xlim([min(roi(:,1)) max(roi(:,1))]); ylim([min(roi(:,2)) max(roi(:,2))]);

%% Matching with like GP tiles first pass
% this section runs the first matching round to find all the tiles within
% the GP range that share a tile borders. This will be inserted into
% GP_v{:,5}
tic
for i = 1:length(GP_v)
    GP_neigh = cell(length(GP_v{i,4}),1,1);
    for j = 1:length(GP_v{i,4})
        if length(intersect(GP_v{GP_v{i,4}(j),1},GP_v{i,1},'rows'))>=2
        GP_neigh{j} = 1;
        end
    end
    x = cellfun('isempty', GP_neigh);
    GP_v{i,5} = GP_v{i,4}(x==0,1);
    clear GP_neigh
end
looptimes(4) = toc
%% taking neighbour indices to complete neighbouring
% this then repeats the neighbouring from the initial neighbouring result
% above until no new combinations/neighbourings of tiles can be combined.
% The final fully neighboured result will appear in compGPV.

[Au, idx ,idx2] = uniquecell(GP_v(:,5));
masterAu = Au;
clear nes Au_ints Au_ind Au_ind1

tic
while 1
    nes = cell(length(Au),1,1);
    for j = 1:length(Au)
        for i = 1:length(Au)
        nes{i} = intersect(Au{j},Au{i});
        end
    x = cellfun('isempty', nes);
    Au_ints{j} = find(x==0);
    Au_ind{j} = unique(cell2mat(Au(Au_ints{j},1)));
    Au_ind1 = uniquecell(Au_ind);
    end
    if length(Au_ind1) == length(Au)
        break
    else
        Au = Au_ind1';
        clear nes Au_ints Au_ind Au_ind1
    end
end
looptimes(5) = toc
compGPv = Au;

%% Calculating new values of combined GP tiles.
% This section uses the indexes of the neighbrous found in the previous
% section and calculates the new total area of the combined tile and the
% mean GP value for that tile. These can be used for further analysis or
% visualization (see bottom).

tic
compGPareas = zeros(length(compGPv),1);
compMeanGP = zeros(length(compGPv),1);
for i = 1:length(compGPv)
    compGPareas(i) = sum(cell2mat(GP_v(compGPv{i},3)));
    compMeanGP(i) = mean(cell2mat(GP_v(compGPv{i},2)));
end 
looptimes(6) = toc

%%Below are some codes which can be run uncommmented to plot the original or combined GP
%%tiles and color code them according to GP value

%%Run below line uncommented to generate color scale for all the tiles if
%%individual, ungrouped, tile plot desired.

% GP_color_all = cell2mat(GP_v(:,2));


%%The below loop will plot the GP tiles color coded according to their GP
%%value.

% for i = 1:length(compGPv)
% for j = cell2mat({compGPv{i}})'
%         patch(GP_v{j}(:,1),GP_v{j}(:,2),compMeanGP(i,1),'EdgeColor','none')
%         %%uncomment below line, and comment above line, to plot
%         %%individual ungrouped tiles
%         %patch(GP_v{j}(:,1),GP_v{j}(:,2),GP_color_all(j,1))
%         colormap cool
%         hold on
%         %%below can be uncommented to render in the GP value in text
%         %text(mean(GP_v{j}(:,1)),mean(GP_v{j}(:,2)),num2str(GP_color_all(j,1)),'HorizontalAlignment', 'center','Color','black');
% end
% hold on
% end
% xlim([min(roi(:,1)) max(roi(:,1))]); ylim([min(roi(:,2)) max(roi(:,2))]);