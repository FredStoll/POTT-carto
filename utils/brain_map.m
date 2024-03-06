function brain_map(gr,area2test,scales,col,fig)

if nargin<3
    scales=[];
    col = 'Blues';
    fig = figure("Color",[1 1 1]);
elseif nargin<4
    col = 'Blues';
elseif nargin<5
    fig = figure("Color",[1 1 1]);
end

if isobject(gr)

    %- normalization for line sizes
    norm = @(data,minW,maxW) (1-0.4)*(  (data-minW)/(maxW-minW)) +0.4;

    %- extract the size for the various lines
    if ~isempty(scales)
        min_w = scales(1); max_w = scales(2);
    else
        min_w = min(gr.Edges.Weight);
        max_w = max(gr.Edges.Weight);
    end
    LWidths = 0.01+(6*((gr.Edges.Weight-min_w)/(max_w-min_w)));

    %- mean weight for every node
    meanWeight = [];
    for i = 1 : length(area2test)
        meanWeight(i) = mean(gr.Edges.Weight(sum(ismember(gr.Edges.EndNodes,area2test(i))')==1));
    end

else
    meanWeight = gr(:,1);

end

%-load the brain images (only boundaries + area colored)
path4pic = 'C:\Users\Fred\Dropbox\Rudebeck Lab\Brain_pic\';
if length(area2test)<=5
    brain = imread([path4pic 'Orbital_view_brain_lines_5areas.bmp']);
    brain_col = imread([path4pic 'Orbital_view_brain_5areas.bmp']);
    %- area location and color on the brain map
    areamapping =  { 'a12mr'    'a12ol'   'a11ml'      'a13ml'      'LAI'      };
    colormapping = [0 255 0 ; 255 0 9 ;  0 200 200   ; 100 0 100 ; 0 0 150 ];
    locXY =        [1073 892 ;1304 1226;   739 742 ; 950 1199  ; 995 1528];

else
    brain = imread([path4pic 'Orbital_view_brain_lines_8areas.bmp']);
    brain_col = imread([path4pic 'Orbital_view_brain_8areas.bmp']);

    %- area location and color on the brain map
    areamapping =  {'a12r'      'a12m'    'a12o'   'a12l'     'a11ml'    'a13l'      'a13m'      'LAI'      };
    colormapping = [200 200 0 ; 0 255 0 ; 255 0 9 ; 0 0 255 ; 0 200 200 ; 0 100 0  ; 100 0 100 ; 0 0 150 ];
    locXY =        [867 630 ;  1073 946 ;1313 1226; 1281 944;   739 742 ; 1070 1293; 832 1225  ; 995 1528];

end

%- transform value to get a white background and a range that match the measure to plot
if ~isempty(scales)
    minScale = scales(1); maxScale = scales(2);
else
    minScale = min(meanWeight)-((max(meanWeight)- min(meanWeight))*0.1); %0.1
    maxScale =  max(meanWeight)+((max(meanWeight)- min(meanWeight))*0.1); %0.7
end
brainsize = size(brain_col);
backgrd=((maxScale- minScale)*0.05);
Perf_map = (minScale+backgrd)*ones(brainsize(1:2));

for i = 1 : length(area2test)
    b = ismember(areamapping,area2test{i}); % match area with brain
    ar1 = (brain_col(:,:,1)==colormapping(b,1) & brain_col(:,:,2)==colormapping(b,2) & brain_col(:,:,3)==colormapping(b,3));

  %  if ~isnan(meanWeight(i))
        Perf_map(ar1)=meanWeight(i);
  %  end
    %  Perf_map(ar1)=i; %- sanity check that color matching is correct!
end
boun = (brain(:,:,1)<100 & brain(:,:,2)<100  & brain(:,:,3)<100 );

Perf_map(boun)=minScale;

colorsJ = cbrewer('seq', col, 98);
colorsJ = [0 0 0 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; 1 1 1 ; colorsJ];
%axes1 = axes('Parent',fig);
imagesc(Perf_map(5:end-5,5:end-5),[minScale maxScale]);colormap(colorsJ); %- 5 to remove some black border
axis off;
if isobject(gr)
    % colorbar('Location','south','Position',[0.34 0.19 0.38 0.024]);
    colorbar
else
    colorbar
end
set(gca,'FontSize',16);

if isobject(gr) %- only for graph plot

    %- add lines, with size = connection strenght
    hold on
    LWidths_color = 1-repmat(norm(gr.Edges.Weight,minScale,maxScale),1,3);
    for i = 1 : length(LWidths)
        el_1 = find(ismember(areamapping,gr.Edges.EndNodes{i,1}));
        el_2 = find(ismember(areamapping,gr.Edges.EndNodes{i,2}));
        vect = [locXY(el_1,:) ; locXY(el_2,:)  ];

        line(vect(:,1),vect(:,2),'LineWidth',LWidths(i),'Color',LWidths_color(i,:));hold on
    end
    plot(locXY(:,1),locXY(:,2),'.','MarkerSize',40,'Color',[253 141 60]/255) %- [66 146 200]/255
    %  for i = 1 : length(locXY)
    %     text(locXY(i,1),locXY(i,2)+50,area2plot{i},'Color',[66 146 200]/255)
    %  end

end