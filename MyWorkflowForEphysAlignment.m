% WORKFLOW TO ALIGN EPHYS ANATOMY PROBE TRACK WITH ALLEN BRAIN ATLAS

% 1) obtain X microns of one image and X pixels, important for the
%    downsampling for the alignment
% 2) export images in tiff from FIJI in order to have both Probe
%    track and DAPI staining in one single-channel RGB image. (rember to
%    flip all the images to have them on the left emisphere
% 3) crop every single anatomy section and make one folder with all images
%    for every Mice (folders separated)

% - - - - - - - - - - - - %
% Set now some parameters %

% write here the folder were are the images (they should be tiff)
image_folder =  'E:\Users\filippo.michelon\Documents\EXPERIMENTS - RUNNING ANALYSIS\ANATOMY_NPX\2Conc5Od_Rec_NPX_Anatomy_0522\M3';
% write here the folder were processed images will be saved (a new folder named processed will be created)
save_folder = image_folder;
% directory of histology
processed_images_folder = append(image_folder,'\processed'); 
% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = 'M3'; 
% root folder
root = 'E:\Users\filippo.michelon\Documents\MATLAB\Allen CCF';
% the real probe length
RealProbeLengths = RealProbeLengths(3); % choose here the mouse from this variable that must be manually created obtaining real length from REC!

% set image info for downsampling
xMicronsImage = 8761.10;
xPixelsImage = 4563;
% decide the basefile name for output slices
save_file_name = probe_save_name_suffix;
plane = 'coronal';
% increase gain if for some reason the images are not bright enough
gain = 1; 

% directory of reference atlas files
annotation_volume_location = 'E:\Users\filippo.michelon\Documents\MATLAB\Allen CCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'E:\Users\filippo.michelon\Documents\MATLAB\Allen CCF\structure_tree_safe.csv';
template_volume_location = 'E:\Users\filippo.michelon\Documents\MATLAB\Allen CCF\template_volume_10um.npy';





%%  SET FILE AND PARAMETERS

% * remember to run one cell at a time, instead of the whole script at once *

% directory of histology images
image_folder = image_folder;

% directory to save the processed images -- can be the same as the above image_folder
% results will be put inside a new folder called 'processed' inside of this image_folder
save_folder = save_folder; 

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep '*.tif']); % get the contents of the image_folder % original was tif
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order

% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas

% note: since x and y are equally sampled in the image it doesn't matter if
% you use microns and pixels of x or y:
xMicronsImage = xMicronsImage;
xPixelsImage = xPixelsImage;

microns_per_pixel = xMicronsImage / xPixelsImage;     
microns_per_pixel_after_downsampling = 10;


% ----------------------
% additional parameters
% ----------------------

% if the images are cropped (image_file_are_individual_slices = false),
% name to save cropped slices as; e.g. the third cropped slice from the 2nd
% image containing many slices will be saved as: save_folder/processed/save_file_name02_003.tif
save_file_name = save_file_name;

% increase gain if for some reason the images are not bright enough
gain = gain; 

% plane to view ('coronal', 'sagittal', 'transverse')
plane = plane;

% size in pixels of reference atlas brain. For coronal slice, this is 800 x 1140
if strcmp(plane,'coronal')
    atlas_reference_size = [800 1140]; 
elseif strcmp(plane,'sagittal')
    atlas_reference_size = [800 1320]; 
elseif strcmp(plane,'transverse')
    atlas_reference_size = [1140 1320];
end


% finds or creates a folder location for processed images -- 
% a folder within save_folder called processed
folder_processed_images = fullfile(save_folder, 'processed');
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end

%% LOAD AND PROCESS SLICE PLATE IMAGES

% close all figures
close all
   
% if the images need to be downsampled to 10um pixels (use_already_downsampled_image = false), 
% this will downsample and allow you to adjust contrast of each channel of each image from image_file_names
%
% if the images are already downsampled (use_already_downsampled_image = true), this will allow
% you to adjust the contrast of each channel
%
% Open Histology Viewer figure
try; figure(histology_figure);
catch; histology_figure = figure('Name','Histology Viewer'); end
warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');

% Function to downsample and adjust histology image
HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, folder_processed_images, image_files_are_individual_slices, ...
            use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)

  
%% GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, SHARPEN, and CHANGE ORDER

% close all figures
close all
            
% this takes images from folder_processed_images ([save_folder/processed]),
% and allows you to rotate, flip, sharpen, crop, and switch their order, so they
% are in anterior->posterior or posterior->anterior order, and aesthetically pleasing
% 
% it also pads images smaller than the reference_size and requests that you
% crop images larger than this size
%
% note -- presssing left or right arrow saves the modified image, so be
% sure to do this even after modifying the last slice in the folder
slice_figure = figure('Name','Slice Viewer');
SliceFlipper(slice_figure, folder_processed_images, atlas_reference_size)

%% (THE ABOVE PART WAS TAKEN FROM Process_Histology.m)
%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
    tv_plot = tv;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
    tv_plot = permute(tv,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
    tv_plot = permute(tv,[2 3 1]);
end

% create Atlas viewer figure
f = figure('Name','Atlas Viewer'); 

% % show histology in Slice Viewer
try; figure(slice_figure_browser); title('');
catch; slice_figure_browser = figure('Name','Slice Viewer'); end
 reference_size = size(tv_plot);
sliceBrowser(slice_figure_browser, processed_images_folder, f, reference_size);

set(gcf,'Position',[10 300 950 700]);

% % use application in Atlas Transform Viewer
% % use this function if you have a processed_images_folder with appropriately processed .tif histology images
f = AtlasTransformBrowser(f, tv_plot, av_plot, st, slice_figure_browser, processed_images_folder, probe_save_name_suffix, plane);

set(gcf,'Position',[950 300 950 700]);

%% (THE ABOVE PART WAS TAKEN FROM Navigate_Atlas_and_Register_Slices.m)
%% %% how to calcuate the actual lenght of the probe in the brain?? > take it during the recording!
% %  Analyze_Clicked_Points
% % 
% % z = max(roi_location(:,2)) - min(roi_location(:,2));
% % x = max(roi_location(:,3)) - min(roi_location(:,3));
% % y = max(roi_location(:,1)) - min(roi_location(:,1));
% % xz = hypot(x,z);
% % xyz = hypot(xz,y); % approximation of probe length

% % cd(processed_images_folder)
% % filename = sprintf(append('probe_lenght', probe_save_name_suffix,'.mat'));
% % save(filename,'xyz')
% % cd(root)

%% THIS PART IS TO VISUALIZE REGIONS TRAVERSED BY THE ELECTRODE

% either set to 'all' or a list of indices from the clicked probes in this file, e.g. [2,3]
probes_to_analyze = 'all';  % [1 2]

% --------------
% key parameters
% --------------
% how far into the brain did you go from the surface, either for each probe or just one number for all -- in mm
probe_lengths = RealProbeLengths; %xyz; 
% % % % % % % NB!!!remeber to remove 100um!!!!!!!!!!!!!!WORK HERE!!!!

% from the bottom tip, how much of the probe contained recording sites -- in mm
active_probe_length = 3.84;

% distance queried for confidence metric -- in um
probe_radius = 100; 

% overlay the distance between parent regions in gray (this takes a while)
show_parent_category = false; 

% plot this far or to the bottom of the brain, whichever is shorter -- in mm
distance_past_tip_to_plot = 0.5;

% set scaling e.g. based on lining up the ephys with the atlas
% set to *false* to get scaling automatically from the clicked points
scaling_factor = false;


% ---------------------
% additional parameters
% ---------------------
% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% probe insertion direction 'down' (i.e. from the dorsal surface, downward -- most common!) 
% or 'up' (from a ventral surface, upward)
probe_insertion_direction = 'down';

% show a table of regions that the probe goes through, in the console
show_region_table = true;

% black brain?
black_brain = true;

% GET AND PLOT PROBE VECTOR IN ATLAS SPACE

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
end

% load probe points
probePoints = load(fullfile(processed_images_folder, ['probe_points' probe_save_name_suffix]));
ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};
fwireframe = [];

% scale active_probe_length appropriately
active_probe_length = active_probe_length*100;

% determine which probes to analyze
if strcmp(probes_to_analyze,'all')
    probes = 1:size(probePoints.pointList.pointList,1);
else
    probes = probes_to_analyze;
end 





% PLOT EACH PROBE -- FIRST FIND ITS TRAJECTORY IN REFERENCE SPACE

% create a new figure with wireframe
fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
hold on; 
fwireframe.InvertHardcopy = 'off';

for selected_probe = probes
    
% get the probe points for the currently analyzed probe 
if strcmp(plane,'coronal')
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [3 2 1]);
elseif strcmp(plane,'sagittal')
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 2 3]);
elseif strcmp(plane,'transverse')
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 3 2]);
end



% get user-defined probe length from experiment
if length(probe_lengths) > 1
    probe_length = probe_lengths(selected_probe);
else
    probe_length = probe_lengths;
end

% get the scaling-factor method to use
if scaling_factor
    use_tip_to_get_reference_probe_length = false;
    reference_probe_length = probe_length * scaling_factor;
    disp(['probe scaling of ' num2str(scaling_factor) ' determined by user input']);    
else
    use_tip_to_get_reference_probe_length = true;
    disp(['getting probe scaling from histology data...']);
end

% get line of best fit through points
% m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
[m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
if isnan(m(1))
    disp(['no points found for probe ' num2str(selected_probe)])
    continue
end

% ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
if p(2)<0
    p = -p;
end

% determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
ann = 10;
out_of_brain = false;
while ~(ann==1 && out_of_brain) % && distance_stepped > .5*active_probe_length)
    m = m-p; % step 10um, backwards up the track
    ann = av_plot(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
    if strcmp(st.safe_name(ann), 'root')
        % make sure this isn't just a 'root' area within the brain
        m_further_up = m - p*20; % is there more brain 200 microns up along the track?
        ann_further_up = av_plot(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
        if strcmp(st.safe_name(ann_further_up), 'root')
            out_of_brain = true;
        end
    end
end

% focus on wireframe plot
figure(fwireframe);

% plot probe points
hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);

% plot brain entry point
plot3(m(1), m(3), m(2), 'r*','linewidth',1)

% use the deepest clicked point as the tip of the probe, if no scaling provided (scaling_factor = false)
if use_tip_to_get_reference_probe_length
    % find length of probe in reference atlas space
    if strcmp(probe_insertion_direction, 'down')
        [depth, tip_index] = max(curr_probePoints(:,2));
    elseif strcmp(probe_insertion_direction, 'up')
        [depth, tip_index] = min(curr_probePoints(:,2));    
    end
    reference_probe_length_tip = sqrt(sum((curr_probePoints(tip_index,:) - m).^2)); 
    
    % and the corresponding scaling factor
    shrinkage_factor = (reference_probe_length_tip / 100) / probe_length;
    
    % display the scaling
    disp(['probe length of ' num2str(reference_probe_length_tip/100) ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
    disp(['probe scaling of ' num2str(shrinkage_factor)]); disp(' ');
    
    % plot line the length of the probe in reference space
    probe_length_histo = round(reference_probe_length_tip);
    
% if scaling_factor is user-defined as some number, use it to plot the length of the probe
else 
    probe_length_histo = round(reference_probe_length * 100); 
end

% find the percent of the probe occupied by electrodes
percent_of_tract_with_active_sites = min([active_probe_length / (probe_length*100), 1.0]);
active_site_start = probe_length_histo*(1-percent_of_tract_with_active_sites);
active_probe_position = round([active_site_start  probe_length_histo]);

% plot line the length of the active probe sites in reference space
plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], m(3)+p(3)*[active_probe_position(1) active_probe_position(2)], m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
    'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);
% plot line the length of the entire probe in reference space
plot3(m(1)+p(1)*[1 probe_length_histo], m(3)+p(3)*[1 probe_length_histo], m(2)+p(2)*[1 probe_length_histo], ...
    'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');


%% ----------------------------------------------------------------
% Get and plot brain region labels along the extent of each probe
% ----------------------------------------------------------------

% convert error radius into mm
error_length = round(probe_radius / 10);

% find and regions the probe goes through, confidence in those regions, and plot them
[borders,borders_table,fD] = plotDistToNearestToTip(m, p, av_plot, st, probe_length_histo, error_length, active_site_start, distance_past_tip_to_plot, show_parent_category, show_region_table, plane); % plots confidence score based on distance to nearest region along probe
title(['Probe ' num2str(selected_probe)],'color',ProbeColors(selected_probe,:))

pause(.05)
end

TheName = append('Regions_Electrode',probe_save_name_suffix);

RegionsElectrode.Mouse = TheName;
RegionsElectrode.borders = borders;
RegionsElectrode.borders_table = borders_table;


% (THE ABOVE PART WAS TAKEN FROM Display_Probe_Track.m)
% % % % % % % % % % % Work here to map channels to anatomical positions

TopChannel = active_probe_position(1);
BottomChannel = active_probe_position(2);
depthsAnatomy = linspace(TopChannel,BottomChannel,192); % 192 is the number of doublets of channels in NPX
depthsAnatomy = [depthsAnatomy; depthsAnatomy];
depthsAnatomy = depthsAnatomy(:);
depthsAnatomy = round(depthsAnatomy*10);
% put location for every channel in a file
for channel = 1 : size(depthsAnatomy,1)
   AnatPosition = size( find(depthsAnatomy(channel) > borders_table.upperBorder) | find(depthsAnatomy(channel) > borders_table.upperBorder),1);
        
   if AnatPosition == 0
       AnatPosition = 1;
   end
   
     ChPos(channel).ChanNum = channel;
     ChPos(channel).ChannelDepth = depthsAnatomy(channel);
     ChPos(channel).Acronym = borders_table.acronym(AnatPosition);
     ChPos(channel).parentID = borders_table.parentID(AnatPosition);
     ChPos(channel).upperBorder = borders_table.upperBorder(AnatPosition);
     ChPos(channel).lowerBorder = borders_table.lowerBorder(AnatPosition);
     ChPos(channel).name = borders_table.name(AnatPosition);
     ChPos(channel).avIndex = borders_table.avIndex(AnatPosition);
     
end

% look for channels that are near parent borders changes, to remove them in further analysis
for channel = 1 : size(depthsAnatomy,1)-2
    if  ChPos(channel).parentID == ChPos(channel+1).parentID & ChPos(channel).parentID == ChPos(channel+2).parentID
       Down(channel) = 0;
    else
        Down(channel) = 1;
    end
end
for channel = 3 : size(depthsAnatomy,1)
    if  ChPos(channel).parentID == ChPos(channel-1).parentID & ChPos(channel).parentID == ChPos(channel-2).parentID
       Up(channel) = 0;
    else
        Up(channel) = 1;
    end
end

Toremove = Up+[Down 0 0];
for channel = 1 : size(depthsAnatomy,1)
ChPos(channel).UncertainPos = Toremove(channel);
end

cd(processed_images_folder)
filename = sprintf(append('ChPos', probe_save_name_suffix,'.mat'));
save(filename,'ChPos')
filename = sprintf(append('Regions_Electrode', probe_save_name_suffix,'.mat'));
save(filename,'RegionsElectrode')
print('fD','-dpng')

cd(root)


%%
% this is a usefull part to plot the channel points in the 3d map of brain, will be potentially usefull
% % % start_of_active_sites = [m(1)+p(1)*active_probe_position(1) m(2)+p(2)*active_probe_position(1) m(3)+p(3)*active_probe_position(1)]; % location in atlas were the active sites start
% % % probe_tip = [m(1)+p(1)*active_probe_position(2) m(2)+p(2)*active_probe_position(2) m(3)+p(3)*active_probe_position(2)]; % tip of the probe in atlas
% % % 
% % % % since there are 384 channels in two columns I will put 384/2 channels evenly spaced from 0 to 1 (0 tip probe / 1 upper part probe)
% % % channelLocation = 0:1/191:1; 
% % % channelLocation = [channelLocation; channelLocation];
% % % channelLocation = channelLocation(:)';
% % % counter = 1;
% % % for channel = channelLocation
% % % Location_Sites(counter,:) = start_of_active_sites * channel + probe_tip * (1 - channel);
% % % counter = counter + 1;
% % % end
% % % % look at top middle and end points
% % % loc = 0.5;
% % % Location_Sites = start_of_active_sites * loc + probe_tip * (1 - loc);
% % % scatter3(Location_Sites(1), Location_Sites(3),Location_Sites(2),'MarkerFaceColor',[0 1 1]);hold on
% % % scatter3(start_of_active_sites(1),start_of_active_sites(3),start_of_active_sites(2),'MarkerFaceColor',[1 0 0])
% % % scatter3(probe_tip(1),probe_tip(3),probe_tip(2),'MarkerFaceColor',[0 1 0])
% % % % look at all the sites
% % % for i = 1:384
% % % scatter3(Location_Sites(i,1), Location_Sites(i,3),Location_Sites(i,2),50,'MarkerFaceColor',[1 .2 .2],'MarkerFacealpha',0,...
% % %     'MarkerEdgeColor',[1 .2 .2],'MarkerEdgealpha',1);hold on
% % % end
% % % load here general coordinates of the NPX probe sites
% % load('E:\Users\filippo.michelon\Documents\MATLAB\Allen CCF\CoordNPX.mat')
% % % on the x axis probe sites has spacing 32um
% % % on the y axis probe sites has spacing 20um
% % coordNPX = [coordNPX(:,1)-(11+16), coordNPX(:,2)]
% % % given the fact that in the x dimension I don't have the possibility to
% % % know precisely is which direction put the offset I will leave them equal,
% % % in any case I don't have the spacial resolution to disambiguate between
% % % two regions that are 32um apart.