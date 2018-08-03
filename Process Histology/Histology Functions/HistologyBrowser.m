function HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, ...
                use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)

% display image and set up user controls for contrast change        
ud_histology.contrast = [0 1]; 
ud_histology.break = 0; 
ud_histology.key = 0; 
ud_histology.show_original = 0; 
ud_histology.contrast_type = 1;
ud_histology.channel = 3;
ud_histology.file_num = 1;
ud_histology.num_files = length(image_file_names);
ud_histology.save_folder = save_folder;
ud_histology.image_folder = image_folder;
ud_histology.microns_per_pixel = microns_per_pixel;
ud_histology.microns_per_pixel_after_downsampling = microns_per_pixel_after_downsampling;
ud_histology.gain = gain;

% load histology image
disp(['loading image ' num2str(ud_histology.file_num) '...'])


% load already processed image
if use_already_downsampled_image
    image = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '_processed.tif']));
else %process image now
    image = imread(fullfile(image_folder,image_file_names{ud_histology.file_num}));
    original_image_size = size(image);

    % resize (downsample) image to 25 micron pixels
    image = imresize(image, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
end
original_image = image*gain;
imshow(original_image);

ud_histology.original_image = original_image;
ud_histology.adjusted_image = original_image;

if ~use_already_downsampled_image
    imwrite(ud_histology.adjusted_image, fullfile(ud_histology.save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '_processed.tif']))
end

set(histology_figure, 'UserData', ud_histology);

set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyHotkeyFcn(histology_figure, keydata, image_file_names, use_already_downsampled_image));
set(histology_figure, 'WindowScrollWheelFcn', @(src,evt)HistologyScrollFcn(histology_figure, evt))

fprintf(1, '\n Controls: \n \n');
fprintf(1, 'scroll: adjust contrast \n');
fprintf(1, 'space: switch btwn adjusting upper and lower saturation points \n');
fprintf(1, 'e: view original version \n');
fprintf(1, 'any other key: return to modified version \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'c: move to next channel \n');
fprintf(1, 's: save image \n');
fprintf(1, 'left/right arrow: save and move to next slide image \n');

fprintf(1, '\n adjusting minimum intensity limit \n');



% --------------------
% Respond to keypress
% --------------------
function HistologyHotkeyFcn(fig, keydata, image_file_names, use_already_downsampled_image)

ud = get(fig, 'UserData');

switch lower(keydata.Key)    
    case 'e' % show original
        figure(fig);
        imshow(ud.original_image)
    case 'r' % return to original
        ud.contrast = [0 1];       
    case 'space'
        contrast_effect = {'minimum intensity limit','maximum intensity limit'};
        if ud.contrast_type==2
            ud.contrast_type = 1;
        elseif ud.contrast_type==1
            ud.contrast_type = 2;
        end
        disp(['switch contrast effect -- now adjusting ' contrast_effect{ud.contrast_type}])
    case 'c' % break
        disp('next channel')
        ud.channel = ud.channel + 1 - (ud.channel==3)*3;
        ud.contrast = [0 1];
        ud.contrast_type = 1;
        disp('adjusting minimum intensity limit')        
    case 's' % save image
        disp('saving downsampled and processed image');
        imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) '_processed.tif']))
        imshow(ud.adjusted_image)
    case 'leftarrow' % save image and move to previous image
    disp('saving downsampled and processed image');
    imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) '_processed.tif']))
    imshow(ud.adjusted_image)           
        if ud.file_num > 1
            ud.file_num = ud.file_num - 1;
            move_on = true;
        else
            move_on = false;
        end
    case 'rightarrow' % save image and move to next image
    disp('saving downsampled and processed image');
    imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) '_processed.tif']))
    imshow(ud.adjusted_image)               
        if ud.file_num < ud.num_files;
            ud.file_num = ud.file_num + 1;
            move_on = true;        
        else
            disp('that''s all, folks; continue to the next cell')
            move_on = false;
        end
    otherwise
        imshow(ud.adjusted_image)
end

if (strcmp(lower(keydata.Key),'leftarrow') || strcmp(lower(keydata.Key),'rightarrow')) && move_on
            
    % load histology image
    disp(['loading image ' num2str(ud.file_num) '...'])

    % reinitialize parameters
    ud.contrast = [0 1];
    ud.contrast_type = 1;
    ud.channel = 1;
    disp('adjusting minimum intensity limit')
    
    % load already processed image
    if use_already_downsampled_image
        image = imread(fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) '_processed.tif']));
    else %process image now
        image = imread(fullfile(ud.image_folder, image_file_names{ud.file_num}));
        original_image_size = size(image);

        % resize (downsample) image to 25 micron pixels
        image = imresize(image, [round(original_image_size(1)*ud.microns_per_pixel/ud.microns_per_pixel_after_downsampling)  NaN]);
    end
    original_image = image*ud.gain;
    imshow(original_image);

    ud.original_image = original_image;
    ud.adjusted_image = original_image;
    if ~use_already_downsampled_image
        imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) '_processed.tif']))
    end
end

ud.key = keydata.Key;
set(fig, 'UserData', ud);



% --------------------
% Scrolling function
% --------------------
function HistologyScrollFcn(fig, evt)

ud = get(fig, 'UserData');
ud.key = 'scroll';


%modify based on scrolling
if ud.contrast(1) < ud.contrast(2)
    ud.contrast(ud.contrast_type) = ud.contrast(ud.contrast_type) + evt.VerticalScrollCount*.05;
else
    disp('contrast limit hit')
    ud.contrast(ud.contrast_type) = ud.contrast(3 - ud.contrast_type) - .1 * (2*(ud.contrast_type==1)-1);
end

% make sure within limit of 0 to 1
if ud.contrast(ud.contrast_type) < 0
    ud.contrast(ud.contrast_type) = 0;
    disp('contrast limit hit')
elseif ud.contrast(ud.contrast_type) > 1
    ud.contrast(ud.contrast_type) = 1;
    disp('contrast limit hit')
end

try
adjusted_image_slice = imadjust(ud.original_image(:,:,ud.channel),ud.contrast);
ud.adjusted_image(:,:,ud.channel) = adjusted_image_slice; 
imshow(ud.adjusted_image)
catch
    disp('parameter out of bounds')
end
    

set(fig, 'UserData', ud);



    