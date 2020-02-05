% Script to take a single slice from a 3d particle to 2d slices and normalize it for relion  

% _rlnImageName #14
% _rlnMagnification #9
% _rlnDetectorPixelSize #3
% _rlnGroupNumber #4
% _rlnOriginX #5
% _rlnOriginY #6


% Read motivelist
%motl = emread('/home/mweber/Miriam/Subtomogram_Averaging/FirstTestProjMatching/NAPS_Proc_SART3/MotlsAllWBP_refAF/motlAllMK_6.em'); 
%a = motl(:,motl(19,:) >= 160 & motl(19,:) <= 200); 
%b= motl(:,motl(19,:) <= 20 & motl(19, :) >= -20); 
%motl_com = [a,b]; 
motl_com = motl(:,motl(19,:) <= 160 & motl(19,:) >= 20);

% Read mask 
mask1 = artia.em.read('/path/to/mask.em');
mask2 = mask1(31,:,:);

relionDir = '/path/to/relion';
% Read all particles
% Initialize particle, planes arrays
parts = cell(size(motl_com, 2), 1);
planes = cell(size(motl_com, 2), 1);
outnames = cell(size(motl_com, 2), 1);
star = struct();
star.ImageName = cell(size(motl_com, 2), 1);
star.Magnification = int32(repmat(64000, size(motl_com, 2), 1));
star.DetectorPixelSize = repmat(28.16, size(motl_com, 2), 1);
star.GroupNumber = int32(ones(size(motl_com, 2), 1));
for i = 1:size(motl_com, 2)
    % Build name of particle
    name = sprintf('/path/to/part_%d_%d.em', motl_com(5, i), motl_com(6, i));
    % Read particle
    parts{i} = artia.em.read(name) .* -1;
    % Move particle
    parts{i} = artia.img.move(parts{i}, -motl_com(11:13, i));
    % Rotate particle
    parts{i} = artia.img.rot(parts{i}, [-motl_com(18, i) -motl_com(17, i) -motl_com(19, i)]);
    % Get plane
    planes{i} = parts{i}(31, :, :);
    % Normalize plane
    mask = logical(1 - artia.mask.sphere([1 64 64], 24, 0));
    p_mean = mean(planes{i}(mask));
    p_std = std(planes{i}(mask));
    planes{i} = (planes{i} - p_mean)./p_std;
    %multiple plane with mask
    planes{i} = planes{i} .* mask2; 
    % rotate plane only when using an y plane
    planes{i} = squeeze(planes{i});
    % uncomment for noise masking
    %noise = randn(64,64);
    %noise = (noise - mean(noise(:)))./std(noise(:));
    %planes{i} = noise .* (1-mask2) + planes{i} .* mask2;
    % Build output name
    outnames{i} = sprintf('%s/proj_31/part_%d_%d.mrc', relionDir, motl_com(5, i), motl_com(6, i));
    star.ImageName{i} = sprintf('proj_31/part_%d_%d.mrc', motl_com(5, i), motl_com(6, i));
    % Output information
    %star.GroupNumber(i) = motl(5, i);
    % Write
   artia.mrc.write(planes{i}, outnames{i}, 'float', 4.4);
end

star.OriginX = -zeros(size(motl_com, 2), 1);
star.OriginY = -zeros(size(motl_com, 2), 1);

% create a star file which can be directly used for the classification in relion  
artia.star.write(star , '/path/to/stalk_slice_31.star');



%% Add ClassNumber to motivelist 

Number = artia.star.read('/path/to/run_it050_data.star', {'ClassNumber'}, {'float'})
for i=1:size(motl_com,2)
    motl_com(20,i) = Number.ClassNumber(i); 
end


%% split motivelist

motl_class1 = motl_com(:,motl_com(20,:) == 1); 
motl_class2 = motl_com(:,motl_com(20,:) == 2); 
motl_class3 = motl_com(:,motl_com(20,:) == 3); 
motl_class4 = motl_com(:,motl_com(20,:) == 4); 
motl_class5 = motl_com(:,motl_com(20,:) == 5); 

mkdir '/path/to/motls_after_stalk_classification' 

artia.em.write(motl_class1,'/path/to/motl_class_1_proj31_1.em'); 
artia.em.write(motl_class2,'/path/to/motl_class_2_proj31_1.em');
artia.em.write(motl_class3,'/path/to/motl_class_3_proj31_1.em');
artia.em.write(motl_class4,'/path/to/motl_class_4_proj31_1.em');
artia.em.write(motl_class5,'/path/to/motl_class_5_proj31_1.em');







