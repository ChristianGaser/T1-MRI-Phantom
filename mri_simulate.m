function mri_simulate(simu, rf)
% MRI_SIMULATE - Simulates MR images with various optional features
%
% Overview:
%   `mri_simulate` generates simulated MRI images, allowing for T1-weighted
%   (T1w) imaging simulations. It employs a high-quality (high-resolution, low
%   noise) T1w image as a base. Users can introduce various artifacts and
%   features like white matter hyperintensities (WMHs), noise, and RF B1
%   inhomogeneities.
%   Noise can be added as either:
%     • Gaussian noise specified as a percentage of the WM mean (simu.pn)
%     • Rician magnitude noise at a target WM SNR (simu.snrWM>0)
%   It supports simulations of atrophy or cortical thickness modifications. 
%   Preprocessing with SPM12 segmentation is required for custom images.
%   Requirements: SPM12 or SPM25 with CAT12 >= 12.10 installed.
%
%   The function also writes JSON sidecars next to the main image
%   containing key simulation metadata (tool/version, voxel size, noise/SNR,
%   RF field parameters, thickness tags), using SPM's spm_jsonwrite.
%
%   Thickness/PVE pipeline (when simu.thickness is set):
%   - A constant cortical thickness (global or region-wise) is synthesized by
%     expanding GM outward from the original WM using a Euclidean distance map.
%   - To obtain partial-volume-like transitions, the label boundary is shifted
%     across 20 sub-voxel offsets in the range [-0.25, 0.25] voxels, simulating
%     realistic boundary uncertainty. Each offset yields a hard label image
%     (CSF=1, GM=2, WM=3), and the results are averaged to produce a smooth
%     PVE-like label map.
%   - The final PVE label map is then used to synthesize a T1 image from the
%     SPM segmentation model by replacing GM/WM/CSF class posteriors and using
%     their Gaussian mixture parameters (means, variances, weights).
%
% Syntax:
%   mri_simulate(simu, rf)
%
% Parameters:
%   simu (struct): Simulation parameters. Defaults are applied for missing fields.
%       - 'name' (char): T1-weighted input image filename. Default: '' (empty),
%         which triggers an interactive file selection dialog (T1 only).
%       - 'pn' (double): Percentage noise level (Gaussian) relative to the WM
%         mean intensity. Ignored if 'snrWM' > 0. Default: 0 (percent of WM).
%       - 'snrWM' (double): If >0, adds Rician noise at a user-defined SNR for
%         white matter. Uses the (noise-free) WM mean to compute the complex
%         noise sigma via sigma = WMmean / snrWM, and generates magnitude
%         Rician noise: sqrt((S + n1).^2 + n2.^2). Default: 30.
%       - 'rng' (double, NaN or []): Seed for the random number generator. 
%         Default: NaN (reproducible noise across runs). Set [] to use MATLAB's 
%         default RNG behavior (non-deterministic across sessions).
%       - 'derivative' (logical): If true, save outputs under a BIDS-style
%         derivatives folder at the dataset root, using the pipeline name
%         'mri_simulate-1.0' and mirroring the subject/session path under it.
%         Example: root/derivatives/mri_simulate-1.0/sub-*/ses-*/[anat|...]/
%         Default: 0 (save next to input file).
%       - 'contrast' (double): Power-law contrast-change exponent applied to the
%         simulated image intensities after noise. The image is normalized to
%         [0,1], transformed as Y.^contrast, and rescaled back to its original
%         min/max range. Use values >1 to increase contrast, <1 to reduce.
%         Default: 1 (no change). Meaningful values to simulate contrast
%         are 0.5 (low contrast) and 1.5 (high contrast).
%       - 'resolution' (double or [x, y, z]): Spatial resolution of the
%         simulated image. Default: NaN (keep original resolution). If scalar,
%         it is applied to all three axes; if a 3-vector, each axis is set individually.
%       - 'closeWMHholes' (logical): Detect and fill WMHs in WM and correct
%         segmenations to obtain a clean simulated image without WMHs (which 
%         later allows to simulate additional WMHs using the WMH option).
%       - 'WMH' (integer or >=1 scalar): Strength of simulated white matter
%         hyperintensities (WMHs).
%           0  -> no WMHs
%           1  -> mild, mainly periventricular patches
%           2  -> medium extent/contrast
%           3  -> strong, widespread WMHs
%         Notes:
%           • Values >=1 are allowed (not just 1/2/3). Higher values broaden
%             the WMH extent by reshaping the WMH prior with an exponent
%             1/(WMH-0.8), and adjust the contribution to the label map by a
%             scaling of ~1/WMH^0.75 to keep intensities in a plausible range.
%           • WMHs are constrained to (eroded) WM and modulated by a random
%             field to produce patchy distributions.
%         Default: 0 (disabled).
%       - 'atrophy' (cell): Specifies regions of interest (ROIs) for simulating
%         atrophy, including atlas name, ROI IDs, and atrophy values. Multiple ROIs 
%         and the respective atrophy values can be defined. An atrophy value of 
%         1.5 leads to a GM reduction of about 5%, while a value of 2 corresponds 
%         to around 10% and 3 to around 15%. Please note, that this option takes 
%         a lot of time because the atlas labels have to be interpolated using 
%         categorical interpolation (i.e. each label separetely). Either thickness 
%         or atrophy can be simulated. Default: [] (disabled).
%       - 'thickness' (double or [double double double]): Specifies the cortical 
%         thickness for simulation. The WM label of the image is used to add a 
%         layer of GM with a defined cortical thickness to have constant thickness. 
%         Unlike all other simulations, we can only use the modified label image 
%         for the MRI simulation and not the tissue probabilities (which give a 
%         more detailed and realistic T1w image).
%         Either a scalar value for global constant thickness or a vector with 
%         3 thickness values can be defined for the occipital and frontal lobes 
%         (1st and 3rd values) and the rest of the brain (2nd value). The 
%         Hammer atlas is used to define these areas, as well as subcortical areas 
%         and the cerebellum, which are excluded from the thickness simulation to 
%         obtain a more realistic MRI. Either thickness or atrophy can be simulated.
%         Default: 0 (disabled).
%   rf (struct): RF bias field parameters.
%       - 'percent' (double): Amplitude of the bias field in percentage.
%         Negative values invert the field. Default: 30.
%       - 'type' (char or [int, int]): Specifies the bias field type:
%           'A' | 'B' | 'C' for predefined MNI fields, or [strength, rngSeed]
%           for a simulated field. Meaningful strength values are 1..4; 3–4
%           resemble stronger inhomogeneities (e.g., 7T without correction).
%         Default: [2 0] (simulated field with moderate strength and RNG seed 0).
%       - 'save' (logical): Save the simulated bias field only when 'type' is
%         numeric (simulated). Ignored for predefined 'A'/'B'/'C'. Default: 0.
%
% Optional Inputs and Defaults:
%   - If 'simu' and/or 'rf' are omitted, default values are applied internally.
%     You can also pass a partial struct; any missing fields are filled with
%     defaults (i.e., user-specified fields override defaults only where set).
%   - Interactive mode: if 'simu.name' is empty (''), the function opens a file
%     selection dialog (spm_select) and runs the simulation for the chosen file(s).
%     This allows you to launch the tool without pre-specifying the input path.
%   See the examples below for constructing minimal 'simu' and 'rf' structs.
%
% Outputs:
%   Simulated MRI image files based on the specified parameters and features:
%     - Main simulated image (full brain)
%     - Masked simulated image (brain only)
%     - Ground-truth PVE label image
%     - Optional: RF bias field (for simulated fields)
%     - JSON sidecars for main images with simulation metadata
%
% Usage:
%   To simulate an MRI, specify the simulation (`simu`) and RF bias field
%   (`rf`) parameters:
%       mri_simulate(simu, rf);
%
% Examples:
%   Example 1 - Basic simulation with specific SNR with 0.5mm voxel size:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 20,...
%                     'resolution', 0.5);
%       rf = struct('percent', 20, 'type', 'A','save',0);
%       mri_simulate(simu, rf);
%
%   Example 2 - Advanced simulation with atrophy (10% in left middle frontal gyrus 
%               and 15% in right middle frontal gyrus based on Hammers atlas), 
%               custom RF field and thicker slices: 
%       simu = struct('name', 'custom_t1.nii', 'snrWM', 20,...
%                     'resolution', [0.5, 0.5, 1.5]);
%       simu.atrophy = {'hammers',[28, 29], [2, 3]};
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
%   Example 3 - Thickness simulation with 3 different thickness values using 
%   original voxel size:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 30,...
%                     'resolution', NaN,...
%                     'thickness', [1.5 2.0 2.5]);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
%   Example 4 - Simulation with custom RF field and added WMHs (medium strength)
%       simu = struct('name', 'custom_t1.nii', 'snrWM', 30,...
%                     'resolution', NaN, 'WMH', 2);
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
%   Example 5 - Apply contrast change (power-law, high contrast)
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 30, ...
%                     'resolution', NaN, 'contrast', 2);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
%
%
% TODO: simulation of motion artefacts using FFT and shift of phase information

version = '0.9.7';

if ~exist('cat_main_LASsimple')
  error('Please update to a newer version >=CAT12.10 to use mri_simulate')
end

% Default simulation parameters
def.name       = '';
def.pn         = 0;
def.resolution = NaN;
def.WMH        = 0;
def.atrophy    = [];
def.thickness  = 0;
def.rng        = 0;
def.snrWM      = 30;
def.contrast   = 1;  % power-law contrast change exponent (1 = unchanged)
def.derivative = 1;  % save outputs into BIDS derivatives
def.closeWMHholes = 1; % close WMHs inside deep WM

if nargin < 1, simu = def;
else, simu = cat_io_checkinopt(simu, def); end

% keep requested output resolution separate from internal thickness resampling
requested_resolution = simu.resolution;

% Default bias field parameters
def.percent    = 30;
def.type       = [2 0];
def.save       = 0;

if nargin < 2, rf = def;
else, rf = cat_io_checkinopt(rf, def); end

% call tool interactively
if isempty(simu.name)
  simu.name = spm_select(Inf,'image','Select T1w image(s) for simulation');
end

n_images = size(simu.name,1);
if n_images > 1
  P = simu.name;
  parfor i=1:n_images
    simu_i = simu;
    simu_i.name = deblank(P(i,:));
    mri_simulate(simu_i, rf);
  end
  return
end

[pth, name, ext, ~] = spm_fileparts(simu.name);
if ~exist(fullfile(pth, [name ext]), 'file')
  fprintf("File %s not found.\n", simu.name);
  return
end
if strcmp(ext,'.gz')
  fname = gunzip(fullfile(pth, [name ext]));
  simu.name = fname{1};
  [pth, name, ext] = spm_fileparts(simu.name);
  is_gz = 1;
else
  is_gz = 0;
end

% keep output names clean if a temporary resample suffix is present
name_out = regexprep(name,'_thicknessRes05$','');

% if simu.rng is not defined we use the filename to create a seed
if isempty(simu.rng) | isnan(simu.rng)
  simu.rng = sum(double(name));
end

pth_root = fileparts(which(mfilename));

% Determine output folder (optionally BIDS derivatives)
out_pth = pth;
if isfield(simu,'derivative') && simu.derivative
  try
    parts = strsplit(pth, filesep);
    if isempty(parts{1}), parts{1} = filesep; end
    idx_sub = find(strncmp(parts,'sub-',4), 1, 'first');
    if ~isempty(idx_sub) && idx_sub > 1
      % Derivatives at dataset root
      root_dir = fullfile(parts{1:idx_sub-1});
      rel_parts = parts(idx_sub:end); % sub-..[/ses-..]/anat/...
      if simu.thickness(1) > 0
        pipeline_dir = fullfile(root_dir, 'derivatives', ['mri_simulate_thickness-' version]);
      else
        pipeline_dir = fullfile(root_dir, 'derivatives', ['mri_simulate-' version]);
      end
      out_pth = fullfile(pipeline_dir, rel_parts{:});
    else
      % Fallback: place outputs under derivatives without mirroring
      root_dir = fileparts(pth);
      if simu.thickness(1) > 0
        pipeline_dir = fullfile(root_dir, 'derivatives', ['mri_simulate_thickness-' version]);
      else
        pipeline_dir = fullfile(root_dir, 'derivatives', ['mri_simulate-' version]);
      end
      out_pth = pipeline_dir;
    end
    if ~exist(out_pth,'dir'), mkdir(out_pth); end
  catch
    % In case of any issue, fall back to input folder
    out_pth = pth;
  end
end

% name of seg8.mat file that contains SPM12 segmentation parameters
mat_name = fullfile(pth, [name '_seg8.mat']);

% CAT12 template dir is later used for defining atrophy atlas
template_dir = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym');

% call SPM segmentation if necessary and only save the seg8.mat file
if ~exist(mat_name,'file')
  matlabbatch{1}.spm.spatial.preproc.channel.vols = {simu.name};
  spm_jobman('run',matlabbatch);
  clear matlabbatch

  % remove native segmentations that were saved
  for i=1:5
    spm_unlink(fullfile(pth,['c' num2str(i) name ext]));
  end
end

% check that only one of these options is used
if ((isfield(simu,'atrophy') && ~isempty(simu.atrophy) && any(simu.atrophy{3} > 1))) && any(simu.thickness)
  fprintf('Option to simulate atrophy cannot be combined with option to simulate thickness images.\n');
  return
end

% check that strength for simulated bias field is at least 1
if isnumeric(rf.type) && rf.type(1) < 1
  fprintf('Strength of simulated bias field should be at least 1.\n');
  return
end

% check that rf.save is not set for predefined MNI bias fields
if ischar(rf.type) && rf.save
  fprintf('Predifined MNI bias fields cannot be saved.\n');
  rf.save = 0;
end

% any parameter to simulate atrophy defined?
if isfield(simu,'atrophy') && ~isempty(simu.atrophy) && any(simu.atrophy{3} > 1)
  simu_atrophy = 1;
  atlas_name = fullfile(template_dir,[simu.atrophy{1} '.nii']);
  
  % check that atlas exists
  if ~exist(atlas_name,'file')
    fprintf('Atlas %s could not be found. Please check that the name is correct.\n', atlas_name);
    return
  end
else
  simu_atrophy = 0;
end

if simu.WMH < 1 && simu.WMH ~= 0
  fprintf('Strength of simulating WMHs should be either 0 or >=1.\n');
  return
end

% load seg8.mat file and define some parameters
res = load(mat_name);

if size(res.mn,1) > 1
  fprintf('Multi-modal segmentation is not recommended! Try again to segment the image using T1w data only.\n');
end

% get means for GM/WM/CSF
mn = zeros(3,1);
for k=1:3
  ind = find(res.lkp==k);
  for j=1:numel(ind)
    mn(k) = mn(k) + res.mg(ind(j))*res.mn(1,ind(j));
  end  
end

% check that it's indeed T1w data by checking CSF < GM < WM
[~, ind] = sort(mn);
if ~isequal(ind(:).',[3 1 2])
  fprintf('Warning: No typical T1w intensities were found. Please note that segmentation quality can be much lower for non T1w data.\n');
  if simu.WMH
    simu.WMH = 0;
    fprintf('Warning: WMHs can be only simulated for T1w data. WMH option was therefore disabled.\n');
  end
end

[~, bname, ext] = spm_fileparts(res.image(1).fname);
res.image(1) = spm_vol(fullfile(pth,[bname ext]));
V = res.image(1);
dim   = V.dim(1:3);
vx = sqrt(sum(V.mat(1:3,1:3).^2));

% keep native geometry for final output if needed
V_native = V;
dim_native = dim;
vx_native = vx;

% inverse deformation field name must follow the original image filename
[idef_pth, idef_bname, idef_ext] = spm_fileparts(res.image(1).fname);
idef_name_orig = fullfile(idef_pth, ['iy_' idef_bname idef_ext]);
idef_name = idef_name_orig;
thickness_resampled = false;
resampled_name = '';

% obtain SPM segmentations and write inverse deformation field
[Ysrc, Ycls, Yy] = cat_spm_preproc_write8(res,zeros(max(res.lkp),4),zeros(2,2),[1 0],0,2);

% get tissue thresholds for CSF/GM/WM (see get_tissue_thresholds for details)
T3th = get_tissue_thresholds(Ysrc, Ycls, res);

% LAS correction and SANLM denoising
Ycorr = cat_main_LASsimple(Ysrc, Ycls, T3th);
cat_sanlm(Ycorr,3,1);

% Replace GM/WM/CSF segmentation by labels using LAS corrected image
Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
Yseg = zeros([dim, 3]);

% Use CAT12 adaptive probability region-growing (APRG) approach for
% skull-stripping (uses T3th anchors; see skull_strip_APRG)
brainmask = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th);
Ycorr = Ycorr.*brainmask;

seg_order = [2 3 1];
for i = 1:3
    Yseg(:,:,:,i) = Yp0toC(3*Ycorr, seg_order(i)); 
end

% optionally close WMHs within deep WM using CAT12 approach
if isfield(simu,'closeWMHholes') && simu.closeWMHholes
  Yseg = close_WM_GM_holes(Yseg, Ysrc, Ycorr, Ycls, Yy, res, vx);
end
clear Ycls Ycorr

% For thickness simulation, resample to 0.5mm before atrophy/thickness
if any(simu.thickness)
  target_res = [0.5 0.5 0.5];
  if any(abs(vx - target_res) > 1e-6)
    V0 = V;
    Vres_tmp = V;
    Vres_tmp.dim = round(V.dim.*vx./target_res);
    Ptmp = spm_imatrix(V.mat);
    Ptmp(7:9) = Ptmp(7:9)./vx.*target_res;
    Ptmp(1:3) = Ptmp(1:3) + vx - target_res;
    Vres_tmp.mat = spm_matrix(Ptmp);

    % resample source image
    volres_tmp = zeros(Vres_tmp.dim, 'single');
    for sl = 1:Vres_tmp.dim(3)
      M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
      M1 = Vres_tmp.mat\V0.mat\M;
      volres_tmp(:,:,sl) = spm_slice_vol(V0, M1, Vres_tmp.dim(1:2), -5);
    end

    resampled_name = fullfile(pth, [name '_thicknessRes05.nii']);
    Vres_tmp.fname = resampled_name;
    Vres_tmp.dt    = [16 0];
    Vres_tmp.pinfo = [1 0 352]';
    spm_write_vol(Vres_tmp, volres_tmp);

    res.image(1) = spm_vol(resampled_name);
    V = res.image(1);
    dim   = V.dim(1:3);
    vx = sqrt(sum(V.mat(1:3,1:3).^2));

    % resample Yseg
    Yseg_res = zeros([Vres_tmp.dim 3], 'single');
    for sl = 1:Vres_tmp.dim(3)
      M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
      M1 = Vres_tmp.mat\V0.mat\M;
      for j = 1:3
        Yseg_res(:,:,sl,j) = spm_slice_vol(Yseg(:,:,:,j), M1, Vres_tmp.dim(1:2), 1);
      end
    end
    Yseg = Yseg_res;

    thickness_resampled = true;
  end
end

% add atrophy to GM by decreasing GM value in ROI and increasing CSF value
if simu_atrophy
  Yseg = simulate_atrophy(simu, Yseg, dim, template_dir, idef_name);
end

% get ground truth label using GM/WM/CSF
label_pve = zeros(dim, 'single');
order = [3 1 2];
for k = 1:3
    label_pve = label_pve + k*(Yseg(:,:,:,order(k)));
end

% load WMH map
if simu.WMH
  [WMH, res, label_pve] = simulate_WMHs(simu, res, label_pve, template_dir, idef_name);
else
  WMH = [];
end

% extend target voxel size if needed
if isscalar(requested_resolution)
  resolution_is_nan = isnan(requested_resolution);
else
  resolution_is_nan = isnan(requested_resolution(1));
end

if resolution_is_nan
  if thickness_resampled
    simu.resolution = vx_native;
    change_resolution = 1;
  else
    simu.resolution = vx;
    change_resolution = 0;
  end
else
  if isscalar(requested_resolution)
    simu.resolution = requested_resolution * ones(1,3);
  else
    simu.resolution = requested_resolution;
  end
  change_resolution = 1;
end

if any(simu.thickness)
  [label_pve, Yseg] = simulate_thickness(label_pve, simu, Yseg, dim, ...
        template_dir, idef_name, vx, V, order);
end

Ysimu = synthesize_from_segmentation(Yseg, name, res, mn, dim, WMH);

% apply either predefined MNI bias field or simulated bias filed before resizing
% to defined output resolution
if rf.percent ~= 0
  if ischar(rf.type)
    [Ysimu, rf_field] = add_bias_field(Ysimu, rf, idef_name, pth_root); % add predefined MNI field
  else
    [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf);
  end
end

mx_vol = max(Ysimu(:));

% output matrix
if resolution_is_nan && thickness_resampled
  Vout = V_native;
  vx_out = vx_native;
else
  Vout = V;
  vx_out = vx;
end

Vres.dim = round(Vout.dim.*vx_out./simu.resolution);
P = spm_imatrix(Vout.mat);
P(7:9) = P(7:9)./vx_out.*simu.resolution;
P(1:3) = P(1:3) + vx_out - simu.resolution;
Vres.mat = spm_matrix(P);

% output in defined resolution
volres   = zeros(Vres.dim);
labelres_pve = zeros(Vres.dim);

if change_resolution
  for sl = 1:Vres.dim(3)
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vres.mat\V.mat\M;
    
    % use sinc interpolation for simulated image
    volres(:,:,sl) = spm_slice_vol(Ysimu,M1,Vres.dim(1:2),-5);
    % and linear interpolation for label image
    labelres_pve(:,:,sl) = spm_slice_vol(label_pve,M1,Vres.dim(1:2),1);
  end
else % we can skip interpolation if voxels size is the same
  volres = Ysimu;
  labelres_pve = label_pve;
end

volres = volres / mx_vol;

% optionally ensure same noise for every trial
rng(simu.rng,'twister');

% Apply contrast change (power-law) if requested: normalize to [0,1], apply
% Y.^x, then rescale back to the original pre-transform min/max intensity.
if simu.contrast ~= 1
  vmin = min(volres(:));
  vmax = max(volres(:));
  if isfinite(vmin) && isfinite(vmax) && vmax > vmin && simu.contrast > 0
    v0 = (volres - vmin) / (vmax - vmin);
    v0 = max(min(v0,1),0);
    v0 = v0 .^ simu.contrast;
    volres = v0 * (vmax - vmin) + vmin;
  end
end

% Add noise:
% - If simu.snrWM>0: add Rician noise with target WM SNR
% - Else: fall back to Gaussian noise using percentage of WM mean
if simu.snrWM > 0
  % Desired SNR is defined for the underlying complex signal amplitude at WM
  % Compute complex noise std using absolute WM mean, then convert to normalized units
  sigma_abs = mn(2) / simu.snrWM;    % absolute-domain sigma (same units as Ysimu)
  sigma = sigma_abs / mx_vol;        % normalized-domain sigma
  % Complex noise components (n1,n2) are Gaussian. Magnitude combination yields
  % Rician noise (approx. Gaussian at high SNR).
  n1 = sigma * randn(size(volres));
  rng(simu.rng+1,'twister');
  n2 = sigma * randn(size(volres));
  volres = sqrt( (volres + n1).^2 + (n2).^2 );
else
  % Gaussian noise using percentage of WM mean intensity
  sigma_abs = (simu.pn/100) * mn(2); % absolute-domain std dev
  sigma = sigma_abs / mx_vol;        % normalized-domain std dev
  noise = sigma * randn(size(volres));
  volres = volres + noise;
end
% Clamp and rescale back
volres = mx_vol * max(min(volres, 1), 0);

mean_resolution = round(10*mean(simu.resolution));

if simu.snrWM > 0
  str1 = sprintf('snr%g',simu.snrWM);
elseif simu.pn > 0
  str1 = sprintf('pn%g',simu.pn);
else
  str1 = '';
end
if rf.percent ~= 0
  if isnumeric(rf.type)
    str1 = sprintf('%s_rf%g_%d', str1, rf.percent, rf.type(1));
  else
    str1 = sprintf('%s_rf%g_%s', str1, rf.percent, rf.type);
  end
end
if simu.contrast ~= 1
  switch  simu.contrast
    case 0.5
      str1 = sprintf('%s_conLow', str1);
    case 1.5
      str1 = sprintf('%s_conHigh', str1);
    otherwise
      str1 = sprintf('%s_con%g', str1, simu.contrast);
  end
end

if simu_atrophy
  if numel(simu.atrophy{2}) > 1
    str2 = sprintf('%s_multi',simu.atrophy{1});
  else
    str2 = sprintf('%s_%d_%g',simu.atrophy{1},simu.atrophy{2},simu.atrophy{3});
  end
else
  str2 = '';
end
if simu.WMH
  str2 = sprintf('%s_WMH%g', str2, simu.WMH);
end
if change_resolution
  str2 = sprintf('%s_res-%02dmm',str2,mean_resolution);
end
if any(simu.thickness)
  thickness = round(10*simu.thickness);
  if isscalar(simu.thickness)
    str2 = sprintf('%s_thickness%02dmm',str2,thickness);
  else
    str2 = sprintf('%s_thickness%02dmm-%02dmm',str2,min(thickness),max(thickness));
  end
end

if ~isempty(str2) & strcmp(str2(1),'_')
  if length(str2)>1
    str2 = str2(2:end);
  else
    str2 = '';
  end
end
str = strrep([str1 '_' str2 '_'],'__','_');
if strcmp(str(1),'_'), str = str(2:end); end
if contains(name_out, '_T1w')
  new_name = strrep(name_out,'_T1w',['_desc-' str 'T1w']);
  new_name_label = strrep(name_out,'_T1w',['_desc-' str2  '_label-seg']);
  new_name_bias = strrep(name_out,'_T1w',['_desc-' str2  '_RFfield']);
else
  new_name = [name_out '_desc-' str];
  new_name_label = [name_out '_desc-' str2 '_label-seg'];
  new_name_bias = [name_out '_desc-' str2 '_RFfield'];
end

% Remove remaining string issues
new_name_label = strrep(new_name_label,'desc-_','desc-');
new_name_bias = strrep(new_name_bias,'desc-_','desc-');

% write simulated image (optionally to derivatives folder)
simu_name = fullfile(out_pth, [new_name '.nii']); simu_name_main = simu_name;
fprintf('Save simulated image %s\n', simu_name);
Vres.fname = simu_name;
Vres.pinfo = [1 0 352]';
Vres.dt    = [spm_type('float32') 0];
spm_write_vol(Vres, single(volres));
if is_gz
  gzip(simu_name);
  spm_unlink(simu_name);
end

% write JSON sidecar with simulation parameters for both images
try
  gen.Name = 'mri_simulate';
  gen.Version = 'unknown';
  gen.SourceDatasets = {sprintf('%s%s', name_out, ext)};

  if simu.snrWM > 0
    SNRval = simu.snrWM;
    NoiseFrac = NaN; % not used when SNR is specified
  else
    SNRval = NaN;
    NoiseFrac = simu.pn/100;
  end

  if any(simu.thickness)
    if isscalar(simu.thickness)
      thickStr = sprintf('%.3gmm', simu.thickness);
    else
      thickStr = sprintf('%.3g-%.3gmm', min(simu.thickness), max(simu.thickness));
    end
  else
    thickStr = '';
  end

  if isfinite(NoiseFrac), simpar.NoiseFraction = NoiseFrac; end
  if isfinite(SNRval), simpar.SNR = SNRval; end
  if isfield(simu,'contrast') && ~isempty(simu.contrast) && simu.contrast ~= 1
    simpar.ContrastChange = simu.contrast;
  end
  simpar.VoxelSize = simu.resolution(:)';
  simpar.BiasFieldStrength = rf.percent;
  if rf.percent ~= 0
      simpar.BiasFieldType = rf.type(1);
  end
  if any(simu.thickness)
    simpar.Thickness = thickStr;
  end
  if simu_atrophy
    if numel(simu.atrophy{2}) > 1
      simpar.atrophy = sprintf('%s_multi',simu.atrophy{1});
    else
      simpar.atrophy = sprintf('%s_%d_%g',simu.atrophy{1},simu.atrophy{2},simu.atrophy{3});
    end
  end
  if simu.WMH
    simpar.WMHs = simu.WMH;
  end

  meta = struct();
  meta.GeneratedBy = {gen};
  meta.SimulationParameters = simpar;

  % write JSON next to main image using SPM's writer (handles NaN/null nicely)
  jsonMain = regexprep(simu_name_main,'\.nii(\.gz)?$','.json');
  spm_jsonwrite(jsonMain, meta);
catch ME
  fprintf('Warning: Failed to write JSON sidecar(s): %s\n', ME.message);
end

% write ground truth label
label_pve_name = fullfile(out_pth, [new_name_label '.nii']);
fprintf('Save %s\n', label_pve_name);
Vres.fname = label_pve_name;
Vres.pinfo = [1/255/3 0 352]';
Vres.dt    = [4 0];
spm_write_vol(Vres, labelres_pve);
if is_gz
  gzip(label_pve_name);
  spm_unlink(label_pve_name);
end

% save simulated bias field if defined
if rf.save
  rf_name = fullfile(out_pth, [new_name_bias '.nii']);
  fprintf('Save %s\n', rf_name);
  Vres.fname = rf_name;
  Vres.pinfo = [1 0 352]';
  Vres.dt    = [16 0];
  spm_write_vol(Vres, rf_field);
  if is_gz
    gzip(rf_name);
    spm_unlink(rf_name);
  end
end

% remove temporary files
if ~isempty(idef_name_orig) && exist(idef_name_orig,'file')
  spm_unlink(idef_name_orig);
end
if is_gz
  spm_unlink(simu.name);
end
if thickness_resampled && exist(resampled_name,'file')
  spm_unlink(resampled_name);
end

fprintf('================================================================================\n');

%==========================================================================
% function t = transf(B1,B2,B3,T)
% from spm_preproc_write8.m
%
% Purpose
%   Reconstruct a low-rank 3D field (e.g., bias field) using separable DCT
%   bases along x/y/z with coefficients T.
%
% Inputs
%   B1,B2,B3 - DCT basis matrices for x, y, and z.
%   T        - Coefficient tensor; if empty, a zero field is returned.
%
% Output
%   t        - Reconstructed field with size [size(B1,1) size(B2,1) size(B3,1)].
%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
  d2 = [size(T) 1];
  t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
  t  = B1*t1*B2';
else
  t = zeros(size(B1,1),size(B2,1),size(B3,1));
end
return


%==========================================================================
% function p = likelihoods(f,bf,mg,mn,vr)
% from spm_preproc_write8.m
%
% Purpose
%   Compute per-voxel likelihoods for a Gaussian mixture model given features
%   (optionally bias-corrected) and component parameters.
%
% Inputs
%   f   - cell array with one feature image; values are vectorized internally.
%   bf  - optional bias field multiplier (same shape as f{1}); [] for none.
%   mg  - mixture weights for each Gaussian component (1 x K).
%   mn  - means per component (N x K, N=#features).
%   vr  - covariance matrices per component (N x N x K).
%
% Output
%   p   - per-voxel likelihoods for each component (numel(f{1}) x K).
%==========================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N
  if isempty(bf)
    cr(:,n) = double(f{1}(:));
  else
    cr(:,n) = double(f{1}(:).*bf{1}(:));
  end
end

p  = ones(numel(f{1}),K);
for k=1:K
  amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
  d      = bsxfun(@minus,cr,mn(:,k)')/chol(vr(:,:,k));
  p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end
return

%==========================================================================
% function mask = is_in_atlas(atlas, regions)
% Purpose
%   Create mask of defined regions 
%
% Inputs
%   atlas - single(dims): atlas with integer values for defined regions.
%   regions - vector of integers that define regions to choose.
%
% Output
%   mask - logic(dim): mask of selected regions.
%==========================================================================
function mask = is_in_atlas(atlas, regions)

atlas = round(atlas);
mask = ismember(atlas, regions);

%==========================================================================
% function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, Vref, order)
%
% Purpose
%   Synthesize a constant cortical thickness by growing GM outward from the
%   original WM using an Euclidean distance map, and convert the resulting
%   hard labels into a partial-volume-like (PVE) segmentation by boundary
%   jittering and averaging.
%
% Inputs
%   label    - single(dims): Current PVE-like label image with values in [1..3]
%              (CSF=1, GM=2, WM=3). Used as baseline and to preserve subcortical
%              and cerebellar regions.
%   simu     - struct: Simulation options. Relevant fields:
%                .thickness: either a scalar (global thickness, in mm) or a
%                            3-element vector [occipital rest frontal], in mm.
%   Yseg     - single(dims,3): Current tissue probability/label volumes in the
%              order specified by 'order'. Will be overwritten by the simulated
%              PVE-like segmentation created here.
%   d        - [nx ny nz]: Volume dimensions.
%   template_dir - char: Path to CAT/SPM template directory (for Hammer atlas).
%   idef_name    - char: Filename of inverse deformation field to warp atlas
%                        into subject/native space with categorical interpolation.
%   vx       - [vx vy vz]: Voxel sizes in mm.
%   Vref     - struct: Target volume definition for resampling (same space as label).
%   order    - [3x1] int: Mapping from class index (CSF/GM/WM) to Yseg order.
%
% Outputs
%   label    - single(dims): New PVE-like label map in [1..3] after averaging
%              across boundary jitters, aligned with Yseg/order.
%   Yseg     - single(dims,3): Updated class volumes (CSF/GM/WM) representing
%              the simulated PVE-like segmentation.
%
% Algorithm
%   1) Atlas masks: Warp Hammer atlas to native space and build masks to exclude
%      subcortical/cerebellar regions from thickness manipulation. Optionally
%      define region masks to apply three distinct thickness values (occipital,
%      rest, frontal) when simu.thickness is a 3-vector.
%   2) Boundary jittering (PVE simulation): To emulate partial volume effects,
%      shift the label boundaries by 20 sub-voxel offsets uniformly spaced in
%      [-0.25, 0.25] voxels. For each offset:
%        a. Threshold labels to obtain a hard WM mask and clean it with simple
%           morphological steps.
%        b. Compute Euclidean distance transform from WM (compensated by 0.5*voxel).
%        c. Set CSF everywhere, then assign GM to voxels where D_WM <= thickness
%           (per region), preserving WM where present.
%        d. Convert the hard labels (1..3)
%   3) Average the 20 accumulated volumes to form a smooth PVE-like segmentation
%      and rebuild the final label image in [1..3] by weighted sum with class IDs.
%
% Notes
%   - Class encoding: CSF=1, GM=2, WM=3 throughout.
%   - When simu.thickness is scalar, the same thickness is applied globally.
%     When it has 3 values, they are applied to occipital, rest, and frontal
%     regions as defined by the Hammer atlas masks.
%   - This function modifies Yseg directly to reflect the new PVE-like tissue
%     maps, which are later used by the synthesis step to generate a T1 image.
%==========================================================================
function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, Vref, order)

Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
csf_val = 1; gm_val = 2; wm_val = 3;

% rescue original label for cerebellum/basal ganglia
label0 = label;

% warp atlas to native space using categorical interpolation
fprintf('Transform atlas to native space. This may take a while...\n');
atlas_name = fullfile(template_dir,'neuromorphometrics.nii');
atlas = cat_vol_defs(struct('field1',{{idef_name}},'images',{{atlas_name}},'interp',-1,'modulate',0));
atlas = single(atlas{1}{1});

% resample atlas to current grid if needed
if any(size(atlas) ~= d)
  Vdef = spm_vol(idef_name);
  Vdef = Vdef(1);
  atlas_res = zeros(Vref.dim, 'single');
  for sl = 1:Vref.dim(3)
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vref.mat\Vdef.mat\M;
    atlas_res(:,:,sl) = spm_slice_vol(atlas, M1, Vref.dim(1:2), 0);
  end
  atlas = atlas_res;
end

% create mask for non-cortical structures (i.e. basal ganglia, cerebellum, 
% hippocampus, amygdala
basal_ganglia = [31:32 36:37 47:48 55:62];
cerebellum = [11 35 38:41 71:73];
mask_orig = is_in_atlas(atlas, [basal_ganglia cerebellum]);
mask_orig = cat_vol_morph(mask_orig,'dd',1);

% soften mask borders to avoid hard transitions when reusing original labels
mask_soft = single(mask_orig);
spm_smooth(mask_soft, mask_soft, 1.5*vx);
mask_soft = min(max(mask_soft,0),1);

% create mask for mainly occipital and frontal areas
region1 = [108:109 114:115 128:129 134:135 144:145 156:161 170:171 176:177 182:183 196:197];
region3 = [31:32 102:105 120:121 132:133 136:137 146:147 152:155 162:165 172:173 178:179 190:191 202:205];

mask_thickness{1} = is_in_atlas(atlas, region1); % mainly occipital
mask_thickness{3} = is_in_atlas(atlas, region3); % mainly frontal
mask_thickness{2} = ~mask_thickness{1} & ~mask_thickness{3}; % remaining parts

mask = round(label) > 0;

% force stronger PVE effects by smoothing
spm_smooth(label,label,2.5*vx);

label1 = cell(numel(simu.thickness),1);

% save original segmentation to later include original cerebellum and basal
% ganglia
Yseg0 = Yseg(:,:,:,1:3);
Yseg(:,:,:,1:3) = 0;

% vary range of PVE from -0.25..0.25 in 15 steps to get more realistic PVE
% effects (optionally weighted)
pve_range = linspace(-0.25,0.25,15);

% apply gray closing to strengthen thin WM structures
label = cat_vol_morph(label,'gc',2);

for pve_step = 1:numel(pve_range)
  % define wm and remove disconnected regions
  wm  = round(label+pve_range(pve_step)) == wm_val;
  wm = cat_vol_morph(wm,'l',1, vx);

  % remove basal ganglia and cerebellum with soft blending
  wm(mask_orig) = 0;

  % euclidean distance to wm (CAT12 function if Image Toolbox is not available)
  if exist('bwdist')
    Dwm = (bwdist(wm) - 0.5) * mean(vx); % also consider voxelsize/2 correction of distance
  else
    Dwm = (cat_vbdist(single(wm)) - 0.5) * mean(vx); % also consider voxelsize/2 correction of distance
  end
  
  for k=1:numel(simu.thickness)
    label1{k} = round(label+pve_range(pve_step));

    label1{k}(~wm) = csf_val;
    label1{k}(~mask) = 0;
    
    % limit dilated gm to defined thickness
    label1{k}(label1{k} == csf_val & Dwm <= simu.thickness(k)) = gm_val; 
  end
    
  % replace tissue maps with modified label
  for j = 1:3
    if isscalar(simu.thickness)
      % only simulate 2mm thickness
      tmp_seg = single(round(label1{1}) == (j));
    else
      tmp_seg = Yseg(:,:,:,order(j));
      for k = 1:3
        tmp_seg(mask_thickness{k}) = single(round(label1{k}(mask_thickness{k})) == (j));
      end
    end

    Yseg(:,:,:,order(j)) = Yseg(:,:,:,order(j)) + tmp_seg/numel(pve_range);
  end
end

% update ground truth label
label = zeros(d, 'single');

for k = 1:3
  Yseg(:,:,:,order(k)) = Yseg(:,:,:,order(k)).*(1-mask_soft) + Yseg0(:,:,:,order(k)).*mask_soft;
  label = label + k*Yseg(:,:,:,order(k));
end


% close remaining holes in CSF
mask = label > 0.5;
label(mask ~= cat_vol_morph(mask,'dc',4)) = 1;


function [label, Yseg] = simulate_thickness_hammers(label, simu, Yseg, d, template_dir, idef_name, vx, Vref, order)

Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
csf_val = 1; gm_val = 2; wm_val = 3;

% rescue original label for cerebellum/basal ganglia
label0 = label;

% warp atlas to native space using categorical interpolation
fprintf('Transform atlas to native space. This may take a while...\n');
hammers_name = fullfile(template_dir,'hammers.nii');
hammers = cat_vol_defs(struct('field1',{{idef_name}},'images',{{hammers_name}},'interp',-1,'modulate',0));
hammers = single(hammers{1}{1});

% resample atlas to current grid if needed
if any(size(hammers) ~= d)
  Vdef = spm_vol(idef_name);
  Vdef = Vdef(1);
  hammers_res = zeros(Vref.dim, 'single');
  for sl = 1:Vref.dim(3)
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vref.mat\Vdef.mat\M;
    hammers_res(:,:,sl) = spm_slice_vol(hammers, M1, Vref.dim(1:2), 0);
  end
  hammers = hammers_res;
end

% create mask for basal ganglia and cerebellum
mask_orig = (hammers>0   & hammers<5)  | hammers==9  | hammers==10 | ...
            (hammers>16  & hammers<20) | (hammers>33 & hammers<50) | ...
             hammers==74 | hammers==75;
mask_orig = cat_vol_morph(mask_orig,'dd',1);

% soften mask borders to avoid hard transitions when reusing original labels
mask_soft = single(mask_orig);
spm_smooth(mask_soft, mask_soft, 1.5*vx);
mask_soft = min(max(mask_soft,0),1);

% create mask for occipital and frontal areas and remaining parts
mask_thickness{1} = (hammers > 63.5 & hammers < 67.5) | (hammers > 21.5 & hammers < 23.5); % occipital
mask_thickness{3} = (hammers > 55.5 & hammers < 59.5) | (hammers > 27.5 & hammers < 29.5) | (hammers > 68.5 & hammers < 71.5); % frontal

mask_thickness{3} = (hammers > 55.5 & hammers < 59.5) | (hammers > 12.5 & hammers < 14.5); % inferior/superior frontal + inferior/middle temporal
mask_thickness{2} = ~mask_thickness{1} & ~mask_thickness{3}; % remaining parts

mask = round(label) > 0;

% force stronger PVE effects by smoothing
spm_smooth(label,label,2.5*vx);

label1 = cell(numel(simu.thickness),1);

% save original segmentation to later include original cerebellum and basal
% ganglia
Yseg0 = Yseg(:,:,:,1:3);
Yseg(:,:,:,1:3) = 0;

% vary range of PVE from -0.25..0.25 in 15 steps to get more realistic PVE
% effects (optionally weighted)
pve_range = linspace(-0.25,0.25,15);

% apply gray closing to strengthen thin WM structures
label = cat_vol_morph(label,'gc',2);

for pve_step = 1:numel(pve_range)
  % define wm and remove disconnected regions
  wm  = round(label+pve_range(pve_step)) == wm_val;
  wm = cat_vol_morph(wm,'l',1, vx);

  % remove basal ganglia and cerebellum with soft blending
  wm(mask_orig) = 0;

  % euclidean distance to wm (CAT12 function if Image Toolbox is not available)
  if exist('bwdist')
    Dwm = (bwdist(wm) - 0.5) * mean(vx); % also consider voxelsize/2 correction of distance
  else
    Dwm = (cat_vbdist(single(wm)) - 0.5) * mean(vx); % also consider voxelsize/2 correction of distance
  end
  
  for k=1:numel(simu.thickness)
    label1{k} = round(label+pve_range(pve_step));

    label1{k}(~wm) = csf_val;
    label1{k}(~mask) = 0;
    
    % limit dilated gm to defined thickness
    label1{k}(label1{k} == csf_val & Dwm <= simu.thickness(k)) = gm_val; 
  end
    
  % replace tissue maps with modified label
  for j = 1:3
    if isscalar(simu.thickness)
      % only simulate 2mm thickness
      tmp_seg = single(round(label1{1}) == (j));
    else
      tmp_seg = Yseg(:,:,:,order(j));
      for k = 1:3
        tmp_seg(mask_thickness{k}) = single(round(label1{k}(mask_thickness{k})) == (j));
      end
    end

    Yseg(:,:,:,order(j)) = Yseg(:,:,:,order(j)) + tmp_seg/numel(pve_range);
  end
end

% update ground truth label
label = zeros(d, 'single');

for k = 1:3
  Yseg(:,:,:,order(k)) = Yseg(:,:,:,order(k)).*(1-mask_soft) + Yseg0(:,:,:,order(k)).*mask_soft;
  label = label + k*Yseg(:,:,:,order(k));
end


% close remaining holes in CSF
mask = label > 0.5;
label(mask ~= cat_vol_morph(mask,'dc',4)) = 1;


%==========================================================================
% function Yseg = simulate_atrophy(simu, Yseg, dims, template_dir, idef_name)
%
% Purpose
%   Apply regional atrophy by increasing CSF (and effectively reducing GM)
%   within specified atlas ROIs warped into native space.
%
% Inputs
%   simu         - struct with .atrophy = {atlasName, roiIds[], factors[]}
%   Yseg         - single(dims,3): tissue maps (CSF/GM/WM order as used).
%   dims         - [nx ny nz] dimensions of the volume.
%   template_dir - path to atlas templates; atlas is warped categorically.
%   idef_name    - inverse deformation field to native space.
%
% Output
%   Yseg         - updated tissue maps with CSF increased in target ROIs.
%==========================================================================
function Yseg = simulate_atrophy(simu, Yseg, dims, template_dir, idef_name)

% warp atlas to native spave using categorical interpolation
fprintf('Transform atlas to native space. This may take a while...\n');
atlas_name = fullfile(template_dir,[simu.atrophy{1} '.nii']);
atlas = cat_vol_defs(struct('field1',{{idef_name}},'images',{{atlas_name}},'interp',-1,'modulate',0));
atlas = single(atlas{1}{1});

% go through all defined ROIs
for i = 1:numel(simu.atrophy{2})
  ind_atlas = round(atlas) == simu.atrophy{2}(i);
  
  % check that ROI is not empty because it was not correct defined
  if isempty(ind_atlas)
    fprintf('ROI #%d seems to be wrong and does not exist in %s.\n', simu.atrophy{2}, simu.atrophy{1});
    return
  end
  
  % create atrophy mask with defined value in mask and otherwise 1
  mod_atlas = zeros(dims, 'single');
  mod_atlas(ind_atlas) = simu.atrophy{3}(i);
  mod_atlas(~ind_atlas) = 1.0;
  
  [~, maxind] = max(Yseg,[],4);
  GM_atlas = sum(maxind(ind_atlas) == 1);

  % increase CSF by factor defined in atrophy mask
  Yseg(:,:,:,3) = mod_atlas.*Yseg(:,:,:,3);
  
  [~, maxind] = max(Yseg,[],4);
  GM_atlas_simu = sum(maxind(ind_atlas) == 1);
  fprintf('GM reduction in region #%d = %3.2f%s\n',simu.atrophy{2}(i),100*(1-GM_atlas_simu/GM_atlas),'%');
end


%==========================================================================
% function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d, WMH)
%
% Purpose
%   Generate a T1-like image from provided CSF/GM/WM maps by inserting them
%   into SPM's mixture model and computing the expected intensity per voxel.
%
% Inputs
%   vol_seg - single(dims,3): tissue maps used in place of SPM posteriors.
%   name    - base name for progress display.
%   res     - struct from SPM segmentation (mg, mn, vr, Tbias, etc.).
%   mn      - 3x1 means for CSF/GM/WM replacing the corresponding mixture means.
%   d       - [nx ny nz] dimensions.
%   WMH     - optionally add white matter hyperintensities (WMHs)
%
% Output
%   Ysimu   - synthesized T1-weighted image volume (single).
%==========================================================================
function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d, WMH)
% go through all peaks that are defined
% mainly copied from spm_preproc_write8.m

K = size(res.mn,2);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3 = 1:d(3);

% prepare DCT parameters for bias correction
chan = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
d3      = [size(res.Tbias{1}) 1];
chan.B3 = spm_dctmtx(d(3),d3(3),x3);
chan.B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
chan.B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
chan.T  = res.Tbias{1};

% output image
Ysimu = zeros(d, 'single');

spm_progress_bar('init',length(x3),['Working on ' name],'Planes completed');
for z = 1:length(x3)

  % Bias corrected image
  f  = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
  bf = exp(transf(chan.B1,chan.B2,chan.B3(z,:),chan.T));
  cr{1} = bf.*f;

  msk = any((f==0) | ~isfinite(f),3);

  % Parametric representation of intensity distributions
  q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
  q1  = reshape(q1,[d(1:2),numel(res.mg)]);
  
  % replace classes and means for GM/WM/CSF by external segmentations
  for k=1:3
    ind = find(res.lkp==k);
    for j=1:numel(ind)
      q1(:,:,ind(j)) = vol_seg(:,:,z,k);
      res.mn(1,ind(j)) = mn(k);
    end
  end

  % add WMHs as last class
  if ~isempty(WMH)
    q1(:,:,K) = WMH(:,:,z);
  end
  
  s   = sum(q1,3);
  tmp = zeros(d(1:2));

  % sum up over first 3 classes and normalize to sum of 1
  for k=1:find(res.lkp<=3, 1, 'last')
    tmp = tmp + res.mg(k)*res.mn(1,k)*q1(:,:,k)./s;
  end

  % add remaining 3 BG classes from bias corrected image
  ind = tmp == 0;
  tmp(ind) = cr{1}(ind);
  tmp(msk) = 1e-3;

  Ysimu(:,:,z) = tmp;
  
  spm_progress_bar('set',z);
end

% Sometimes huge values occur due to bias correction in the noisy
% background and we have to limit values to 98% percentile
th = prctile(Ysimu(:),98);
Ysimu(Ysimu>th) = th;
spm_progress_bar('clear');


%==========================================================================
% function [Ysimu, rf_field]  = add_bias_field(Ysimu, rf, idef_name, pth)
%
% Purpose
%   Apply a predefined RF bias field (A/B/C) from template space to native
%   space, scaled by rf.percent, optionally invert for negative percent.
%
% Inputs
%   Ysimu     - simulated image.
%   rf        - struct with .percent (signed), .type ('A'|'B'|'C').
%   idef_name - inverse deformation field for warping the RF template.
%   pth       - path to directory containing rf100_*.nii fields.
%
% Outputs
%   Ysimu     - modulated image.
%   rf_field  - applied RF field in native space.
%==========================================================================
function [Ysimu, rf_field] = add_bias_field(Ysimu, rf, idef_name, pth)

fprintf('Transform RF field to native space.\n');
% warp defined rf field to native space
rf_name = fullfile(pth,['rf100_' rf.type '.nii']);
rf_field = cat_vol_defs(struct('field1',{{idef_name}},'images',{{rf_name}},'interp',1,'modulate',0));
rf_field = single(rf_field{1}{1});

% apply defined percent and strength
rf_field = abs(rf.percent)/100 * (single(rf_field));

% invert field for neg. values
if rf.percent < 0
  rf_field = 1./rf_field;
end

ind = isfinite(rf_field);
rf_field = 1 + rf_field - mean(rf_field(ind));
Ysimu(ind) = rf_field(ind).*Ysimu(ind);


%==========================================================================
% function [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf)
%
% Purpose
%   Create a smooth, random RF bias field via FFT-domain filtering with a
%   strength-dependent grid size, interpolate to image size, scale by rf.percent,
%   and apply to the simulated image.
%
% Inputs
%   Ysimu  - image to modulate.
%   rf     - struct: .type = [strength, rngSeed], .percent amplitude (signed).
%
% Outputs
%   Ysimu    - modulated image.
%   rf_field - generated RF field after smoothing and scaling.
%==========================================================================
function [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf)

dim = size(Ysimu);

% set seeds to defined value
rng(rf.type(2),'twister')

N = 2^round(rf.type(1)); % Define the size of the 3D field w.r.t. defined strength
fwhm = 30; % Smoothing size
pad  = 3*fwhm; % pad border to prevent smoothing issues at image borders

% Generate a random 3D field
field = rand(N, N, N);

% Apply 3D FFT to the field
fieldFFT = fftn(field);

% Create a 3D frequency filter to manipulate smoothness
[x, y, z] = meshgrid(-N/2:N/2-1, -N/2:N/2-1, -N/2:N/2-1);
radius = sqrt(x.^2 + y.^2 + z.^2); % Distance from the center in 3D
smoothnessFilter = exp(-radius./(N/8)); % Gaussian filter for smoothness

% Shift the FFT for filtering
fieldFFTShifted = fftshift(fieldFFT);

% Apply the smoothness filter in the frequency domain
filteredFFT = fieldFFTShifted .* smoothnessFilter;

% Shift back and apply inverse 3D FFT
filteredField = ifftn(ifftshift(filteredFFT));

% Original grid
[x, y, z] = ndgrid(1:N, 1:N, 1:N);

% New grid dimensions
[xq, yq, zq] = ndgrid(linspace(1, N, dim(1)), ...
                      linspace(1, N, dim(2)), ...
                      linspace(1, N, dim(3)));

% Interpolate using interpn
filteredField = interpn(x, y, z, filteredField, xq, yq, zq, 'linear');

% resize field and scale it to 0..1
filteredField = filteredField - min(filteredField(:));
filteredField = filteredField/max(filteredField(:));

% pad border and create background 0f 0.5 to prevent smoothing issues at image border
vol = 0.5*ones(dim(1)+2*pad,dim(2)+2*pad,dim(3)+2*pad);
vol(pad+1:dim(1)+pad,pad+1:dim(2)+pad,pad+1:dim(3)+pad) = filteredField;

% smooth and resize back
spm_smooth(vol,vol,[fwhm fwhm fwhm]);
rf_field = vol(pad+1:dim(1)+pad,pad+1:dim(2)+pad,pad+1:dim(3)+pad);

% scale to an amplitude and mean of 1
rf_field = rf_field - min(rf_field(:));
rf_field = rf_field/max(rf_field(:));

% apply defined percent
rf_field = abs(rf.percent)/100 * rf_field;

% invert field for neg. values
if rf.percent < 0
  rf_field = 1./rf_field;
end
rf_field = 1 + rf_field - mean(rf_field(:));

% and finally apply bias field
Ysimu = rf_field.*Ysimu;


%==========================================================================
% function [WMH, res, label_pve] = simulate_WMHs(simu, res, label_pve, template_dir, idef_name)
%
% Purpose
%   Add simulated white matter hyperintensities (WMHs) consistent with a
%   prior WMH probability map and subject anatomy. The routine warps a WMH
%   atlas to native space, modulates it by a random field and a user-defined
%   strength, and injects the resulting WMH class into the PVE label image
%   and mixture model used for synthesis.
%
% Inputs
%   simu         - struct with field .WMH (double): strength parameter (>1
%                  increases WMH expression; values close to 1 keep it mild).
%   res          - SPM segmentation result struct (contains image header and
%                  GMM parameters mg/mn/vr/lkp; mn is updated for WMH class).
%   label_pve    - single(dims): current PVE-like label image in [1..3] (CSF=1,
%                  GM=2, WM=3). Will be extended to [1..4] where 4 encodes WMH.
%   template_dir - path containing the WMH prior map 'cat_wmh_miccai2017.nii'.
%   idef_name    - inverse deformation field to warp the WMH atlas to native
%                  space (continuous interpolation for intensities).
%
% Outputs
%   WMH          - single(dims): simulated WMH map in [0..1].
%   res          - struct: mixture model with an added WMH class (mn updated).
%   label_pve    - single(dims): updated label image where WMH contributes as
%                  an additional class (value 4 contribution), later used for
%                  synthesis.
%
% Algorithm
%   1) Warp the MICCAI2017 WMH prior map to native space and lightly smooth it.
%   2) Create a random 3D field, resample it to image dimensions, threshold to
%      form spatial support for patchy WMH distribution.
%   3) Erode WM to ensure spacing from GM and constrain WMHs to deep WM.
%   4) Combine WMH prior^(1/(strength-0.8)) with random field and WM mask,
%      smooth and normalize to [0,1].
%   5) Update label_pve by adding 2*WMH and clipping to max class label 4,
%      thereby introducing an additional WMH class contribution.
%   6) Extend the GMM by adding a WMH class intensity equal to the GM mean
%
% Notes
%   - Class encoding after this step: CSF=1, GM=2, WM=3, WMH=4.
%   - The strength parameter shapes the prior map nonlinearly; larger values
%     emphasize regions with higher WMH probability.
%   - The WM erosion step helps avoid spuriously labeling GM/CSF boundaries
%     as WMH.
%==========================================================================
function [WMH, res, label_pve] = simulate_WMHs(simu, res, label_pve, template_dir, idef_name)

Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
strength = simu.WMH;
vx = sqrt(sum(res.image(1).mat(1:3,1:3).^2));

% warp WMH map to native space
fprintf('Transform WMH map to native space.\n');
WMH_name = fullfile(template_dir,'cat_wmh_miccai2017.nii');
WMH = cat_vol_defs(struct('field1',{{idef_name}},'images',{{WMH_name}},'interp',1,'modulate',0));
WMH = single(WMH{1}{1});
dim = size(WMH);

% slightly smooth WMH atlas
spm_smooth(WMH, WMH, 2./vx);

rng(simu.rng, 'twister')

N = 2^round(5); % Define the size of the initial 3D field

% Generate a random 3D field
fieldN = rand(N, N, N);

% Original grid
[x, y, z] = ndgrid(1:N, 1:N, 1:N);

% New grid dimensions
[xq, yq, zq] = ndgrid(linspace(1, N, dim(1)), ...
                      linspace(1, N, dim(2)), ...
                      linspace(1, N, dim(3)));

% Interpolate using interpn
field = interpn(x, y, z, fieldN, xq, yq, zq, 'linear');
field = single(field>0.7);

% find WM in original image and erode it to ensure some space to GM
WM = cat_vol_morph(Yp0toC(label_pve, 3) > 0.5,'de',3,vx);

% apply strength parameter to WM atlas and combine it with field and WM
WMH = WMH.^(1/(strength-0.8)).*field.^2.*WM;
spm_smooth(WMH, WMH, 2./vx);
WMH = WMH/max(WMH(:));

% correct PVE label and add additional WMH class (which should be a
% increased in intensities)
label_pve = label_pve + 5*WMH/strength.^0.75;
label_pve(label_pve > 4) = 4;

% use mean of GM for additional class
K   = numel(res.mg);
intensity_WMH = mean(res.mn(1,res.lkp==1));
K = K + 1;
res.mn(1,K) = intensity_WMH;


% function Yb = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th)
% SKULL_STRIP_APRG - Brain extraction via CAT12 APRG
%
% Purpose
%   Compute a robust brain mask using CAT12's Adaptive Probability Region-Growing
%   (APRG) method, leveraging SPM/CAT tissue posteriors and data-driven
%   intensity thresholds from the current image. This improves downstream PVE
%   label generation by restricting to brain voxels only.
%
% Inputs
%   Ysrc - single(dims): Source image used for threshold estimation (same space
%          as SPM/CAT outputs from preproc). Typically the first channel.
%   Ycls - 1x6 cell of uint8(dims): Tissue class posteriors (0..255) from
%          cat_spm_preproc_write8, in SPM/CAT convention (c1..c6).
%   res  - struct: SPM/CAT segmentation result with fields used here:
%              .mn  (means per class/component)
%              .mg  (mixture weights)
%              .lkp (lookup of class index per component, 1..6)
%   dim  - [nx ny nz]: Volume dimensions of the images.
%   T3th - 1x3 double: Intensity thresholds [CMth, CSFth, WMth] used by APRG
%          (derived from get_tissue_thresholds). CMth is a central-matter
%          threshold between GM and WM, CSFth approximates CSF mean, WMth is a
%          robust WM anchor.
%
% Output
%   Yb   - logical/single(dims): Brain mask estimated by APRG (1=brain, 0=non-brain).
%
% Algorithm (summary)
%   1) Estimate class-specific intensity anchors via weighted means of the GMM
%      parameters (clsint for CSF/GM/WM).
%   2) Build a 4D stack P(:,:,:,1..6) from Ycls (uint8 posteriors), as required
%      by CAT12’s APRG routine.
%   3) Derive robust thresholds for WM and a central matter (CM) value using the
%      WM posterior median within Ysrc, which mitigates issues in presence of WMHs.
%   4) Call cat_main_APRG(Ysrc, P, res, T3th) to obtain the brain mask.
%
% Notes
%   - Assumes T1-like contrast (WM > GM > CSF) for threshold reasoning; a median
%     WM intensity is used to be robust to WMHs.
%   - Requires CAT12 on the MATLAB path; uses cat_main_APRG.
%   - The output mask matches the input volume dimensions and orientation.
%==========================================================================
function Yb = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th)

% we need a 4D array for cat_main_APRG
P = zeros([dim 6],'uint8');
for i=1:6
  P(:,:,:,i) = Ycls{i};
end

res.isMP2RAGE = 0;
Yb = cat_main_APRG(Ysrc, P, res, T3th);


function Yseg = close_WM_GM_holes(Yseg, Ysrc, Ycorr, Ycls, Yy, res, vx_vol)
% CLOSE_WM_GM_HOLES - Fill WMHs in WM and correct WM and GM segmentation
%
% Purpose
%   Detect and close WMHs using CAT12 WMH detection
%
% Inputs
%   Yseg   - single(dims,3): current class volumes (CSF/GM/WM-like) derived 
%            from LAS-corrected intensities.
%   Ycls   - 1x6 cell of uint8(dims): SPM tissue posteriors (0..255).
%   res    - segmentation structure from SPM segmentation.
%   vx_vol - voxel size in mm
%
% Output
%   Yseg   - class volumes with WMHs reassigned to WM, leaving cortex and
%            and boundaries untouched as much as possible.

% we have to prepare some parameters for cat_main_updateSPM1639
global cat; cat_defaults;
job = cat;
job.extopts.inv_weighting = 0; job.extopts.verb = 0;
P = zeros([size(Ycls{1}) numel(Ycls)],'uint8');
for i=1:numel(Ycls), P(:,:,:,i) = Ycls{i}; end
clear Ycls;
tpm.dat = cell(6,1);
tpm.V = res.tpm;
for i=1:6
  tpm.dat{i} = spm_read_vols(res.tpm(i));
end
noise = 0.03;
res.image0 = res.image; res.ppe.affreg.skullstripped = 0; res.ppe.affreg.highBG = 0;
stime  = cat_io_cmd('');
stime2 = cat_io_cmd('');

[~,Ycls,Yb] = cat_main_updateSPM1639(Ysrc,P,Yy,tpm,job,res,stime,stime2);
[Yl1,Ycls] = cat_vol_partvol(Ycorr,Ycls,Yb,Yy,vx_vol,job.extopts,tpm.V,noise,job,false(size(Yb)));

NS = @(Ys,s) Ys==s | Ys==s+1; 
LAB = job.extopts.LAB;
Ynwmh = NS(Yl1,LAB.TH) | NS(Yl1,LAB.BG) | NS(Yl1,LAB.HC) | NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS);
Ynwmh = cat_vol_morph(cat_vol_morph( Ynwmh, 'dd', 8 , vx_vol),'dc',12 , vx_vol) & ...
        ~cat_vol_morph( NS(Yl1,LAB.VT), 'dd', 4 , vx_vol); 
Ywmh  = Ycls{7}>0 & ~Ynwmh;

% reassign these WMHs to WM in Yseg
Yseg(:,:,:,1) = min(Yseg(:,:,:,1), (1 - Ywmh));
Yseg(:,:,:,2) = max(Yseg(:,:,:,2), Ywmh);
Yseg(:,:,:,3) = min(Yseg(:,:,:,3), (1 - Ywmh));

% renormalize CSF/GM/WM per voxel to sum <=1 (simple clamp)
S = sum(Yseg,4);
mask = S>1;
if any(mask(:))
  for k=1:3
    tmp = Yseg(:,:,:,k);
    tmp(mask) = tmp(mask) ./ S(mask);
    Yseg(:,:,:,k) = tmp;
  end
end



function T3th = get_tissue_thresholds(Ysrc, Ycls, res)
% GET_TISSUE_THRESHOLDS - Estimate robust CSF/GM/WM thresholds
%
% Purpose
%   Compute three intensity anchors used by LAS/APRG steps:
%     - WMth: a robust white-matter threshold derived from the weighted WM
%       mean (from the GMM) and the median intensity within high WM posterior
%       voxels to mitigate WMH effects.
%     - CSFth: an approximation of the CSF intensity anchor from the GMM.
%     - CMth: a central-matter threshold between GM and WM that adapts to the
%       current image contrast.
%
% Inputs
%   Ysrc - single(dims): Source image used to compute medians within WM.
%   Ycls - 1x6 cell of uint8(dims): SPM/CAT posteriors (0..255), where
%          Ycls{2} corresponds to WM.
%   res  - struct: SPM/CAT segmentation result providing Gaussian mixture
%          parameters mn, mg, and lookup lkp.
%
% Output
%   T3th - 1x3 double: [CMth, CSFth, WMth] thresholds used by skull stripping
%          and LAS normalization routines.
%
% Algorithm
%   1) clsint(k): weighted class mean from the GMM for class k.
%   2) WMth: max(clsint(WM), median(Ysrc within high WM posterior)).
%   3) CMth: if intensities are inverted (CSF > WM), use CSF anchor; otherwise
%      use 2*GMmean - WMth, clipped not to exceed CSF anchor.
%   4) Return [CMth, CSFmean, WMth].
%
% Notes
%   - Class indices follow SPM/CAT convention (typically: 1=GM, 2=WM, 3=CSF).
%   - High WM posterior is defined with a conservative threshold (Ycls{2}>192
%     on 0..255 scale) to obtain a stable median.
%   - These anchors are designed for T1-like contrast but include a simple
%     inversion check (CSF > WM) for robustness.

clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

% Use median for WM threshold estimation to avoid problems in case of WMHs
WMth = double(max(clsint(2), cat_stat_nanmedian(Ysrc(Ycls{2}>192)))); 
if clsint(3)>clsint(2) % invers
  CMth = clsint(3); 
else
  CMth = min([clsint(1) - diff([clsint(1),WMth]), clsint(3)]);
end

T3th = double([CMth, clsint(1), WMth]);