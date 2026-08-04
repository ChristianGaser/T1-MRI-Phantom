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
%   Preprocessing with SPM12 segmentation is required for custom images. It
%   runs automatically with the Blaiotta head and neck TPM and its batch job
%   BlaiottaSegmentJob.m, both stored next to this file, and is cached in
%   <name>_seg8.mat. An existing seg8.mat is reused as it is, thus delete it
%   to segment data again that were segmented with another TPM.
%   Requirements: SPM12 or SPM25 with CAT >= 26 installed (cat_main_LASsimple
%   is required and is checked for at startup). No MATLAB toolboxes beyond
%   base MATLAB are needed: distances use CAT's cat_bwdist (and cat_vbdist
%   where the index of the nearest voxel is needed) so that the result does
%   not depend on the installed toolboxes, and parpool is only used for
%   multiple input images when the Parallel Computing Toolbox is installed.
%
%   The function also writes a JSON sidecar next to the main image
%   containing key simulation metadata (tool/version, voxel size, noise/SNR,
%   RF field parameters, thickness tags), using SPM's spm_jsonwrite.
%
%   Thickness/PVE pipeline (when simu.thickness is set):
%   - A constant cortical thickness (global or region-wise) is synthesized by
%     expanding GM outward from the original WM using a Euclidean distance map.
%   - To obtain partial-volume-like transitions, the label boundary is shifted
%     across 15 sub-voxel offsets in the range [-0.25, 0.25] voxels, simulating
%     realistic boundary uncertainty. Each offset yields a hard label image
%     (CSF=1, GM=2, WM=3), and the results are averaged to produce a smooth
%     PVE-like label map.
%   - Deep inside the WM, i.e. at a safe distance from the GM/WM boundary, the
%     original tissue fractions are blended back in. The hard labels alone
%     would give a WM fraction of exactly 1 there and thus a perfectly flat
%     WM, whereas the original fractions carry the local variation that the
%     simulation without thickness manipulation shows.
%   - The final PVE label map is then used to synthesize a T1 image from the
%     SPM segmentation model by replacing GM/WM/CSF class posteriors and using
%     their Gaussian mixture parameters (means, variances, weights).
%
% Syntax:
%   mri_simulate(simu, rf)
%
% Parameters:
%   simu (struct): Simulation parameters. Defaults are applied for missing fields.
%       - 'name' (char): T1-weighted input image filename, either .nii or
%         .nii.gz. A compressed image is uncompressed next to the original
%         before the segmentation is run, because SPM would otherwise name
%         the segmentation parameters <name>.nii_seg8.mat, and the copy is
%         removed again afterwards. The outputs of a compressed input are
%         written compressed as well. Default: '' (empty), which triggers an
%         interactive file selection dialog (T1 only).
%       - 'pn' (double): Percentage noise level (Gaussian) relative to the WM
%         mean intensity. Ignored if 'snrWM' > 0. Default: 0 (percent of WM).
%       - 'snrWM' (double): If >0, adds Rician noise at a user-defined SNR for
%         white matter. Uses the (noise-free) WM mean to compute the complex
%         noise sigma via sigma = WMmean / snrWM, and generates magnitude
%         Rician noise: sqrt((S + n1).^2 + n2.^2). Default: 40.
%       - 'rng' (double, NaN or []): Seed for the random number generator.
%         A fixed number gives the same noise for every image, which is useful
%         to compare simulations but means that a whole dataset shares one
%         noise pattern. NaN or [] derive the seed from the filename instead,
%         so that every image gets its own reproducible noise.
%         Default: 0.
%       - 'derivative' (logical): If true, save outputs under a BIDS-style
%         derivatives folder at the dataset root, using the pipeline name
%         'mri_simulate-<version>' ('mri_simulate_thickness-<version>' for
%         thickness simulations) and mirroring the subject/session path.
%         Example: root/derivatives/mri_simulate-0.10.1/sub-*/ses-*/[anat|...]/
%         Default: 1 (save into derivatives).
%       - 'contrast' (double): Power-law contrast-change exponent applied to the
%         simulated image intensities before the noise is added. It is
%         normalized to
%         [0,1], transformed as Y.^contrast, and rescaled back to its original
%         min/max range. Use values >1 to increase contrast, <1 to reduce.
%         Default: 1 (no change). Meaningful values to simulate contrast
%         are 0.5 (low contrast) and 1.5 (high contrast).
%       - 'resolution' (double or [x, y, z]): Spatial resolution of the
%         simulated image. Default: NaN (keep original resolution). If scalar,
%         it is applied to all three axes; if a 3-vector, each axis is set individually.
%       - 'closeWMHholes' (logical): Detect and fill WMHs in WM and correct
%         segmentations to obtain a clean simulated image without WMHs (which
%         later allows to simulate additional WMHs using the WMH option).
%         Default: 0 (disabled).
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
%         Neuromorphometrics atlas is used to define these areas, as well as
%         subcortical areas and the cerebellum, which are excluded from the
%         thickness simulation to obtain a more realistic MRI. Either thickness
%         or atrophy can be simulated.
%         Default: 0 (disabled).
%       - 'parpool' (integer): Number of workers (processors) used by the
%         parpool command if the Parallel Computing Toolbox is available and
%         more than one image is given. It is limited to the number of images.
%         Default: half of the available cores.
%
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
%   Simulated MRI image files based on the specified parameters and features.
%   Names follow the BIDS filename grammar, i.e. entity-value pairs followed
%   by a suffix, where all option tags are collected in one alphanumeric
%   desc label and the output resolution uses the res entity. <opts> covers
%   every option, i.e. the noise, the bias field, the contrast and the
%   anatomy, while <anatOpts> covers the anatomy alone, so that runs which
%   differ only in noise or bias field share one ground truth:
%     - Main simulated image:   <entities>[_res-<vx>mm]_desc-<opts>_T1w.nii
%     - Ground-truth PVE label: <entities>[_res-<vx>mm][_desc-<anatOpts>]_dseg.nii
%     - Optional RF bias field (simulated fields only, rf.save=1):
%                               <entities>[_res-<vx>mm]_desc-<anatOpts>Biasfield_T1w.nii
%     - JSON sidecar next to the simulated image and the label image
%   <vx> is the voxel size in mm with 'p' as decimal separator, e.g.
%   res-0p5mm. Anisotropic voxels are listed per axis (res-0p5x0p5x1p5mm).
%   When simu.derivative is set, a dataset_description.json is written to the
%   root of the derivatives pipeline folder, as BIDS requires for a valid
%   derivative dataset. Fully valid BIDS names are only possible for a BIDS
%   input; for any other input the entities and the suffix added here still
%   follow the specification, but the input basename is kept as it is.
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
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 40,...
%                     'resolution', NaN,...
%                     'thickness', [1.5 2.0 2.5]);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
%   Example 4 - Simulation with custom RF field and added WMHs (medium strength)
%       simu = struct('name', 'custom_t1.nii', 'snrWM', 40,...
%                     'resolution', NaN, 'WMH', 2);
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
%   Example 5 - Apply contrast change (power-law, high contrast)
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 40, ...
%                     'resolution', NaN, 'contrast', 2);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
%
%
% TODO: simulation of motion artefacts using FFT and shift of phase information

% named tool_version and not version to not shadow the MATLAB builtin version()
tool_version = '0.10.1';

if ~exist('cat_main_LASsimple','file')
  error('Please update to a newer version >=CAT26 to use mri_simulate')
end

% Default simulation parameters
def.name       = '';
def.pn         = 0;
def.resolution = NaN;
def.WMH        = 0;
def.atrophy    = [];
def.thickness  = 0;
def.rng        = 0;
def.snrWM      = 40;
def.contrast   = 1;    % power-law contrast change exponent (1 = unchanged)
def.derivative = 1;    % save outputs into BIDS derivatives
def.closeWMHholes = 0; % don't close WMHs inside deep WM
def.parpool = feature('numcores')/2; % use half of the available processors

% An empty rng means the same as NaN here, i.e. seed from the filename. It has
% to be translated before cat_io_checkinopt, which drops empty fields and would
% thus silently replace it by the default seed.
if nargin > 0 && isstruct(simu) && isfield(simu,'rng') && isempty(simu.rng)
  simu.rng = NaN;
end

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
  % set number of workers/processors if Parallel Computing Toolbox is available
  if exist('parpool','file')
    parpool(min(simu.parpool, n_images));
  end
  P = simu.name;
  parfor i=1:n_images
    simu_i = simu;
    simu_i.name = deblank(P(i,:));
    mri_simulate(simu_i, rf);
  end
  if exist('parpool','file')
    delete(gcp('nocreate'))
  end
  return
end

[pth, name, ext, ~] = spm_fileparts(simu.name);
if ~exist(fullfile(pth, [name ext]), 'file')
  fprintf("File %s not found.\n", simu.name);
  return
end

% The interactive selection returns the frame number of the image, i.e.
% 'sub-01_T1w.nii,1'. spm_fileparts splits it off above, thus rebuilding the
% name from its parts removes it from everything that follows.
simu.name = fullfile(pth, [name ext]);

% Name of the input as it was given, thus still with .gz. It is only needed
% for the BIDS URI of the source in the JSON sidecars, which has to name the
% file of the dataset and not the uncompressed copy created below.
source_name = [name ext];

% Uncompress a gzipped input before anything else.
%
% SPM can read a .nii.gz, spm_vol uncompresses it into a temporary file, but
% the segmentation must not run on it: SPM names the segmentation parameters
% after the input file without its extension, so a compressed input is saved
% as <name>.nii_seg8.mat, and the inverse deformation field and everything
% else that is derived from the file name is malformed in the same way.
% Everything below therefore works on an uncompressed copy, and the outputs
% are compressed again at the end.
%
% An uncompressed file that is already there is used as it is and is never
% overwritten. It can be the file of the user and not a leftover of a former
% run, and gunzip would replace it without a word.
is_gz = strcmpi(ext,'.gz');   % spm_vol accepts .GZ as well, thus so do we
if is_gz
  nii_name = fullfile(pth, name);   % name still carries the .nii of a .nii.gz
  if exist(nii_name,'file')
    fprintf('Uncompressed %s exists already and is used as it is.\n', nii_name);
    simu.name = nii_name;
  else
    fname = gunzip(fullfile(pth, [name ext]));
    simu.name = fname{1};

    % Only a copy that was created here may be removed again. That has to
    % happen on every way out of the function, including the parameter checks
    % below that return early and an error inside the simulation, thus it is
    % left to onCleanup instead of to the cleanup at the end.
    gz_temp = simu.name;
    cleanup_gz = onCleanup(@() spm_unlink(gz_temp));
  end
  [pth, name, ext] = spm_fileparts(simu.name);
end

% keep output names clean if a temporary resample suffix is present
name_out = regexprep(name,'_thicknessRes05$','');

% If simu.rng is not defined we use the filename to create a seed. The test
% has to short-circuit: for an empty rng, 'isempty(x) | isnan(x)' is the empty
% array and not true, because isnan([]) is empty, and 'if []' does not run.
if isempty(simu.rng) || isnan(simu.rng)
  simu.rng = sum(double(name));
end

pth_root = fileparts(which(mfilename));

% Determine output folder (optionally BIDS derivatives)
out_pth = pth;
source_uri = '';        % BIDS URI of the input image, if it can be derived
if isfield(simu,'derivative') && simu.derivative
  try
    if any(simu.thickness)
      pipeline_name = ['mri_simulate_thickness-' tool_version];
    else
      pipeline_name = ['mri_simulate-' tool_version];
    end

    parts = strsplit(pth, filesep);
    if isempty(parts{1}), parts{1} = filesep; end
    idx_sub = find(strncmp(parts,'sub-',4), 1, 'first');
    if ~isempty(idx_sub) && idx_sub > 1
      % Derivatives at dataset root
      root_dir = fullfile(parts{1:idx_sub-1});
      rel_parts = parts(idx_sub:end); % sub-..[/ses-..]/anat/...
      pipeline_dir = fullfile(root_dir, 'derivatives', pipeline_name);
      out_pth = fullfile(pipeline_dir, rel_parts{:});
      % BIDS URI of the source image relative to the raw dataset root, which
      % always uses '/' as separator independent of the platform
      source_uri = ['bids::' strjoin([rel_parts, {source_name}], '/')];
    else
      % Fallback: place outputs under derivatives without mirroring
      root_dir = fileparts(pth);
      pipeline_dir = fullfile(root_dir, 'derivatives', pipeline_name);
      out_pth = pipeline_dir;
    end
    if ~exist(out_pth,'dir'), mkdir(out_pth); end

    % A BIDS derivative dataset requires a dataset_description.json at the
    % root of the pipeline directory, otherwise the dataset is not valid.
    write_dataset_description(pipeline_dir, pipeline_name, tool_version);
  catch
    % In case of any issue, fall back to input folder
    out_pth = pth;
    source_uri = '';
  end
end

% name of seg8.mat file that contains SPM12 segmentation parameters
mat_name = fullfile(pth, [name '_seg8.mat']);

% CAT12 template dir is later used for defining atrophy atlas
template_dir = fullfile(spm('dir'),'toolbox','CAT','templates_MNI152NLin2009cAsym');

% Call SPM segmentation if necessary and only keep the seg8.mat file.
%
% The segmentation uses the Blaiotta head and neck TPM and the batch job that
% belongs to it, both stored next to this file. That TPM is considerably more
% reliable than the SPM default one. The job is taken as it is instead of
% being rebuilt here, because it carries the settings of the validation work
% of the paper describing the TPM, and in particular the number of Gaussians
% per tissue, which describes that TPM and not an arbitrary one.
if ~exist(mat_name,'file')
  job_name = fullfile(pth_root,'BlaiottaSegmentJob.m');
  if ~exist(job_name,'file')
    error('Segmentation job %s was not found.', job_name);
  end
  matlabbatch = load_batch_job(job_name);

  matlabbatch{1}.spm.spatial.preproc.channel.vols = {simu.name};

  % The job defines the TPM itself, and a missing one would otherwise only
  % show up as an error of spm_vol deep inside the segmentation.
  tpm_name = spm_file(matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm{1},'number','');
  if ~exist(tpm_name,'file')
    error('TPM %s used by %s was not found.', tpm_name, job_name);
  end

  n_tissue = numel(matlabbatch{1}.spm.spatial.preproc.tissue);
  spm_jobman('run',matlabbatch);
  clear matlabbatch

  % remove native segmentations in case the job was changed to write them
  for i=1:n_tissue
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
  fprintf('Predefined MNI bias fields cannot be saved.\n');
  rf.save = 0;
end

% no bias field is created for an amplitude of 0 and thus nothing can be saved
if rf.percent == 0 && rf.save
  fprintf('No bias field is simulated for rf.percent = 0 and it cannot be saved.\n');
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

% The seg8.mat stores the TPM with the absolute path of the machine it was
% created on, so it is not found if the data were segmented elsewhere. Look
% for the same file name in the local SPM and CAT template folders instead.
res.tpm = find_missing_tpm(res.tpm);

if size(res.mn,1) > 1
  fprintf('Multi-modal segmentation is not recommended! Try again to segment the image using T1w data only.\n');
end

% Intensities of GM/WM/CSF from the Gaussian mixture of the segmentation.
%
% GM and WM use the mixing weighted mean over all their Gaussians. CSF must
% not: the segmentation models CSF with two Gaussians, and only the darker one
% describes CSF, while the brighter one quite often covers GM. Their weighted
% mean is therefore far too bright, and because this value is both the CSF
% intensity of the synthesis and the low anchor of the LAS correction and the
% skull stripping (see get_tissue_thresholds), it destabilizes the whole CSF
% segmentation. Only the darkest CSF Gaussian is used instead, which for a
% single Gaussian is the class mean and thus leaves that case unchanged.
csf_cls = 3;   % SPM class order, i.e. 1 = GM, 2 = WM, 3 = CSF
mn = zeros(3,1);
for k=1:3
  ind = find(res.lkp==k);
  if k == csf_cls
    mn(k) = min(res.mn(1,ind));
  else
    mg_k  = reshape(res.mg(ind),1,[]);   % row, whatever orientation mg has
    mn(k) = sum(mg_k .* res.mn(1,ind));
  end
end

if sum(res.lkp==csf_cls) < 2
  fprintf(['Warning: the segmentation describes CSF by a single Gaussian, thus the ' ...
           'darker CSF peak cannot be separated from the brighter one that often covers ' ...
           'GM.\n  Delete %s and run again to segment with %s, which uses two ' ...
           'Gaussians for CSF.\n'], mat_name, fullfile(pth_root,'BlaiottaSegmentJob.m'));
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
T3th = get_tissue_thresholds(Ysrc, Ycls, mn);

% slightly correct CSF-threshold since it's often underestimated
T3th(1) = 1.1*T3th(1);

% LAS correction and SANLM denoising
Ycorr = cat_main_LASsimple(Ysrc, Ycls, T3th);
cat_sanlm(Ycorr,3,1);

% Replace GM/WM/CSF segmentation by labels using LAS corrected image
Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
% single is sufficient for tissue probabilities in [0,1] and halves the memory,
% which matters most for parfor where every worker holds its own copy (a 0.5mm
% volume needs ~0.7GB instead of ~1.4GB here)
Yseg = zeros([dim, 3], 'single');

% Use CAT12 adaptive probability region-growing (APRG) approach for
% skull-stripping (uses T3th anchors; see skull_strip_APRG)
brainmask = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th);
Ycorr = Ycorr.*brainmask;

seg_order = [2 3 1];
for i = 1:3
    Yseg(:,:,:,i) = Yp0toC(3*Ycorr, seg_order(i)); 
end

% ensure that sum of Yseg is 1
Ysum = sum(Yseg,4)+eps;
for i = 1:3
    Yseg(:,:,:,i) = Yseg(:,:,:,i)./Ysum; 
end
clear Ysum

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
  Yseg = simulate_atrophy(simu, Yseg, dim, template_dir, idef_name, V);
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

% apply either predefined MNI bias field or simulated bias field before resizing
% to defined output resolution
if rf.percent ~= 0
  if ischar(rf.type)
    [Ysimu, rf_field] = add_bias_field(Ysimu, rf, idef_name, pth_root); % add predefined MNI field
  else
    [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf, vx);
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
  if rf.save, rfres = zeros(Vres.dim); end
  for sl = 1:Vres.dim(3)
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vres.mat\V.mat\M;

    % use sinc interpolation for simulated image
    volres(:,:,sl) = spm_slice_vol(Ysimu,M1,Vres.dim(1:2),-5);
    % and linear interpolation for label image
    labelres_pve(:,:,sl) = spm_slice_vol(label_pve,M1,Vres.dim(1:2),1);
    % the bias field has to follow the same grid to be saved with Vres
    if rf.save
      rfres(:,:,sl) = spm_slice_vol(rf_field,M1,Vres.dim(1:2),1);
    end
  end
else % we can skip interpolation if voxels size is the same
  volres = Ysimu;
  labelres_pve = label_pve;
  if rf.save, rfres = rf_field; end
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
  % both components are drawn from the same stream. Reseeding in between with
  % simu.rng+1 would make the n1 of one run identical to the n2 of the run
  % with the preceding seed, i.e. the two runs would be correlated.
  n1 = sigma * randn(size(volres));
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

%--------------------------------------------------------------------------
% Build BIDS-compatible output names.
%
% A BIDS filename is a sequence of entity-value pairs followed by a suffix:
%   <entity>-<label>[_<entity>-<label>...]_<suffix>.<extension>
% Entity labels must consist of letters and digits only, thus all option
% tags are concatenated into a single camelCase label for the desc entity
% (bids_label removes/replaces everything else). The output resolution uses
% the standard res entity instead and therefore stays outside of desc.
% Suffixes are T1w for the simulated image and dseg for the label image.
%--------------------------------------------------------------------------

% tags describing the acquisition: noise, bias field and contrast
desc_acq = '';
if simu.snrWM > 0
  desc_acq = sprintf('snr%g',simu.snrWM);
elseif simu.pn > 0
  desc_acq = sprintf('pn%g',simu.pn);
end
if rf.percent ~= 0
  % the sign of rf.percent selects an inverted field and has to be kept,
  % but a minus sign is not allowed inside a BIDS label
  if rf.percent < 0, rf_sign = 'Neg'; else, rf_sign = ''; end
  if isnumeric(rf.type)
    % the T separates the strength from the amplitude, which would otherwise
    % be concatenated into an ambiguous number (e.g. 30 and 2 -> 302)
    rf_str = sprintf('T%d', rf.type(1));
  else
    rf_str = rf.type;
  end
  desc_acq = sprintf('%sRf%s%g%s', desc_acq, rf_sign, abs(rf.percent), rf_str);
end
if simu.contrast ~= 1
  switch  simu.contrast
    case 0.5
      desc_acq = sprintf('%sConLow', desc_acq);
    case 1.5
      desc_acq = sprintf('%sConHigh', desc_acq);
    otherwise
      desc_acq = sprintf('%sCon%g', desc_acq, simu.contrast);
  end
end

% tags describing the anatomy: atrophy, WMHs and thickness. These do not
% depend on the noise or the bias field and are reused for the label image.
desc_anat = '';
if simu_atrophy
  if numel(simu.atrophy{2}) > 1
    desc_anat = sprintf('%sMulti',simu.atrophy{1});
  else
    % Roi/F separate the region id from the atrophy factor, which would
    % otherwise be concatenated into an ambiguous number
    desc_anat = sprintf('%sRoi%dF%g',simu.atrophy{1},simu.atrophy{2},simu.atrophy{3});
  end
end
if simu.WMH
  desc_anat = sprintf('%sWmh%g', desc_anat, simu.WMH);
end
if any(simu.thickness)
  thickness = round(10*simu.thickness);
  if isscalar(simu.thickness)
    desc_anat = sprintf('%sThickness%02dmm',desc_anat,thickness);
  else
    desc_anat = sprintf('%sThickness%02dto%02dmm',desc_anat,min(thickness),max(thickness));
  end
end

% the simulated image is described by both groups of tags
desc_main = bids_label([desc_acq desc_anat]);
desc_anat = bids_label(desc_anat);

% Without any option (no noise, no bias field, no contrast change and no
% anatomical modification) the desc entity would be empty and the simulated
% image could end up with the very name of the input image, so a minimal
% label is used to always keep the two apart.
if isempty(desc_main), desc_main = 'simu'; end

% Entities of the input name without the T1w suffix. A desc entity of the
% input is dropped, because desc must not occur twice and the new one
% describes this simulation.
if endsWith(name_out,'_T1w')
  bids_prefix = regexprep(name_out(1:end-4),'_desc-[a-zA-Z0-9]+$','');
else
  % not a BIDS input, so the whole name is kept as a prefix and only the
  % entities and the suffix added here can follow the specification
  bids_prefix = name_out;
end

% res is a standard BIDS derivative entity and precedes desc. The voxel size
% is given in mm with 'p' as decimal separator, because a BIDS label must be
% alphanumeric (0.5mm -> res-0p5mm). Anisotropic voxels are listed per axis
% instead of being averaged, which would both hide the anisotropy and report
% a size that no axis actually has.
if change_resolution
  res_lab = arrayfun(@(x) bids_label(sprintf('%g',x)), simu.resolution(:)', ...
                     'UniformOutput', false);
  if all(abs(simu.resolution - simu.resolution(1)) < 1e-6)
    ent_res = ['_res-' res_lab{1} 'mm'];
  else
    ent_res = ['_res-' strjoin(res_lab,'x') 'mm'];
  end
else
  ent_res = '';
end

if isempty(desc_main), ent_main = ''; else, ent_main = ['_desc-' desc_main]; end
if isempty(desc_anat), ent_anat = ''; else, ent_anat = ['_desc-' desc_anat]; end

% there is no registered BIDS suffix for a bias field, so it is written as a
% described variant of the T1w image
if isempty(desc_anat)
  ent_bias = '_desc-biasfield';
else
  ent_bias = ['_desc-' desc_anat 'Biasfield'];
end

new_name       = [bids_prefix ent_res ent_main '_T1w'];
new_name_label = [bids_prefix ent_res ent_anat '_dseg'];
new_name_bias  = [bids_prefix ent_res ent_bias '_T1w'];

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

% write JSON sidecar with simulation parameters
try
  % GeneratedBy only carries the tool itself; the input image is referenced
  % with the top level Sources field, which BIDS expects to hold BIDS URIs
  gen = struct();
  gen.Name = 'mri_simulate';
  gen.Version = tool_version;

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

  simpar = struct();
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
  if ~isempty(source_uri), meta.Sources = {source_uri}; end
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

% sidecar for the label image: the dseg suffix normally implies integer
% labels, so the partial volume encoding has to be documented here
try
  lmeta = struct();
  lmeta.GeneratedBy = {gen};
  if ~isempty(source_uri), lmeta.Sources = {source_uri}; end
  lmeta.Manual = false;
  lmeta.Description = ['Ground truth partial volume label image. Values are ' ...
                       'continuous and interpolate between the tissue labels, ' ...
                       'e.g. 2.5 is an equal mixture of GM and WM.'];
  if simu.WMH
    lmeta.Labels = struct('CSF',1,'GM',2,'WM',3,'WMH',4);
  else
    lmeta.Labels = struct('CSF',1,'GM',2,'WM',3);
  end
  spm_jsonwrite(regexprep(label_pve_name,'\.nii(\.gz)?$','.json'), lmeta);
catch ME
  fprintf('Warning: Failed to write label JSON sidecar: %s\n', ME.message);
end

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
  spm_write_vol(Vres, single(rfres));
  if is_gz
    gzip(rf_name);
    spm_unlink(rf_name);
  end
end

% Remove temporary files. The uncompressed copy of a gzipped input is not
% among them, it is removed by its onCleanup above, which also covers the ways
% out of the function that never arrive here.
if ~isempty(idef_name_orig) && exist(idef_name_orig,'file')
  spm_unlink(idef_name_orig);
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
% function Vtpm = find_missing_tpm(Vtpm)
%
% Purpose
%   Repair the TPM file names of a seg8.mat. The segmentation stores the TPM
%   as spm_vol structs with the absolute path of the machine it was created
%   on, so the file is not found when the data were segmented on a different
%   system or when SPM was installed elsewhere. In that case the same file
%   name is looked up in the local template folders.
%
% Inputs
%   Vtpm - spm_vol struct array of the TPM, one entry per tissue class, all
%          usually pointing to the same 4d file with a different volume
%          index n.
%
% Output
%   Vtpm - the same structure, with the file names of the local copy where
%          the original ones do not exist. Entries that cannot be resolved
%          are returned unchanged, so that the later read fails with its own
%          message rather than silently using a wrong template.
%
% Note
%   The volume index n is preserved per entry, and the structs are rebuilt
%   with spm_vol so that the cached nifti handle in the private field points
%   to the new file as well.
%==========================================================================
function Vtpm = find_missing_tpm(Vtpm)

if isempty(Vtpm) || ~isstruct(Vtpm) || ~isfield(Vtpm,'fname'), return; end

% search order: the tpm folder of SPM, then the CAT templates, finally
% mri_simulate folder
tpm_dirs = { fullfile(spm('dir'),'tpm'), ...
             fullfile(spm('dir'),'toolbox','CAT','templates_MNI152NLin2009cAsym'), ...;
             spm_fileparts(which('mri_simulate')) };
resolved = struct('old',{},'new',{});   % cache, all entries share one file

for i = 1:numel(Vtpm)
  if exist(Vtpm(i).fname,'file'), continue; end

  % reuse the result for a file name that was already looked up
  j = find(strcmp({resolved.old}, Vtpm(i).fname), 1);
  if ~isempty(j)
    newfile = resolved(j).new;
  else
    [~, nam, ext] = spm_fileparts(Vtpm(i).fname);
    newfile = '';
    for k = 1:numel(tpm_dirs)
      candidate = fullfile(tpm_dirs{k}, [nam ext]);
      if exist(candidate,'file'), newfile = candidate; break; end
    end
    resolved(end+1) = struct('old', Vtpm(i).fname, 'new', newfile); %#ok<AGROW>
    if isempty(newfile)
      fprintf(['Warning: TPM %s of the seg8.mat was not found and no copy ' ...
               'exists in\n  %s\n  %s\n'], Vtpm(i).fname, tpm_dirs{:});
    else
      fprintf('TPM %s not found, using %s instead.\n', Vtpm(i).fname, newfile);
    end
  end
  if isempty(newfile), continue; end

  % rebuild the entry and keep its volume index
  try
    Vtpm(i) = spm_vol(sprintf('%s,%d', newfile, Vtpm(i).n(1)));
  catch ME
    fprintf('Warning: could not read %s: %s\n', newfile, ME.message);
  end
end

%==========================================================================
% function matlabbatch = load_batch_job(job_name)
%
% Purpose
%   Return the batch job that a job script defines. The script is executed in
%   the workspace of this function and not in the one of the caller, thus the
%   variables it needs besides matlabbatch cannot collide with anything there.
%
% Inputs
%   job_name - char: filename of the batch job script.
%
% Output
%   matlabbatch - cell: the batch job that the script defines.
%==========================================================================
function matlabbatch = load_batch_job(job_name)

matlabbatch = {};
run(job_name);

%==========================================================================
% function label = bids_label(str)
%
% Purpose
%   Turn a tag into a valid BIDS entity label. BIDS only allows letters and
%   digits inside an entity value, so that a label like 'thickness15mm-25mm'
%   or 'con1.3' would make the filename invalid.
%
% Inputs
%   str - char: arbitrary tag.
%
% Output
%   label - char: str with decimal points replaced by 'p' (1.3 -> 1p3) and
%           all remaining non-alphanumeric characters removed.
%==========================================================================
function label = bids_label(str)

label = strrep(str, '.', 'p');
label = regexprep(label, '[^a-zA-Z0-9]', '');

%==========================================================================
% function write_dataset_description(pipeline_dir, pipeline_name, tool_version)
%
% Purpose
%   Write the dataset_description.json that BIDS requires at the root of a
%   derivative dataset. Without it the derivatives folder is not a valid
%   BIDS dataset. An existing file is never overwritten, so that a manually
%   edited description survives.
%
% Inputs
%   pipeline_dir  - char: root folder of the derivative pipeline.
%   pipeline_name - char: name of the pipeline (used as dataset Name).
%   tool_version  - char: version of mri_simulate.
%==========================================================================
function write_dataset_description(pipeline_dir, pipeline_name, tool_version)

dd_name = fullfile(pipeline_dir, 'dataset_description.json');
if exist(dd_name,'file'), return; end
if ~exist(pipeline_dir,'dir'), mkdir(pipeline_dir); end

dd = struct();
dd.Name        = pipeline_name;
dd.BIDSVersion = '1.8.0';
dd.DatasetType = 'derivative';
dd.GeneratedBy = {struct('Name','mri_simulate','Version',tool_version)};

spm_jsonwrite(dd_name, dd);

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
%   template_dir - char: Path to CAT/SPM template directory (for Neuromorphometrics atlas).
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
%   1) Atlas masks: Warp Neuromorphometrics atlas to native space and build masks to exclude
%      subcortical/cerebellar regions from thickness manipulation. Optionally
%      define region masks to apply three distinct thickness values (occipital,
%      rest, frontal) when simu.thickness is a 3-vector.
%   2) Boundary jittering (PVE simulation): To emulate partial volume effects,
%      shift the label boundaries by 15 sub-voxel offsets uniformly spaced in
%      [-0.25, 0.25] voxels. For each offset:
%        a. Threshold labels to obtain a hard WM mask and clean it with simple
%           morphological steps.
%        b. Compute Euclidean distance transform from WM (compensated by 0.5*voxel).
%        c. Set CSF everywhere, then assign GM to voxels where D_WM <= thickness
%           (per region), preserving WM where present.
%        d. Convert the hard labels (1..3)
%   3) Average the 15 accumulated volumes to form a smooth PVE-like segmentation
%      and rebuild the final label image in [1..3] by weighted sum with class IDs.
%   4) Blend the original tissue fractions back in, both inside the excluded
%      structures and inside the WM interior. The latter restores the local
%      fluctuations of the original segmentation that the hard labels and the
%      grey closing removed, and thus the intensity texture of the WM, without
%      touching the GM/WM boundary that defines the simulated thickness.
%
% Notes
%   - Class encoding: CSF=1, GM=2, WM=3 throughout.
%   - When simu.thickness is scalar, the same thickness is applied globally.
%     When it has 3 values, they are applied to occipital, rest, and frontal
%     regions as defined by the Neuromorphometrics atlas masks.
%   - This function modifies Yseg directly to reflect the new PVE-like tissue
%     maps, which are later used by the synthesis step to generate a T1 image.
%==========================================================================
function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, Vref, order)

csf_val = 1; gm_val = 2; wm_val = 3;

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

% Remove wm and fill it with the label of the nearest remaining structure.
% cat_vbdist and not the faster cat_bwdist, because the index of the nearest
% object voxel is needed here and only the vector propagation of cat_vbdist
% can provide it.
atlas(round(label) == wm_val) = 0;
[~, Yind] = cat_vbdist( single(atlas>0) );
atlas = atlas(Yind);

% In the neuromorphometrics atlas every cortical parcel has an id >= 100,
% while all non-cortical structures (ventricles, subcortical grey matter,
% cerebellum, brainstem, vessels, basal forebrain) use ids below 100. Both
% masks below are built from that split instead of from hand-kept id lists,
% so no structure can be forgotten.

% Non-cortical structures are excluded from the thickness simulation and keep
% their original labels. The ventricle/CSF labels and the two generic cerebral
% white matter labels are deliberately not part of this mask. They do not
% describe a structure that has to be preserved, and after the vbdist fill
% above they also cover the deep and periventricular WM, which has to stay WM
% and has to keep seeding the distance map below.
ventricles = [4 11 49:52];
mask_orig = atlas > 0 & atlas < 100 & ~is_in_atlas(atlas, [ventricles 44 45 46]);
mask_orig = cat_vol_morph(mask_orig,'dd',1);

% Region in which no cortical band may be grown. The distance map below has no
% notion of cortex, so without this a GM band is also added around every other
% WM surface, which is anatomically wrong for the corpus callosum and for the
% periventricular WM: the CSF facing them is ventricle, not sulcus.
%
% This has to be an exclusion mask. Allowing the band only inside the cortical
% parcels (atlas>=100) does not work, because the band is grown to a constant
% thickness: wherever the original cortex was thinner than the requested
% thickness, the band has to extend beyond the parcel into the surrounding CSF
% (id 46) and across the WM boundary (ids 44/45). An inclusive cortex mask
% truncates the band there and leaves large CSF spaces instead.
mask_noband = mask_orig | is_in_atlas(atlas, ventricles);

% Soften the border between the simulated cortex and the original labels that
% are kept inside mask_orig. The two label sources can differ considerably
% right at that border (the simulated side may carry a grown GM band where the
% original side has WM), so a narrow transition shows up as a seam.
%
% The weights come from a smoothed binary mask. A linear ramp over the signed
% distance to the border was tried and is clearly worse: it has a constant
% width, but it is only piecewise linear, and the kinks where it reaches 0 and
% 1 are themselves visible as edges (its second derivative there is two orders
% of magnitude larger than anywhere in a Gaussian profile). Smoothing has no
% such kinks, so widen border_fwhm if the transition is still too sharp.
%
% spm_smooth expects the FWHM in voxels for array input, thus the width in mm
% has to be divided by the voxel size.
border_fwhm = 2;  % width of the transition zone in mm

mask_soft = single(mask_orig);
spm_smooth(mask_soft, mask_soft, border_fwhm./vx);
mask_soft = min(max(mask_soft,0),1);

% create mask for mainly occipital and frontal areas (cortical parcels only,
% 31:32 is the amygdala and was a leftover in the frontal list)
region1 = [108:109 114:115 128:129 134:135 144:145 148:149 156:161 170:171 176:177 196:197];
region3 = [102:105 120:121 132:133 136:137 146:147 152:155 162:165 172:173 178:179 190:191 202:205];

mask_thickness{1} = is_in_atlas(atlas, region1); % mainly occipital
mask_thickness{3} = is_in_atlas(atlas, region3); % mainly frontal
mask_thickness{2} = ~mask_thickness{1} & ~mask_thickness{3}; % remaining parts

mask = round(label) > 0;

% force stronger PVE effects by smoothing (FWHM in mm, converted to voxels)
spm_smooth(label,label,1.25./vx);

label1 = cell(numel(simu.thickness),1);

% save original segmentation to later include original cerebellum and basal
% ganglia
Yseg0 = Yseg(:,:,:,1:3);
Yseg(:,:,:,1:3) = 0;

% apply gray closing to strengthen thin WM structures
label = cat_vol_morph(label,'gc',2);

% vary range of PVE from -0.25..0.25 in 15 steps to get more realistic PVE
% effects (optionally weighted)
pve_range = linspace(-0.25,0.25,15);

for pve_step = 1:numel(pve_range)
  % label of this PVE step before the cortical band is built, used below to
  % keep the original tissue inside the excluded structures
  label_step = round(label+pve_range(pve_step));

  % define wm and remove disconnected regions
  wm  = label_step == wm_val;
  wm = cat_vol_morph(wm,'l',1, vx);

  % the excluded structures must not seed the cortical band
  wm(mask_orig) = 0;

  % Euclidean distance to the WM surface, with the usual voxelsize/2
  % correction that puts the surface on the voxel face.
  %
  % cat_bwdist is CAT's separable distance transform. It is used and not
  % bwdist of the Image Processing Toolbox, because this distance defines the
  % cortical band and a toolbox dependent branch would make the simulated
  % thickness depend on the installation rather than only on the parameters.
  % It is also exact, whereas cat_vbdist propagates the vector to the nearest
  % object voxel and accumulates a small error with increasing distance.
  % The voxel size is passed, so the distance is returned in mm directly and
  % anisotropic voxels are handled correctly.
  Dwm = cat_bwdist(single(wm), vx) - 0.5*mean(vx);

  for k=1:numel(simu.thickness)
    label1{k} = label_step;

    label1{k}(~wm) = csf_val;
    label1{k}(~mask) = 0;

    % limit dilated gm to defined thickness, but not into the ventricles or
    % the non-cortical structures (corpus callosum, periventricular WM)
    label1{k}(label1{k} == csf_val & Dwm <= simu.thickness(k) & ~mask_noband) = gm_val;

    % Keep the original tissue inside the excluded structures. The line
    % 'label1{k}(~wm) = csf_val' above turned all of them into CSF, because
    % their WM was removed from the seed and their GM is not WM either. The
    % soft blend at the end would then average a simulated CSF voxel with an
    % original WM or GM voxel, and in a T1w image that average is exactly a
    % GM-like intensity, which appears as a rim around every excluded
    % structure. Widening the blend only widens that rim, so the two label
    % sources have to agree here instead.
    label1{k}(mask_orig) = label_step(mask_orig);
  end

  % replace tissue maps with modified label
  for j = 1:3
    if isscalar(simu.thickness)
      % only simulate 2mm thickness
      tmp_seg = single(round(label1{1}) == (j));
    else
      % the three region masks tile the whole volume, thus every voxel is
      % overwritten below and the map only has to be preallocated here
      tmp_seg = zeros(d, 'single');
      for k = 1:3
        tmp_seg(mask_thickness{k}) = single(round(label1{k}(mask_thickness{k})) == (j));
      end
    end

    Yseg(:,:,:,order(j)) = Yseg(:,:,:,order(j)) + tmp_seg/numel(pve_range);
  end
end

% restore the original tissue fractions inside the excluded structures
for k = 1:3
  Yseg(:,:,:,order(k)) = Yseg(:,:,:,order(k)).*(1-mask_soft) + Yseg0(:,:,:,order(k)).*mask_soft;
end

% Restore the original tissue fractions deep inside the WM.
%
% Everything above works on hard labels, and the grey closing that repairs
% thin WM structures also fills every small hole in the WM. Those holes are
% exactly the local GM fractions of the original segmentation, and in a
% simulation without thickness manipulation they are what makes the WM look
% like tissue: the synthesis turns a WM fraction of 0.95 into an intensity
% between the WM and the GM mean, so the small fluctuations of the original
% segmentation carry over into the image. After the closing and the averaging
% over the PVE steps the WM fraction is exactly 1 everywhere in the interior,
% and the synthesized WM is therefore perfectly flat.
%
% Blending the original fractions back in for the WM interior restores that
% texture, and because the intensity is a function of the fractions alone
% (see synthesize_from_segmentation) the result is identical to the intensity
% of the untouched simulation there.
%
% The restore has to stay away from the GM/WM boundary, otherwise it would
% move the boundary and change the simulated thickness. Two things keep it
% away: the core is eroded by wm_core_dist, and the weight is smoothed with
% the same width as the border blend above, so it has decayed to zero well
% before the boundary is reached. The erosion also excludes the gyral WM
% blades, which are too thin to hold a core at all.
wm_core_dist = 2;  % distance to the WM surface that is left untouched, in mm

Ywm = Yseg(:,:,:,order(3));
wm_soft = single(cat_vol_morph(Ywm > 0.99,'de',wm_core_dist,vx));
spm_smooth(wm_soft, wm_soft, border_fwhm./vx);
wm_soft = min(max(wm_soft,0),1);

% Limit the restored deviation so that the ground truth stays WM, i.e. the
% blended WM fraction cannot drop below 0.5 and the label cannot drop below
% 2.5. Without this, a hypointense lesion of the original segmentation would
% reappear as an island of GM in the middle of the WM of the ground truth,
% which is not what a phantom with a defined cortical thickness should carry.
% The cap is a per-voxel weight and not a clamp of the fractions, so the three
% classes still sum to one. It leaves the ordinary fluctuations untouched,
% they are an order of magnitude smaller than the 0.5 it allows.
wm_soft = wm_soft .* min(1, max(Ywm - 0.5,0) ./ ...
                            max(Ywm - Yseg0(:,:,:,order(3)), eps('single')));
clear Ywm

for k = 1:3
  Yseg(:,:,:,order(k)) = Yseg(:,:,:,order(k)).*(1-wm_soft) + Yseg0(:,:,:,order(k)).*wm_soft;
end

% update ground truth label
label = zeros(d, 'single');

for k = 1:3
  label = label + k*Yseg(:,:,:,order(k));
end


% close remaining holes in CSF
mask = label > 0.5;
label(mask ~= cat_vol_morph(mask,'dc',4)) = 1;


%==========================================================================
% function Yseg = simulate_atrophy(simu, Yseg, dims, template_dir, idef_name, Vref)
%
% Purpose
%   Apply regional atrophy by increasing CSF (and effectively reducing GM)
%   within specified atlas ROIs warped into native space.
%
% Inputs
%   simu         - struct with .atrophy = {atlasName, roiIds[], factors[]}
%   Yseg         - single(dims,3): tissue maps in GM/WM/CSF order, i.e. the
%                  SPM class order, thus Yseg(:,:,:,3) is CSF.
%   dims         - [nx ny nz] dimensions of the volume.
%   template_dir - path to atlas templates; atlas is warped categorically.
%   idef_name    - inverse deformation field to native space.
%   Vref         - struct: target volume definition of the current grid, used
%                  to resample the warped atlas if it does not match dims.
%
% Output
%   Yseg         - updated tissue maps with CSF increased in target ROIs and
%                  all three classes renormalized to a sum of 1 per voxel.
%==========================================================================
function Yseg = simulate_atrophy(simu, Yseg, dims, template_dir, idef_name, Vref)

% warp atlas to native space using categorical interpolation
fprintf('Transform atlas to native space. This may take a while...\n');
atlas_name = fullfile(template_dir,[simu.atrophy{1} '.nii']);
atlas = cat_vol_defs(struct('field1',{{idef_name}},'images',{{atlas_name}},'interp',-1,'modulate',0));
atlas = single(atlas{1}{1});

% The deformation field is always defined for the original image, thus the
% warped atlas has to be resampled if the current grid differs. Nearest
% neighbour keeps the labels categorical.
if any(size(atlas) ~= dims)
  Vdef = spm_vol(idef_name);
  Vdef = Vdef(1);
  atlas_res = zeros(dims, 'single');
  for sl = 1:dims(3)
    M  = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vref.mat\Vdef.mat\M;
    atlas_res(:,:,sl) = spm_slice_vol(atlas, M1, dims(1:2), 0);
  end
  atlas = atlas_res;
end

% go through all defined ROIs
for i = 1:numel(simu.atrophy{2})
  ind_atlas = round(atlas) == simu.atrophy{2}(i);

  % check that the ROI exists, an id that is not in the atlas gives an
  % all-false mask (and not an empty array, which never triggered here)
  if ~any(ind_atlas(:))
    fprintf('ROI #%d does not exist in %s and is skipped.\n', simu.atrophy{2}(i), simu.atrophy{1});
    continue
  end


  % create atrophy mask with defined value in mask and otherwise 1
  mod_atlas = zeros(dims, 'single');
  mod_atlas(ind_atlas) = simu.atrophy{3}(i);
  mod_atlas(~ind_atlas) = 1.0;
  
  [~, maxind] = max(Yseg,[],4);
  GM_atlas = sum(maxind(ind_atlas) == 1);

  % increase CSF by factor defined in atrophy mask
  Yseg(:,:,:,3) = mod_atlas.*Yseg(:,:,:,3);
  % renormalize tissue posteriors so GM+WM+CSF = 1 per voxel, keeping the
  % weighted-sum label bounded in [0,3] (argmax and % reduction unchanged)
  Ysum = Yseg(:,:,:,1) + Yseg(:,:,:,2) + Yseg(:,:,:,3);
  Ysum(Ysum==0) = 1;                       % guard background
  for c = 1:3, Yseg(:,:,:,c) = Yseg(:,:,:,c)./Ysum; end
  
  [~, maxind] = max(Yseg,[],4);
  GM_atlas_simu = sum(maxind(ind_atlas) == 1);
  fprintf('GM reduction in region #%d = %3.2f%s\n',simu.atrophy{2}(i),100*(1-GM_atlas_simu/GM_atlas),'%');
end


%==========================================================================
% function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d, WMH)
%
% Purpose
%   Generate a T1-like image from provided GM/WM/CSF maps. The intensity of a
%   voxel is the probability weighted mixture of the tissue means, where the
%   probabilities are the external segmentation and the remaining probability
%   describes everything that is not brain.
%
% Inputs
%   vol_seg - single(dims,3): tissue probabilities in GM/WM/CSF order (the SPM
%             class order given by res.lkp), summing to at most 1 per voxel.
%   name    - base name for progress display.
%   res     - struct from SPM segmentation. Used are image(1) and Tbias for the
%             bias corrected image. If WMHs are simulated, res.mn carries one
%             additional column with the WMH intensity (see simulate_WMHs).
%   mn      - 3x1 means for GM/WM/CSF.
%   d       - [nx ny nz] dimensions.
%   WMH     - optionally add white matter hyperintensities (WMHs)
%
% Output
%   Ysimu   - synthesized T1-weighted image volume (single).
%
% Algorithm
%   Pbrain = sum of the given tissue probabilities (plus the WMH map). With
%     s = max(Pbrain,1):
%       I = ( sum_k mn(k)*seg_k [+ mn_WMH*WMH] + Ibias*(s-Pbrain) ) / s
%   where Ibias is the bias corrected input image. Inside the brain Pbrain is 1
%   and this is the linear partial volume mixture of the tissue means. Outside
%   the brain Pbrain is 0 and the bias corrected image is kept, so skull and
%   background stay as they are. In between the two blend continuously.
%
% Notes
%   - The SPM mixture model is not evaluated here. Its Gaussian likelihoods for
%     the non-brain classes are unnormalized densities on an arbitrary scale,
%     so mixing them into the same normalization sum as the tissue
%     probabilities darkened the brain boundary by an ill-defined amount.
%   - WMHs are added on top of the WM class, so Pbrain can reach 2 there and
%     the result is the mean of the WM and the WMH intensity.
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

% intensity of the additional WMH class, stored by simulate_WMHs as the last
% entry of res.mn (it is not part of res.lkp)
if ~isempty(WMH)
  intensity_WMH = res.mn(1,K);
end

spm_progress_bar('init',length(x3),['Working on ' name],'Planes completed');
for z = 1:length(x3)

  % Bias corrected image. It provides the intensities of everything that is
  % not GM/WM/CSF, i.e. skull, soft tissue and background.
  f  = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
  bf = exp(transf(chan.B1,chan.B2,chan.B3(z,:),chan.T));
  cr = bf.*f;

  msk = (f==0) | ~isfinite(f);

  % The expected intensity is the probability weighted mixture of the tissue
  % means, where the probabilities are given by the external segmentation.
  % The remaining probability (1-Pbrain) is non-brain and keeps the intensity
  % of the bias corrected image.
  tmp    = zeros(d(1:2));
  Pbrain = zeros(d(1:2));
  for k = 1:3
    seg    = double(vol_seg(:,:,z,k));
    tmp    = tmp + mn(k)*seg;
    Pbrain = Pbrain + seg;
  end
  if ~isempty(WMH)
    seg    = double(WMH(:,:,z));
    tmp    = tmp + intensity_WMH*seg;
    Pbrain = Pbrain + seg;
  end

  % WMHs are added on top of the WM class, thus Pbrain can exceed 1. The
  % normalization then averages the contributions instead of adding a
  % non-brain part, and s>=1 keeps the division safe everywhere.
  s   = max(Pbrain, 1);
  tmp = (tmp + cr.*(s - Pbrain)) ./ s;

  tmp(msk) = 1e-3;

  Ysimu(:,:,z) = tmp;

  spm_progress_bar('set',z);
end

% Sometimes huge values occur due to bias correction in the noisy
% background and we have to limit values to 98% percentile
th = get_percentile(Ysimu,98);
Ysimu(Ysimu>th) = th;
spm_progress_bar('clear');


%==========================================================================
% function th = get_percentile(Y, p)
%
% Purpose
%   Percentile of the finite values of Y. This replaces prctile, which
%   requires the Statistics and Machine Learning Toolbox, to keep the
%   toolbox requirements of mri_simulate limited to SPM and CAT.
%
% Inputs
%   Y - numeric array (of any size); non-finite values are ignored.
%   p - percentile in the range 0..100.
%
% Output
%   th - the p-th percentile of Y, using the same convention as prctile
%        (sorted value i represents the 100*(i-0.5)/n percentile, with
%        linear interpolation in between and clamping at both ends).
%==========================================================================
function th = get_percentile(Y, p)

Y = sort(Y(isfinite(Y)));
n = numel(Y);

if n == 0
  th = 0;
  return
end

pos = p/100*n + 0.5;
if pos <= 1
  th = double(Y(1));
elseif pos >= n
  th = double(Y(n));
else
  lo = floor(pos);
  w  = pos - lo;
  th = (1-w)*double(Y(lo)) + w*double(Y(lo+1));
end


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

% invert field for neg. values by changing the sign of the modulation, which
% swaps bright and dark areas but keeps the defined amplitude (the previous
% reciprocal 1./rf_field resulted in Inf for zero-valued voxels and in an
% amplitude that was no longer related to rf.percent)
if rf.percent < 0
  rf_field = -rf_field;
end

ind = isfinite(rf_field);
rf_field = 1 + rf_field - mean(rf_field(ind));
Ysimu(ind) = rf_field(ind).*Ysimu(ind);


%==========================================================================
% function [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf, vx)
%
% Purpose
%   Create a smooth, random RF bias field via FFT-domain filtering with a
%   strength-dependent grid size, interpolate to image size, scale by rf.percent,
%   and apply to the simulated image.
%
% Inputs
%   Ysimu  - image to modulate.
%   rf     - struct: .type = [strength, rngSeed], .percent amplitude (signed).
%   vx     - [vx vy vz]: voxel size in mm of Ysimu, used to define the
%            smoothness of the field in mm rather than in voxels.
%
% Outputs
%   Ysimu    - modulated image.
%   rf_field - generated RF field after smoothing and scaling.
%
% Notes
%   - The field is smooth by construction (it originates from an NxNxN random
%     field), and is therefore built and smoothed on a coarse grid that is
%     interpolated to the image size afterwards. Smoothing at full resolution
%     would need a border padding of 3*FWHM voxels in every direction, which
%     for a 0.5mm image means a temporary volume of more than 1GB without
%     changing the resulting field noticeably.
%   - The smoothing FWHM is defined in mm, so that the same field is obtained
%     for any input resolution. The value corresponds to the 30 voxels that
%     were previously used for the 0.5mm reference phantom.
%==========================================================================
function [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf, vx)

dim = size(Ysimu);

% set seeds to defined value
rng(rf.type(2),'twister')

N = 2^round(rf.type(1)); % Define the size of the 3D field w.r.t. defined strength
fwhm_mm = 15;            % Smoothing size in mm

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

% Shift back and apply inverse 3D FFT. The filter is symmetric, thus the
% result is real apart from rounding errors, but ifftn only drops the
% imaginary part for an exactly conjugate symmetric input. Taking the real
% part explicitly avoids that min/max below compare complex numbers by their
% magnitude, which would silently return a different field.
filteredField = real(ifftn(ifftshift(filteredFFT)));

% Coarse grid on which the field is smoothed: at least 6 samples per FWHM are
% used, which is more than enough for a Gaussian and keeps the padded volume
% small. The field is never sampled finer than the image itself.
res_c = max(vx, fwhm_mm/6);              % coarse voxel size in mm
dim_c = max(4, ceil(dim.*vx./res_c));    % coarse dimensions
res_c = dim.*vx./dim_c;                  % exact voxel size after rounding
fwhm_c = fwhm_mm./res_c;                 % FWHM in coarse voxels
pad = ceil(3*max(fwhm_c));               % pad border to prevent smoothing issues

% Original grid
[x, y, z] = ndgrid(1:N, 1:N, 1:N);

% Coarse grid dimensions
[xq, yq, zq] = ndgrid(linspace(1, N, dim_c(1)), ...
                      linspace(1, N, dim_c(2)), ...
                      linspace(1, N, dim_c(3)));

% Interpolate using interpn
filteredField = interpn(x, y, z, filteredField, xq, yq, zq, 'linear');

% resize field and scale it to 0..1
filteredField = filteredField - min(filteredField(:));
filteredField = filteredField/max(filteredField(:));

% pad border and create background of 0.5 to prevent smoothing issues at image border
vol = 0.5*ones(dim_c(1)+2*pad,dim_c(2)+2*pad,dim_c(3)+2*pad);
vol(pad+1:dim_c(1)+pad,pad+1:dim_c(2)+pad,pad+1:dim_c(3)+pad) = filteredField;

% smooth and remove padding
spm_smooth(vol,vol,fwhm_c);
vol = vol(pad+1:dim_c(1)+pad,pad+1:dim_c(2)+pad,pad+1:dim_c(3)+pad);

% interpolate the smoothed coarse field to the image size
[xc, yc, zc] = ndgrid(1:dim_c(1), 1:dim_c(2), 1:dim_c(3));
[xq, yq, zq] = ndgrid(linspace(1, dim_c(1), dim(1)), ...
                      linspace(1, dim_c(2), dim(2)), ...
                      linspace(1, dim_c(3), dim(3)));
rf_field = interpn(xc, yc, zc, vol, xq, yq, zq, 'linear');

% scale to an amplitude and mean of 1
rf_field = rf_field - min(rf_field(:));
rf_field = rf_field/max(rf_field(:));

% apply defined percent
rf_field = abs(rf.percent)/100 * rf_field;

% invert field for neg. values by changing the sign of the modulation, which
% swaps bright and dark areas but keeps the defined amplitude (the previous
% reciprocal 1./rf_field divided by zero here, because the field was scaled
% to a range of exactly 0..1 above, and returned a NaN/-Inf field)
if rf.percent < 0
  rf_field = -rf_field;
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
%   1) Warp the MICCAI2017 WMH prior map to native space, resample it to the
%      current grid if necessary, and lightly smooth it.
%   2) Create a random 3D field, resample it to image dimensions, threshold to
%      form spatial support for patchy WMH distribution.
%   3) Erode WM to ensure spacing from GM and constrain WMHs to deep WM.
%   4) Combine WMH prior^(1/(strength-0.8)) with random field and WM mask,
%      smooth and normalize to [0,1].
%   5) Update label_pve by adding 5*WMH/strength^0.75 and clipping to max class
%      label 4, thereby introducing an additional WMH class contribution.
%   6) Extend res.mn by a WMH class intensity equal to the (weighted) GM mean.
%
% Notes
%   - Class encoding after this step: CSF=1, GM=2, WM=3, WMH=4.
%   - The strength parameter shapes the prior map nonlinearly; larger values
%     emphasize regions with higher WMH probability.
%   - The WM erosion step helps avoid spuriously labeling GM/CSF boundaries
%     as WMH.
%   - Only res.mn is extended by the WMH class, because mg/lkp keep describing
%     the Gaussians of the original SPM segmentation and are still used to
%     compute weighted class means. synthesize_from_segmentation addresses the
%     WMH class by its index K = numel(res.mg)+1.
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

% The deformation field is always defined for the original image, thus the
% warped map has to be resampled if the current grid differs (which is the
% case after the internal resampling to 0.5mm for thickness simulation).
Vref = res.image(1);
dim  = Vref.dim(1:3);
if any(size(WMH) ~= dim)
  Vdef = spm_vol(idef_name);
  Vdef = Vdef(1);
  WMH_res = zeros(dim, 'single');
  for sl = 1:dim(3)
    M  = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vref.mat\Vdef.mat\M;
    WMH_res(:,:,sl) = spm_slice_vol(WMH, M1, dim(1:2), 1);
  end
  WMH = WMH_res;
end

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

% Use mean of GM for the additional class. The mixing proportions mg sum to 1
% within each class, thus the weighted mean is used here in the same way as
% the GM/WM/CSF means are estimated in the main function.
% Only res.mn is extended (and not mg/lkp), because those keep describing the
% Gaussians of the original SPM segmentation. The WMH class is addressed by
% its index K in synthesize_from_segmentation.
ind_GM = res.lkp == 1;
mg_GM  = res.mg(ind_GM);
mn_GM  = res.mn(1,ind_GM);
intensity_WMH = sum(mg_GM(:) .* mn_GM(:));
K = numel(res.mg) + 1;
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
%   Ycls - cell of uint8(dims): Tissue class posteriors (0..255) from
%          cat_spm_preproc_write8, one per tissue of the TPM and in SPM/CAT
%          convention, i.e. GM, WM, CSF, head tissues, background.
%   res  - struct: SPM/CAT segmentation result with fields used here:
%              .mn  (means per class/component)
%              .mg  (mixture weights)
%              .lkp (lookup of class index per component, 1..6)
%   dim  - [nx ny nz]: Volume dimensions of the images.
%   T3th - 1x3 double: Intensity anchors [CMth, GMth, WMth] as returned by
%          get_tissue_thresholds and passed on to APRG in the position where
%          it expects [CSF GM WM]. CMth is an extrapolated low anchor below
%          the GM mean, GMth is the GM mean and WMth a robust WM anchor.
%
% Output
%   Yb   - logical/single(dims): Brain mask estimated by APRG (1=brain, 0=non-brain).
%
% Algorithm (summary)
%   1) Build a 4D stack P(:,:,:,1..6) from Ycls (uint8 posteriors), as required
%      by CAT12's APRG routine.
%   2) Call cat_main_APRG(Ysrc, P, res, T3th) to obtain the brain mask. The
%      intensity anchors themselves are estimated beforehand by the caller in
%      get_tissue_thresholds.
%
% Notes
%   - Assumes T1-like contrast (WM > GM > CSF) for threshold reasoning; a median
%     WM intensity is used to be robust to WMHs.
%   - Requires CAT on the MATLAB path; uses cat_main_APRG.
%   - The output mask matches the input volume dimensions and orientation.
%==========================================================================
function Yb = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th)

% We need a 4D array for cat_main_APRG. All classes of the TPM are passed and
% not only the first six of the SPM default one: cat_main_APRG uses the last
% class as background and everything from the fourth on as head tissue, so
% dropping the last classes of a TPM with more tissues (the Blaiotta TPM has
% seven) would hand it soft tissue in the place of the background.
P = zeros([dim numel(Ycls)],'uint8');
for i=1:numel(Ycls)
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
%   Yseg   - single(dims,3): current class volumes in GM/WM/CSF order, i.e.
%            the SPM class order, derived from LAS-corrected intensities.
%   Ysrc   - single(dims): source image passed on to cat_main_updateSPM1639.
%   Ycorr  - single(dims): LAS-corrected image used by cat_vol_partvol.
%   Ycls   - cell of uint8(dims): SPM tissue posteriors (0..255), one per
%            tissue of the TPM.
%   Yy     - deformation field to the TPM, required by the CAT functions.
%   res    - segmentation structure from SPM segmentation.
%   vx_vol - voxel size in mm
%
% Output
%   Yseg   - class volumes with WMHs reassigned to WM, leaving cortex and
%            boundaries untouched as much as possible.
%
% Note
%   This step costs minutes because cat_vol_partvol solves several region
%   growings with cat_vol_laplace3R, and for an image without WMHs all of it
%   ends in an empty mask that leaves Yseg unchanged. A cheap geometric
%   pre-test to skip it was tried and rejected: WMHs are hypointense islands
%   inside WM, but so are the thalamus and the basal ganglia, and the two
%   cannot be told apart by size or intensity. Filling the WM mask with a
%   closing does separate them, yet only for lesions smaller than the closing
%   radius, so lesions above roughly 7mm (and any periventricular rim) become
%   invisible while deep grey matter stays excluded. Separating the two needs
%   the atlas information that cat_vol_partvol already uses. Set
%   simu.closeWMHholes = 0 to skip this step for inputs known to be free of
%   WMHs, such as the Colin27 template.

% we have to prepare some parameters for cat_main_updateSPM1639
global cat; cat_defaults;
job = cat;
job.extopts.inv_weighting = 0; job.extopts.verb = 0;

% Internal working resolution of cat_vol_partvol. It defaults to 0.7mm, which
% for a 0.5mm phantom means that the iterative region growing in
% cat_vol_laplace3R (the dominating cost of this step, marked as bottleneck in
% cat_vol_partvol itself) runs on ~12 million voxels and takes minutes.
% Only Yl1 and the WMH class are used below, both as coarse masks that are
% afterwards dilated by 8mm and closed by 12mm, so a finer grid brings nothing
% here. cat_vol_partvol restores the native resolution of its outputs, thus
% this only affects the internal computation. Lowering it towards 0.7 gives a
% more detailed WMH mask at a steeply increasing cost (roughly linear in the
% number of voxels, i.e. ~3x from 1.0 to 0.7 and ~12x from 1.5 to 0.7).
job.extopts.uhrlim = max(1.0, max(vx_vol));
P = zeros([size(Ycls{1}) numel(Ycls)],'uint8');
for i=1:numel(Ycls), P(:,:,:,i) = Ycls{i}; end
clear Ycls;
tpm.dat = cell(numel(res.tpm),1);
tpm.V = res.tpm;
for i=1:numel(res.tpm)
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

% Mask of structures where WMHs must not be corrected. The dilation by 8mm and
% the closing by 12mm are deliberately coarse, but at native resolution they
% are distance transforms over the full volume (~20s for a 0.5mm image). The
% mask is therefore built on a 1.5mm grid, where the same radii cost a
% fraction of that, and is then resampled back.
Ynwmh = NS(Yl1,LAB.TH) | NS(Yl1,LAB.BG) | NS(Yl1,LAB.HC) | NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS);
Yvt   = NS(Yl1,LAB.VT);

[Ynwmh_r, Yvt_r, resTr] = cat_vol_resize({single(Ynwmh), single(Yvt)}, ...
                                         'reduceV', vx_vol, 1.5, 32);
vx_r = resTr.vx_volr;
Ynwmh_r = cat_vol_morph(cat_vol_morph( Ynwmh_r>0.5, 'dd', 8, vx_r), 'dc', 12, vx_r) & ...
          ~cat_vol_morph( Yvt_r>0.5, 'dd', 4, vx_r);
Ynwmh = cat_vol_resize(single(Ynwmh_r), 'dereduceV', resTr) > 0.5;

% The WMH class of cat_vol_partvol is only present when WMH correction ran.
%
% It is written to Ycls{7}, which for a TPM with seven tissues (the Blaiotta
% TPM) is also the index of the background class of the segmentation. Both are
% told apart by the brain mask, because WMHs lie inside the brain and the
% background does not.
if numel(Ycls) >= 7 && ~isempty(Ycls{7})
  Ywmh = Ycls{7}>0 & ~Ynwmh & Yb>0.5;
else
  Ywmh = false(size(Ynwmh));
end

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



function T3th = get_tissue_thresholds(Ysrc, Ycls, mn)
% GET_TISSUE_THRESHOLDS - Estimate robust CSF/GM/WM thresholds
%
% Purpose
%   Compute the three intensity anchors that LAS and APRG expect in the
%   [CSF GM WM] positions:
%     - CMth: the low anchor. It is not the CSF mean but the GM mean mirrored
%       away from WM (2*GMmean - WMth), clipped so that it never exceeds the
%       CSF mean. This puts it below the GM mean and adapts to the current
%       image contrast. For inverted contrast (CSF > WM) the CSF mean is used.
%     - GMth: the GM intensity of the mixture model (mn(1)).
%     - WMth: a robust white-matter anchor derived from the WM intensity of
%       the mixture model and the median intensity within high WM posterior
%       voxels to mitigate WMH effects.
%
% Inputs
%   Ysrc - single(dims): Source image used to compute medians within WM.
%   Ycls - cell of uint8(dims): SPM/CAT posteriors (0..255), one per tissue of
%          the TPM, where Ycls{2} corresponds to WM.
%   mn   - 3x1 double: intensities of the tissues in SPM class order, i.e.
%          [GM WM CSF], as estimated by the caller from the Gaussian mixture
%          of the segmentation. They are passed and not derived here, because
%          the CSF entry is not the class mean but the darkest CSF Gaussian
%          and both users of the intensities have to see the same value.
%
% Output
%   T3th - 1x3 double: [CMth, GMth, WMth] anchors used by skull stripping
%          and LAS normalization routines.
%
% Algorithm
%   1) WMth: max(mn(WM), median(Ysrc within high WM posterior)).
%   2) CMth: if intensities are inverted (CSF > WM), use the CSF anchor;
%      otherwise use 2*GMmean - WMth, clipped not to exceed the CSF anchor.
%   3) Return [CMth, GMmean, WMth].
%
% Notes
%   - Class indices follow SPM/CAT convention (1=GM, 2=WM, 3=CSF), so mn(1)
%     is GM and mn(3) is CSF.
%   - High WM posterior is defined with a conservative threshold (Ycls{2}>192
%     on 0..255 scale) to obtain a stable median.
%   - These anchors are designed for T1-like contrast but include a simple
%     inversion check (CSF > WM) for robustness.

% Use median for WM threshold estimation to avoid problems in case of WMHs
WMth = double(max(mn(2), cat_stat_nanmedian(Ysrc(Ycls{2}>192))));
if mn(3)>mn(2) % invers
  CMth = mn(3);
else
  CMth = min([mn(1) - diff([mn(1),WMth]), mn(3)]);
end

T3th = double([CMth, mn(1), WMth]);