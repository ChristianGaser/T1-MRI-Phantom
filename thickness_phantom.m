function [Yp0, Ysimu, info] = thickness_phantom(opt)
% THICKNESS_PHANTOM - Folded sphere phantom with exactly known GM thickness
%
% Overview:
%   `thickness_phantom` creates a geometric phantom to validate cortical
%   thickness measures. The white matter is a sphere with regular folds, the
%   grey matter is the band of constant thickness around it, and a CSF layer
%   surrounds the grey matter. The result is written as a partial volume
%   label image (CSF=1, GM=2, WM=3, background=0) and as an ideal T1w image,
%   both with a JSON sidecar.
%
%   The phantom is free of artefacts by design and it is fully
%   deterministic: there is no noise, no bias field and no random component
%   anywhere, so that the error a thickness measure shows on it is its own
%   theoretical error and nothing else. Use mri_simulate for a simulation
%   that carries noise, bias fields and a realistic anatomy.
%
%   Why the thickness is exact:
%   The phantom is not built by offsetting the radius of the folded sphere.
%   A radial offset of t gives a normal thickness of t*cos(alpha), where
%   alpha is the angle between the radius and the surface normal, so on the
%   flanks of the folds the true thickness would be smaller than requested.
%   Instead the euclidean distance D to the white matter is computed once
%   with cat_bwdist, and both cortical boundaries are level sets of that one
%   distance map: the WM/GM boundary is D = wm_offset and the GM/CSF boundary
%   is D = wm_offset + thickness. Because D is 1-Lipschitz and its gradient
%   has unit length, two of its level sets are exactly thickness apart
%   everywhere, whatever the shape of the folds. This is the same
%   construction that mri_simulate uses to give a real brain a constant
%   cortical thickness, only here it starts from an analytic surface instead
%   of a segmentation.
%
%   The distance map is computed on a grid that is supersampled by
%   `supersample`, so the discretization of the folded surface stays well
%   below the voxel size of the phantom.
%
%   Partial volume:
%   The PVE follows the idea of mri_simulate, i.e. the tissue boundary is
%   shifted across a set of sub-voxel offsets, each offset yields a hard
%   label image and the results are averaged. Here the offsets are applied to
%   the distance map directly, which shifts both cortical boundaries along
%   their normal by a known sub-voxel amount. The offsets are the midpoints
%   of `pve_steps` equal intervals covering `pve_range` voxels, so for a
%   locally flat boundary the average reproduces the exact linear partial
%   volume ramp of one voxel width.
%
%   Requirements: SPM12 or SPM25 with CAT >= 26 installed. Distances use
%   CAT's cat_bwdist so that the result does not depend on the installed
%   MATLAB toolboxes. No MATLAB toolboxes beyond base MATLAB are needed.
%
% Syntax:
%   [Yp0, Ysimu, info] = thickness_phantom(opt)
%
% Parameters:
%   opt (struct): Phantom parameters. Defaults are applied for missing fields.
%       - 'dim' ([nx ny nz]): Dimensions of the phantom in voxels.
%         Default: [128 128 128].
%       - 'vx' (double): Isotropic voxel size in mm. Default: 0.5.
%       - 'radius' (double): Mean radius of the WM surface in mm. Default: 22.
%       - 'amplitude' (double): Amplitude of the folds in mm, i.e. half of
%         their peak to peak height. Default: 2.5.
%       - 'fold' (char): Folding pattern.
%           'cartesian' - amplitude * cos(kx)*cos(ky)*cos(kz) with
%             k = 2*pi/wavelength. The feature size is the same everywhere,
%             which is what a phantom for a thickness measure should have.
%           'spherical' - amplitude * sin(n1*theta)*cos(n2*phi), the pattern
%             of internal/thickness_phantom.py of T1Prep. Note that its folds
%             get arbitrarily fine towards the poles, where the circumference
%             of a parallel goes to zero while the number of folds along it
%             stays n2, so a part of the phantom is not resolved by any voxel
%             size. Use it to reproduce that phantom, not to measure bias.
%         Default: 'cartesian'.
%       - 'wavelength' (double): Fold wavelength in mm for 'cartesian'.
%         Default: 12.
%       - 'folds' ([n1 n2]): Number of folds in theta and phi for
%         'spherical'. Default: [6 6].
%       - 'thickness' (double or vector): GM thickness in mm. A vector
%         creates one phantom per value; they share the distance map, which
%         is by far the most expensive step. Default: 2.5.
%       - 'csf' (double): Thickness of the CSF layer around the GM in mm.
%         Default: 2.
%       - 'wm_offset' (double): Level of the distance map that defines the
%         WM/GM boundary, in mm. It must stay clearly above the sampling of
%         the distance map, because the level sets of a distance map close to
%         its object are scalloped at the scale of that sampling. The mean
%         radius of the WM surface is corrected for it, so `radius` stays the
%         radius of the WM surface. Default: NaN, i.e. half a voxel.
%       - 'supersample' (integer): Factor by which the grid is supersampled
%         for the distance map. Even values are increased by one, since an
%         odd factor puts a sample of the fine grid exactly on every voxel
%         centre of the phantom. Memory and time grow with its third power.
%         Default: 3.
%       - 'pve_steps' (integer): Number of sub-voxel offsets that are
%         averaged for the partial volume effect. Default: 15.
%       - 'pve_range' (double): Total width of those offsets in voxels. The
%         default of 1 gives the exact linear partial volume ramp of a flat
%         boundary; mri_simulate uses 0.5 because it applies the offsets to a
%         label map that was smoothed beforehand. Default: 1.
%       - 'intensities' ([csf gm wm]): Tissue intensities of the synthesized
%         T1w image, background is 0. Default: [30 100 150].
%       - 'name' (char): Basename of the output files.
%         Default: 'thicknessPhantom'.
%       - 'outdir' (char): Output folder. Default: pwd.
%       - 'write' (logical): Write the images and the JSON sidecars.
%         Default: 1.
%       - 'verbose' (logical): Print progress. Default: 1.
%
% Outputs:
%   Yp0   - single(dim): partial volume label image, CSF=1, GM=2, WM=3 and 0
%           outside. A cell array with one entry per thickness if several
%           thickness values are given.
%   Ysimu - single(dim): synthesized T1w image, same convention for several
%           thickness values.
%   info  - struct (array): geometry, tissue volumes and output filenames.
%           info.volume_ref holds the tissue volumes of the supersampled
%           grid, which is the reference the partial volume label should
%           reproduce (see thickness_phantom_eval).
%
%   Files (for write=1), following the naming of mri_simulate:
%     <name>_res-<vx>mm_desc-Thickness<t>mm_dseg.nii  label image
%     <name>_res-<vx>mm_desc-Thickness<t>mm_T1w.nii   ideal T1w image
%     plus a JSON sidecar next to each of them
%
% Examples:
%   Example 1 - Default phantom with 2.5mm thickness at 0.5mm voxel size,
%               and its evaluation:
%       [~, ~, info] = thickness_phantom;
%       thickness_phantom_eval(info);
%
%   Example 2 - A series of thickness values in one call. They share the
%               distance map, thus this is barely slower than a single one:
%       opt = struct('thickness', 1.5:0.5:3.5);
%       thickness_phantom(opt);
%
%   Example 3 - The folding pattern of T1Prep/internal/thickness_phantom.py:
%       opt = struct('fold','spherical', 'folds',[6 6], 'radius',25, ...
%                    'amplitude',2.5, 'thickness',3, 'csf',1);
%       thickness_phantom(opt);
%
%   Example 4 - Phantom on a 1mm grid, i.e. the sampling of a normal scan:
%       opt = struct('vx',1, 'dim',[64 64 64]);
%       thickness_phantom(opt);
%
% See also: thickness_phantom_eval, mri_simulate

% named tool_version and not version to not shadow the MATLAB builtin version()
tool_version = '0.1.0';

if ~exist('cat_bwdist','file')
  error('CAT >= 26 has to be in the MATLAB path to use thickness_phantom.')
end

def.dim         = [128 128 128];
def.vx          = 0.5;
def.radius      = 22;
def.amplitude   = 2.5;
def.fold        = 'cartesian';
def.wavelength  = 12;
def.folds       = [6 6];
def.thickness   = 2.5;
def.csf         = 2;
def.wm_offset   = NaN;
def.supersample = 3;
def.pve_steps   = 15;
def.pve_range   = 1;
def.intensities = [30 100 150];
def.name        = 'thicknessPhantom';
def.outdir      = pwd;
def.write       = 1;
def.verbose     = 1;

if nargin < 1, opt = def;
else, opt = cat_io_checkinopt(opt, def); end

dim = round(opt.dim(:)');
if isscalar(dim), dim = dim*[1 1 1]; end
vx  = opt.vx;
thickness = opt.thickness(:)';

% An odd supersampling factor puts one sample of the fine grid exactly on
% every voxel centre of the phantom, so the distance map can be read out
% without any interpolation.
ss = max(1, round(opt.supersample));
if mod(ss,2) == 0
  ss = ss + 1;
  fprintf('Supersampling has to be odd and was increased to %d.\n', ss);
end
vxf  = vx/ss;
dimf = dim*ss;

if any(thickness <= 0)
  error('opt.thickness has to be positive.')
end

% Level of the distance map that carries the WM/GM boundary. Half a voxel is
% far enough from the object for the scalloping of the level sets (see the
% header) to stay below a tenth of a voxel.
if ~isfinite(opt.wm_offset), wm_offset = 0.5*vx; else, wm_offset = opt.wm_offset; end
opt.wm_offset = wm_offset;   % so that the sidecars record the resolved value
if wm_offset <= vxf
  fprintf(['Warning: opt.wm_offset of %g mm is not larger than the sampling of the ' ...
           'distance map (%g mm),\n  so the WM surface will be visibly scalloped. ' ...
           'Increase wm_offset or supersample.\n'], wm_offset, vxf);
end

% The whole phantom has to fit into the field of view with a little margin,
% otherwise the CSF layer is cut open at the border and the label image no
% longer describes the geometry that is documented here.
r_max = opt.radius + opt.amplitude + max(thickness) + opt.csf;
if r_max + vx > min(dim)*vx/2
  error(['The phantom needs a radius of %g mm but the field of view only offers ' ...
         '%g mm. Increase opt.dim or reduce the geometry.'], r_max, min(dim)*vx/2);
end

% Mean radius of the analytic surface. The object of the distance map is the
% set of voxels inside that surface, thus its level set at wm_offset lies
% wm_offset outside of it, and the voxelisation moves it back in by half a
% sample on average. Both are corrected here so that opt.radius stays the
% mean radius of the WM surface of the phantom.
R0 = opt.radius - wm_offset + 0.5*vxf;

%--------------------------------------------------------------------------
% Distance to the folded sphere on the supersampled grid
%--------------------------------------------------------------------------
if opt.verbose
  fprintf('Build folded sphere on a %dx%dx%d grid (%.3f mm).\n', dimf, vxf);
end

% Coordinates of the fine grid in mm, centred on the volume. The coarse grid
% below uses the same origin, thus every coarse voxel centre coincides with
% the fine sample (i-1)*ss + (ss+1)/2.
cf = (dimf+1)/2;
xf = ((1:dimf(1)) - cf(1))*vxf;
yf = ((1:dimf(2)) - cf(2))*vxf;
zf = ((1:dimf(3)) - cf(3))*vxf;

% The in plane terms do not depend on the slice and are computed once, the
% volume itself is filled slice by slice. A full set of coordinate volumes on
% the fine grid would need several times the memory of the distance map.
[Xg, Yg] = ndgrid(single(xf), single(yf));
Rxy2 = Xg.^2 + Yg.^2;
switch lower(opt.fold)
  case 'cartesian'
    kw = 2*pi/opt.wavelength;
    CxCy = cos(kw*Xg).*cos(kw*Yg);
  case 'spherical'
    Phi = atan2(Yg, Xg);
  otherwise
    error('Unknown opt.fold ''%s'', use ''cartesian'' or ''spherical''.', opt.fold);
end
clear Xg Yg

Ywm = false(dimf);
for k = 1:dimf(3)
  z = single(zf(k));
  r = sqrt(Rxy2 + z^2);
  switch lower(opt.fold)
    case 'cartesian'
      Rs = R0 + opt.amplitude * CxCy * cos(kw*z);
    case 'spherical'
      % theta is undefined at the centre, where r is 0. Any value does, the
      % centre is deep inside the object.
      th = acos(z ./ max(r, eps('single')));
      Rs = R0 + opt.amplitude * sin(opt.folds(1)*th) .* cos(opt.folds(2)*Phi);
  end
  Ywm(:,:,k) = r <= Rs;
end
clear Rxy2 CxCy Phi r th Rs

if ~any(Ywm(:))
  error('The folded sphere is empty, please check opt.radius and opt.amplitude.')
end

if opt.verbose, fprintf('Compute euclidean distance map.\n'); end

% cat_bwdist is CAT's separable distance transform. It is exact and returns
% the distance in mm when the voxel size is given, and it is used instead of
% bwdist so that the geometry of the phantom does not depend on which MATLAB
% toolboxes are installed.
D = cat_bwdist(single(Ywm), vxf*[1 1 1]);
clear Ywm

% Tissue volumes of the fine grid. They are the reference for the partial
% volume label below, which has to reproduce them although it samples the
% same geometry at a much coarser voxel size.
vol_ref = zeros(numel(thickness),3);
for i = 1:numel(thickness)
  vol_ref(i,3) = nnz(D <= wm_offset);
  vol_ref(i,2) = nnz(D >  wm_offset & D <= wm_offset + thickness(i));
  vol_ref(i,1) = nnz(D >  wm_offset + thickness(i) & ...
                     D <= wm_offset + thickness(i) + opt.csf);
end
vol_ref = vol_ref * vxf^3 / 1000;   % ml

% Read the distance map on the grid of the phantom. No interpolation is
% needed, the fine samples below sit exactly on the coarse voxel centres.
idx = @(n) ((1:n)-1)*ss + (ss+1)/2;
Dc = D(idx(dim(1)), idx(dim(2)), idx(dim(3)));
clear D

%--------------------------------------------------------------------------
% Partial volume label and T1w synthesis, one per thickness value
%--------------------------------------------------------------------------

% Sub-voxel offsets of the boundary. They are the midpoints of pve_steps
% equal intervals and not linspace(-r/2, r/2, n): the midpoint rule
% integrates the uniform distribution of the boundary position exactly,
% while linspace weights the two extreme offsets like all the others and
% biases the resulting ramp.
delta = opt.pve_range*vx*( ((1:opt.pve_steps)-0.5)/opt.pve_steps - 0.5 );

if ~exist(opt.outdir,'dir') && opt.write, mkdir(opt.outdir); end

Yp0   = cell(1,numel(thickness));
Ysimu = cell(1,numel(thickness));
info  = struct('thickness',{},'vx',{},'dim',{},'radius',{},'fold',{}, ...
               'wm_offset',{},'volume_ref',{},'volume_pve',{},'files',{});

for i = 1:numel(thickness)
  t = thickness(i);
  if opt.verbose
    fprintf('Thickness %.2f mm: average %d sub-voxel offsets.\n', t, opt.pve_steps);
  end

  Fcsf = zeros(dim,'single'); Fgm = zeros(dim,'single'); Fwm = zeros(dim,'single');
  for j = 1:opt.pve_steps
    % Subtracting the offset from the distance moves all three boundaries
    % outwards by delta(j) along their normal, i.e. it is the sub-voxel shift
    % of the tissue boundary that the averaging below turns into partial
    % volume. The three hard labels are disjoint by construction.
    d = Dc - delta(j);
    Fwm  = Fwm  + single(d <= wm_offset);
    Fgm  = Fgm  + single(d >  wm_offset     & d <= wm_offset + t);
    Fcsf = Fcsf + single(d >  wm_offset + t & d <= wm_offset + t + opt.csf);
  end
  Fcsf = Fcsf/opt.pve_steps; Fgm = Fgm/opt.pve_steps; Fwm = Fwm/opt.pve_steps;

  Yp0{i} = Fcsf + 2*Fgm + 3*Fwm;

  % T1w synthesis: the intensity is the mixture of the tissue intensities
  % weighted by the partial volume fractions, background contributes 0.
  Ysimu{i} = opt.intensities(1)*Fcsf + opt.intensities(2)*Fgm + opt.intensities(3)*Fwm;

  vol_pve = [sum(Fcsf(:)) sum(Fgm(:)) sum(Fwm(:))] * vx^3 / 1000;   % ml
  clear Fcsf Fgm Fwm

  info(i).thickness  = t;
  info(i).vx         = vx;
  info(i).dim        = dim;
  info(i).radius     = opt.radius;
  info(i).fold       = opt.fold;
  info(i).wm_offset  = wm_offset;
  info(i).volume_ref = struct('csf',vol_ref(i,1),'gm',vol_ref(i,2),'wm',vol_ref(i,3));
  info(i).volume_pve = struct('csf',vol_pve(1), 'gm',vol_pve(2), 'wm',vol_pve(3));
  info(i).files      = struct('label','','t1w','');

  if opt.write
    info(i).files = write_phantom(Yp0{i}, Ysimu{i}, opt, t, tool_version, dim, vx);
  end
end

if isscalar(thickness), Yp0 = Yp0{1}; Ysimu = Ysimu{1}; end

if opt.verbose
  fprintf('================================================================================\n');
end


%==========================================================================
% function files = write_phantom(Yp0, Ysimu, opt, t, tool_version, dim, vx)
%
% Purpose
%   Write the label image, the T1w image and their JSON sidecars. The names
%   follow the BIDS filename grammar in the same way as the output of
%   mri_simulate, i.e. entity-value pairs followed by a suffix, with the
%   voxel size in the res entity and all options collected in one
%   alphanumeric desc label.
%
% Inputs
%   Yp0          - single(dim): partial volume label image.
%   Ysimu        - single(dim): synthesized T1w image.
%   opt          - struct: the parameters of the phantom.
%   t            - double: thickness of this phantom in mm.
%   tool_version - char: version written into the sidecars.
%   dim          - [nx ny nz] dimensions of the volume.
%   vx           - double: isotropic voxel size in mm.
%
% Output
%   files - struct with the filenames of the label and the T1w image.
%==========================================================================
function files = write_phantom(Yp0, Ysimu, opt, t, tool_version, dim, vx)

% The world origin is placed on the centre of the volume, which is also the
% centre of the sphere.
mat = [vx 0 0 -vx*(dim(1)+1)/2; ...
       0 vx 0 -vx*(dim(2)+1)/2; ...
       0 0 vx -vx*(dim(3)+1)/2; ...
       0 0 0  1];

ent_desc = ['_desc-' bids_label(sprintf('Thickness%02dmm', round(10*t)))];
ent_res  = ['_res-' bids_label(sprintf('%g',vx)) 'mm'];
base     = [opt.name ent_res ent_desc];

files.label = fullfile(opt.outdir, [base '_dseg.nii']);
files.t1w   = fullfile(opt.outdir, [base '_T1w.nii']);

V = struct('fname','', 'dim',dim, 'mat',mat, 'pinfo',[1 0 352]', ...
           'dt',[spm_type('float32') 0], 'descrip','thickness_phantom');

fprintf('Save %s\n', files.label);
V.fname = files.label;
spm_write_vol(V, single(Yp0));

fprintf('Save %s\n', files.t1w);
V.fname = files.t1w;
spm_write_vol(V, single(Ysimu));

% The sidecars document the geometry, since the whole point of the phantom
% is that its thickness is known, and the evaluation reads them back.
gen = struct('Name','thickness_phantom', 'Version',tool_version);

par = struct();
par.Thickness   = t;
par.CSFLayer    = opt.csf;
par.Radius      = opt.radius;
par.Amplitude   = opt.amplitude;
par.Fold        = opt.fold;
if strcmpi(opt.fold,'cartesian')
  par.FoldWavelength = opt.wavelength;
else
  par.FoldNumber = opt.folds(:)';
end
par.WMOffset    = opt.wm_offset;
par.VoxelSize   = vx*[1 1 1];
par.Supersample = opt.supersample;
par.PVESteps    = opt.pve_steps;
par.PVERange    = opt.pve_range;

try
  lmeta = struct('GeneratedBy',{{gen}});
  lmeta.Manual = false;
  lmeta.Description = ['Ground truth partial volume label image of a folded sphere ' ...
                       'with a constant GM thickness. Values are continuous and ' ...
                       'interpolate between the tissue labels, e.g. 2.5 is an equal ' ...
                       'mixture of GM and WM.'];
  lmeta.Labels = struct('CSF',1,'GM',2,'WM',3);
  lmeta.PhantomParameters = par;
  spm_jsonwrite(regexprep(files.label,'\.nii$','.json'), lmeta);

  imeta = struct('GeneratedBy',{{gen}});
  imeta.Sources = {['bids::' spm_file(files.label,'filename')]};
  ipar = par;
  ipar.Intensities = opt.intensities(:)';
  imeta.Description = ['Ideal T1w image of the phantom, i.e. the tissue intensities ' ...
                       'mixed by the partial volume fractions. It carries no noise, ' ...
                       'no bias field and no other artefact, use mri_simulate for those.'];
  imeta.PhantomParameters = ipar;
  spm_jsonwrite(regexprep(files.t1w,'\.nii$','.json'), imeta);
catch ME
  fprintf('Warning: Failed to write JSON sidecar(s): %s\n', ME.message);
end


%==========================================================================
% function label = bids_label(str)
%
% Purpose
%   Turn a tag into a valid BIDS entity label, which only allows letters and
%   digits. Same convention as in mri_simulate, so that the names of the two
%   tools match.
%
% Inputs
%   str - char: arbitrary tag.
%
% Output
%   label - char: str with decimal points replaced by 'p' (0.5 -> 0p5) and
%           all remaining non-alphanumeric characters removed.
%==========================================================================
function label = bids_label(str)

label = strrep(str, '.', 'p');
label = regexprep(label, '[^a-zA-Z0-9]', '');
