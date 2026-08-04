function res = thickness_phantom_eval(P, opt)
% THICKNESS_PHANTOM_EVAL - Accuracy of a thickness phantom and of PBT on it
%
% Overview:
%   `thickness_phantom_eval` checks a phantom of `thickness_phantom` in two
%   independent ways.
%
%   1. Fidelity of the phantom. The label image is the only thing a thickness
%      pipeline ever sees, so the question is not whether the geometry that
%      was requested is correct, but whether the label image at its voxel
%      size still carries it. Two things are measured:
%      - The distance between the two cortical boundaries, recovered from the
%        label image alone. The label is supersampled, thresholded at 2.5 and
%        at 1.5, and cat_bwdist gives the euclidean distance from the WM
%        surface, which is read out on the GM/CSF surface (and the other way
%        round). This is a purely geometric measure and does not model
%        anything, thus a deviation is an error of the phantom and not of a
%        thickness method.
%      - The tissue volumes, compared against the volumes of the supersampled
%        grid the phantom was built on. They test the partial volume label as
%        a whole, including the voxels that no boundary passes through.
%
%   2. Accuracy of a thickness measure on that phantom, optionally with CAT's
%      projection based thickness cat_vol_pbtsimple, which takes exactly this
%      kind of partial volume label as input. Its thickness map is compared
%      against the known thickness, both over all GM voxels and on the
%      central surface (percentage position 0.5), where CAT samples it.
%
%   The two parts answer different questions and both are needed: without the
%   first one an error of the thickness measure cannot be separated from an
%   error of the phantom.
%
% Syntax:
%   res = thickness_phantom_eval(P, opt)
%
% Parameters:
%   P: The phantom, given as one of
%       - char/cellstr: filename(s) of the label image (..._dseg.nii). The
%         thickness and the voxel size are read from the JSON sidecar that
%         thickness_phantom wrote, so nothing else has to be passed.
%       - single/double array: the label image itself, then opt.thickness and
%         opt.vx are required.
%       - struct: the info struct of thickness_phantom, which names the files.
%
%   opt (struct): Options. Defaults are applied for missing fields.
%       - 'thickness' (double): Known thickness in mm. Only needed if it
%         cannot be read from a sidecar. Default: NaN.
%       - 'vx' (double): Isotropic voxel size in mm, same. Default: NaN.
%       - 'supersample' (integer): Factor for the geometric measurement. The
%         boundaries are recovered from a thresholded label image, so the
%         accuracy of that measurement is about half a sample of this grid.
%         Even values are increased by one. Default: 3.
%       - 'pbt' (logical): Also run cat_vol_pbtsimple. Default: 1.
%       - 'pbtopt' (struct): Options passed on to cat_vol_pbtsimple.
%         Default: struct() (its own defaults).
%       - 'fig' (char): Filename of a QC figure (PNG) with a central slice
%         and the thickness histograms. Default: '' (no figure).
%       - 'verbose' (logical): Print the report. Default: 1.
%
% Outputs:
%   res - struct array with one entry per phantom:
%       .thickness              known thickness in mm
%       .vx                     voxel size in mm
%       .geom.outer             stats of the distance from the WM surface,
%                               read on the GM/CSF surface
%       .geom.inner             stats of the distance from the GM/CSF
%                               surface, read on the WM surface
%       .volume                 tissue volumes of the label and, if the
%                               sidecar is available, of the reference grid
%       .pbt.gm / .pbt.central  stats of cat_vol_pbtsimple over all GM voxels
%                               and on the central surface
%     Every stats entry holds n, mean, median, sd, bias (mean - thickness),
%     rmse, and the 5th/95th percentile.
%
% Examples:
%   Example 1 - Create a phantom and evaluate it. Note that info is the third
%               output of thickness_phantom, the first two are the images:
%       [~, ~, info] = thickness_phantom;
%       res = thickness_phantom_eval(info);
%
%   Example 2 - Evaluate a series of thickness values and write a QC figure:
%       [~, ~, info] = thickness_phantom(struct('thickness', 1.5:0.5:3.5));
%       res = thickness_phantom_eval(info, struct('fig','phantom_qc.png'));
%
%   Example 3 - Only the geometry of an existing label image:
%       res = thickness_phantom_eval('thicknessPhantom_res-0p5mm_desc-Thickness25mm_dseg.nii', ...
%                                    struct('pbt',0));
%
% See also: thickness_phantom, cat_vol_pbtsimple, mri_simulate

if ~exist('cat_bwdist','file')
  error('CAT >= 26 has to be in the MATLAB path to use thickness_phantom_eval.')
end

def.thickness   = NaN;
def.vx          = NaN;
def.supersample = 3;
def.pbt         = 1;
def.pbtopt      = struct();
def.fig         = '';
def.verbose     = 1;

if nargin < 2, opt = def;
else, opt = cat_io_checkinopt(opt, def); end

% Collect the phantoms to evaluate. Everything is turned into a list of
% entries that hold the label image, the voxel size and the known thickness.
job = collect_input(P, opt);

ss = max(1, round(opt.supersample));
if mod(ss,2) == 0, ss = ss + 1; end

res = struct('name',{},'thickness',{},'vx',{},'geom',{},'volume',{},'pbt',{});

for i = 1:numel(job)
  Yp0 = job(i).Yp0;
  vx  = job(i).vx;
  t   = job(i).thickness;

  if ~isfinite(t)
    error(['The thickness of %s is unknown. Pass it as opt.thickness or keep the ' ...
           'JSON sidecar next to the label image.'], job(i).name);
  end

  res(i).name      = job(i).name;
  res(i).thickness = t;
  res(i).vx        = vx;

  %------------------------------------------------------------------------
  % 1a. Geometry: distance between the two boundaries of the label image
  %------------------------------------------------------------------------
  if opt.verbose
    fprintf('\n%s\n', repmat('=',1,80));
    fprintf('%s (thickness %.2f mm, %.2f mm voxels)\n', job(i).name, t, vx);
    fprintf('%s\n', repmat('=',1,80));
    fprintf('Recover the boundaries at %.3f mm and measure their distance.\n', vx/ss);
  end

  % Both boundaries are recovered from the label image alone: 2.5 is the
  % WM/GM boundary of a partial volume label and 1.5 the GM/CSF boundary.
  % The supersampling is a plain trilinear interpolation, so it does not add
  % any geometry, it only lets the distance map below resolve the boundary
  % better than the voxel size of the phantom would allow.
  [Yf, vxf] = upsample_label(Yp0, vx, ss);

  % Distance from the WM surface. The 0.5 sample correction puts the surface
  % on the face of the last object voxel instead of on its centre, the same
  % convention that mri_simulate uses for its cortical band.
  Dwm = cat_bwdist(single(Yf >= 2.5), vxf*[1 1 1]) - 0.5*vxf;

  % Read it on the GM/CSF surface. That surface is not taken as a layer of
  % voxels: a layer has a thickness of its own and its centres sit somewhere
  % inside it, which biases the result by a good part of a sample. The
  % surface is instead located between the samples, at the crossings of the
  % isolevel 1.5 along the three axes, and the distance is interpolated
  % there, so the readout has no bias of its own.
  d_outer = sample_at_isolevel(Yf, 1.5, Dwm);
  res(i).geom.outer = stats(d_outer, t);
  % the samples themselves are only kept for the histogram of the QC figure,
  % they are far too many to return them with the result
  if ~isempty(opt.fig), job(i).d_outer = d_outer; end
  clear Dwm d_outer

  % The same in the other direction, which is not redundant. The distance
  % from the WM outwards is exactly the thickness by construction, but the
  % distance from the GM/CSF surface back to the WM is only the same as long
  % as there is CSF on the other side of the sulcus. Where the folds are
  % narrower than twice the thickness the two banks of GM touch, the CSF is
  % squeezed out, and this direction measures how far the nearest remaining
  % CSF is. It is the situation that makes a thickness measure overestimate,
  % thus the phantom should report how much of it it contains.
  Dcsf = cat_bwdist(single(Yf < 1.5), vxf*[1 1 1]) - 0.5*vxf;
  d_inner = sample_at_isolevel(Yf, 2.5, Dcsf);
  res(i).geom.inner  = stats(d_inner, t);
  res(i).geom.buried = mean(d_inner > t + 0.5*vx);
  clear Dcsf Yf d_inner

  %------------------------------------------------------------------------
  % 1b. Tissue volumes of the partial volume label
  %------------------------------------------------------------------------
  % Yp0 mixes the tissues linearly, thus the fractions can be read back from
  % it: a value of 2.3 is 0.7 GM and 0.3 WM. Only three classes and the
  % background are involved and they are ordered, so the fraction of class k
  % is the usual triangular weight.
  Yp0toC = @(Y,c) max(0, 1-abs(Y-c));
  res(i).volume.label = struct('csf', sum(sum(sum(Yp0toC(Yp0,1))))*vx^3/1000, ...
                               'gm',  sum(sum(sum(Yp0toC(Yp0,2))))*vx^3/1000, ...
                               'wm',  sum(sum(sum(Yp0toC(Yp0,3))))*vx^3/1000);
  res(i).volume.ref = job(i).volume_ref;

  %------------------------------------------------------------------------
  % 2. Thickness measure
  %------------------------------------------------------------------------
  res(i).pbt = [];
  if opt.pbt
    if ~exist('cat_vol_pbtsimple','file')
      fprintf('Warning: cat_vol_pbtsimple was not found, the thickness measure is skipped.\n');
    else
      if opt.verbose, fprintf('Run cat_vol_pbtsimple.\n'); end
      [Ygmt, Ypp] = cat_vol_pbtsimple(single(Yp0), vx*[1 1 1], opt.pbtopt);

      % Over the whole GM band and, separately, on the central surface, which
      % is where CAT samples the thickness for the surface based analysis.
      gm      = Ygmt > 0 & Yp0 > 1.5 & Yp0 < 2.5;
      central = Ygmt > 0 & Ypp > 0.45 & Ypp < 0.55;
      res(i).pbt.gm      = stats(double(Ygmt(gm)), t);
      res(i).pbt.central = stats(double(Ygmt(central)), t);
      job(i).Ygmt = Ygmt;
      clear Ygmt Ypp gm central
    end
  end

  if opt.verbose, print_report(res(i)); end
end

if ~isempty(opt.fig)
  make_figure(job, res, opt.fig);
  fprintf('\nSave %s\n', opt.fig);
end


%==========================================================================
% function job = collect_input(P, opt)
%
% Purpose
%   Turn the different ways of naming a phantom into one list of entries that
%   all hold the label image, the voxel size, the known thickness and, where
%   it is available, the reference volumes of the supersampled grid.
%
% Inputs
%   P   - char/cellstr filename(s), numeric label image, or the info struct
%         of thickness_phantom.
%   opt - struct with the fallback fields thickness and vx.
%
% Output
%   job - struct array with the fields name, Yp0, vx, thickness, volume_ref.
%==========================================================================
function job = collect_input(P, opt)

job = struct('name',{},'Yp0',{},'vx',{},'thickness',{},'volume_ref',{},'Ygmt',{}, ...
             'd_outer',{});

if isstruct(P) && isfield(P,'files')
  % info struct of thickness_phantom
  names = arrayfun(@(x) x.files.label, P, 'UniformOutput', false);
  if any(cellfun(@isempty, names))
    error(['The phantom was created with write=0, thus there is no label image to ' ...
           'read.\n  Pass the label image itself instead, i.e. the first output of ' ...
           'thickness_phantom,\n  together with opt.vx and opt.thickness.'], '');
  end
  for i = 1:numel(names)
    job(i) = read_label(names{i}, opt);
    % the info struct knows the reference volumes even without a sidecar
    if isfield(P(i),'volume_ref'), job(i).volume_ref = P(i).volume_ref; end
    if isfield(P(i),'thickness') && isfinite(P(i).thickness)
      job(i).thickness = P(i).thickness;
    end
  end
  return
end

if isnumeric(P) || islogical(P)
  job(1).name       = 'label image';
  job(1).Yp0        = single(P);
  job(1).vx         = opt.vx;
  job(1).thickness  = opt.thickness;
  job(1).volume_ref = [];
  job(1).Ygmt       = [];
  job(1).d_outer    = [];
  if ~isfinite(job(1).vx)
    error('opt.vx is required when the label image is passed directly.')
  end
  return
end

if ischar(P), P = cellstr(P); end

% A series may be given as a vector of thickness values, one per entry, which
% is what thickness_phantom is called with for a series of phantoms.
tt = opt.thickness;
if ~isscalar(tt) && numel(tt) ~= numel(P)
  error('opt.thickness has %d values but %d phantoms were given.', numel(tt), numel(P));
end

for i = 1:numel(P)
  if isscalar(tt), opt_i = opt; else, opt_i = opt; opt_i.thickness = tt(i); end
  if isnumeric(P{i}) || islogical(P{i})
    % label images that are passed directly, e.g. the first output of
    % thickness_phantom for several thickness values
    job(i)      = collect_input(P{i}, opt_i);
    job(i).name = sprintf('label image %d', i);
  else
    job(i) = read_label(deblank(P{i}), opt_i);
  end
end


%==========================================================================
% function e = read_label(fname, opt)
%
% Purpose
%   Read a label image and everything the JSON sidecar of thickness_phantom
%   knows about it.
%
% Inputs
%   fname - char: filename of the label image.
%   opt   - struct with the fallback fields thickness and vx.
%
% Output
%   e - struct entry for the job list of collect_input.
%==========================================================================
function e = read_label(fname, opt)

V = spm_vol(spm_file(fname,'number',''));
e.name       = spm_file(fname,'filename');
e.Yp0        = single(spm_read_vols(V));
vx           = sqrt(sum(V.mat(1:3,1:3).^2));
e.vx         = mean(vx);
e.thickness  = opt.thickness;
e.volume_ref = [];
e.Ygmt       = [];
e.d_outer    = [];

if any(abs(vx - e.vx) > 1e-6)
  fprintf('Warning: %s is not isotropic (%g %g %g mm), the mean is used.\n', e.name, vx);
end

% The sidecar carries the geometry of the phantom, i.e. the thickness that
% the evaluation compares against. An explicit opt.thickness wins over it.
json = regexprep(fname,'\.nii(\.gz)?$','.json');
if exist(json,'file')
  try
    meta = spm_jsonread(json);
    if ~isfinite(e.thickness) && isfield(meta,'PhantomParameters') && ...
       isfield(meta.PhantomParameters,'Thickness')
      e.thickness = meta.PhantomParameters.Thickness;
    end
  catch ME
    fprintf('Warning: could not read %s: %s\n', json, ME.message);
  end
end

if isfinite(opt.thickness), e.thickness = opt.thickness; end


%==========================================================================
% function [Yf, vxf] = upsample_label(Yp0, vx, ss)
%
% Purpose
%   Resample the partial volume label on a grid that is supersampled by ss,
%   so that the distance map built from it resolves the tissue boundaries
%   better than the voxel size of the phantom does.
%
% Inputs
%   Yp0 - single(dim): partial volume label, CSF=1, GM=2, WM=3.
%   vx  - double: isotropic voxel size in mm.
%   ss  - integer: odd supersampling factor.
%
% Output
%   Yf  - single: the label on the fine grid.
%   vxf - double: voxel size of the fine grid in mm.
%==========================================================================
function [Yf, vxf] = upsample_label(Yp0, vx, ss)

vxf = vx/ss;
if ss == 1, Yf = Yp0; return; end

d = size(Yp0);
% Coordinates of the fine grid in voxels of the coarse grid, using the same
% centring as thickness_phantom, i.e. the fine sample (i-1)*ss+(ss+1)/2 sits
% exactly on the coarse voxel i.
f = @(n) ( (0:n*ss-1) - (ss-1)/2 )/ss + 1;
[Xi, Yi, Zi] = ndgrid(single(f(d(1))), single(f(d(2))), single(f(d(3))));

% Outside of the volume interp3 returns the fill value, which is background
% here. The phantom keeps a margin to the border, so nothing is lost.
Yf = interp3(Yp0, Yi, Xi, Zi, 'linear', 0);   % note the ndgrid/meshgrid swap


%==========================================================================
% function v = sample_at_isolevel(Yf, level, Dmap)
%
% Purpose
%   Read a map on an isosurface of a second map with sub-voxel accuracy.
%   Every pair of neighbouring samples along the three axes whose values
%   straddle the isolevel holds one crossing. Its position between the two
%   samples follows from linear interpolation of Yf, and Dmap is evaluated
%   there by linear interpolation as well.
%
%   This is used instead of a layer of voxels around the isosurface, because
%   such a layer has a thickness of its own: its centres sit a fraction of a
%   sample inside the surface, and that fraction depends on the orientation
%   of the surface, which biases the distance that is read on it.
%
% Inputs
%   Yf    - single: the map whose isosurface defines where to sample.
%   level - double: the isolevel.
%   Dmap  - single: the map to read, same size as Yf.
%
% Output
%   v - double column vector, one entry per crossing.
%==========================================================================
function v = sample_at_isolevel(Yf, level, Dmap)

d = size(Yf);
B = Yf >= level;
stride = [1 d(1) d(1)*d(2)];
v = cell(1,3);

for k = 1:3
  % Mark the sample of every neighbouring pair along axis k that straddles
  % the isolevel. The mask is built on the full grid so that its linear
  % indices point into Yf directly.
  C = false(d);
  switch k
    case 1, C(1:end-1,:,:) = xor(B(1:end-1,:,:), B(2:end,:,:));
    case 2, C(:,1:end-1,:) = xor(B(:,1:end-1,:), B(:,2:end,:));
    case 3, C(:,:,1:end-1) = xor(B(:,:,1:end-1), B(:,:,2:end));
  end
  i1 = find(C);
  i2 = i1 + stride(k);

  % The two values straddle the level, thus they differ and the weight is
  % well defined and lies in [0,1).
  a = double(Yf(i1)); b = double(Yf(i2));
  w = (level - a)./(b - a);
  v{k} = double(Dmap(i1)) + w.*(double(Dmap(i2)) - double(Dmap(i1)));
end

v = vertcat(v{:});


%==========================================================================
% function s = stats(v, t)
%
% Purpose
%   Summarize a set of thickness values against the known thickness.
%
% Inputs
%   v - double vector: measured values in mm.
%   t - double: known thickness in mm.
%
% Output
%   s - struct with n, mean, median, sd, bias, rmse, p5 and p95.
%==========================================================================
function s = stats(v, t)

v = v(isfinite(v));
s.n = numel(v);
if isempty(v)
  [s.mean, s.median, s.sd, s.bias, s.rmse, s.p5, s.p95] = deal(NaN);
  return
end
sv       = sort(v);
s.mean   = mean(v);
s.median = sv(max(1,round(0.500*numel(sv))));
s.sd     = std(v);
s.bias   = s.mean - t;
s.rmse   = sqrt(mean((v - t).^2));
s.p5     = sv(max(1,round(0.05*numel(sv))));
s.p95    = sv(max(1,round(0.95*numel(sv))));


%==========================================================================
% function print_report(r)
%
% Purpose
%   Print the result of one phantom as a table.
%
% Inputs
%   r - struct: one entry of the result of thickness_phantom_eval.
%==========================================================================
function print_report(r)

fprintf('\n  Known thickness %.3f mm\n\n', r.thickness);
fprintf('  %-34s %6s %8s %8s %8s %8s %8s\n', ...
        'measure','n','mean','sd','bias','rmse','p5-p95');
fprintf('  %s\n', repmat('-',1,78));

row = @(nam,s) fprintf('  %-34s %6d %8.3f %8.3f %+8.3f %8.3f %6.2f-%.2f\n', ...
                       nam, s.n, s.mean, s.sd, s.bias, s.rmse, s.p5, s.p95);

row('geometry: WM surf -> GM/CSF surf', r.geom.outer);
row('geometry: GM/CSF surf -> WM surf', r.geom.inner);
if ~isempty(r.pbt)
  row('cat_vol_pbtsimple, all GM',      r.pbt.gm);
  row('cat_vol_pbtsimple, central surf',r.pbt.central);
end
fprintf('\n  %.1f%% of the WM surface faces a sulcus without CSF (buried sulcus).\n', ...
        100*r.geom.buried);

fprintf('\n  %-10s %10s %10s %10s\n', 'volume/ml', 'CSF', 'GM', 'WM');
fprintf('  %-10s %10.2f %10.2f %10.2f\n', 'label', ...
        r.volume.label.csf, r.volume.label.gm, r.volume.label.wm);
if ~isempty(r.volume.ref)
  fprintf('  %-10s %10.2f %10.2f %10.2f\n', 'reference', ...
          r.volume.ref.csf, r.volume.ref.gm, r.volume.ref.wm);
  fprintf('  %-10s %9.2f%% %9.2f%% %9.2f%%\n', 'error', ...
          100*(r.volume.label.csf/r.volume.ref.csf - 1), ...
          100*(r.volume.label.gm /r.volume.ref.gm  - 1), ...
          100*(r.volume.label.wm /r.volume.ref.wm  - 1));
end


%==========================================================================
% function make_figure(job, res, fname)
%
% Purpose
%   Write a QC figure with a central slice of the label image, of the
%   thickness map where one exists, and the histogram of the measured
%   thickness values against the known thickness.
%
% Inputs
%   job   - struct array of collect_input, with the label and thickness maps.
%   res   - struct array of results.
%   fname - char: filename of the PNG to write.
%==========================================================================
function make_figure(job, res, fname)

n = numel(job);
h = figure('Visible','off','Color','w','Position',[0 0 1000 260*n]);

% Crop all panels to the bounding box of the largest phantom, so that the
% slices of a series stay comparable and no space is spent on background.
d  = size(job(1).Yp0);
sl = round(d(3)/2);
bb = [d(1) 1 d(2) 1];
for i = 1:n
  S = squeeze(job(i).Yp0(:,:,sl)) > 0;
  [r1, c1] = find(S);
  bb = [min(bb(1),min(r1)) max(bb(2),max(r1)) min(bb(3),min(c1)) max(bb(4),max(c1))];
end
bb = [max(1,bb(1)-2) min(d(1),bb(2)+2) max(1,bb(3)-2) min(d(2),bb(4)+2)];
crop = @(Y) rot90(Y(bb(1):bb(2), bb(3):bb(4)));

for i = 1:n
  t = res(i).thickness;

  subplot(n,3,(i-1)*3+1);
  imagesc(crop(squeeze(job(i).Yp0(:,:,sl))), [0 3]); axis image off; colormap(gca,gray);
  title(sprintf('label, known %.2f mm', t), 'FontSize', 9);

  subplot(n,3,(i-1)*3+2);
  if isempty(job(i).Ygmt)
    axis off; title('no thickness map', 'FontSize', 9);
  else
    % Show the thickness only where it was estimated, i.e. inside the GM
    % band. Everything else stays white, otherwise the background and the WM
    % interior dominate the colour scale and hide the variation that matters.
    G = crop(squeeze(job(i).Ygmt(:,:,sl)));
    M = crop(squeeze(job(i).Yp0(:,:,sl))) > 1.25 & ...
        crop(squeeze(job(i).Yp0(:,:,sl))) < 2.75 & G > 0;
    im = imagesc(G, t + 0.4*[-1 1]);
    set(im,'AlphaData',double(M)); set(gca,'Color','w');
    axis image off; colormap(gca,jet); colorbar;
    title('cat\_vol\_pbtsimple', 'FontSize', 9);
  end

  subplot(n,3,(i-1)*3+3); hold on; box on;
  edges = t + linspace(-0.8,0.8,80);
  ctr   = (edges(1:end-1)+edges(2:end))/2;
  leg   = {};
  plot(ctr, histcounts(job(i).d_outer, edges, 'Normalization','probability'), ...
       'k', 'LineWidth', 1);
  leg{end+1} = 'geometry'; %#ok<AGROW>
  if ~isempty(job(i).Ygmt)
    v = job(i).Ygmt(job(i).Ygmt>0 & job(i).Yp0>1.5 & job(i).Yp0<2.5);
    plot(ctr, histcounts(double(v(:)), edges, 'Normalization','probability'), 'b', 'LineWidth', 1);
    leg{end+1} = 'pbtsimple'; %#ok<AGROW>
  end
  yl = ylim; plot(t*[1 1], yl, 'r--'); ylim(yl); leg{end+1} = 'known'; %#ok<AGROW>
  xlim([edges(1) edges(end)]);
  xlabel('thickness (mm)'); ylabel('fraction');
  legend(leg,'Location','northwest','FontSize',7); legend boxoff
  if isempty(res(i).pbt)
    title(sprintf('geometry bias %+.3f mm', res(i).geom.outer.bias), 'FontSize', 9);
  else
    title(sprintf('geometry %+.3f, pbt %+.3f mm', ...
                  res(i).geom.outer.bias, res(i).pbt.gm.bias), 'FontSize', 9);
  end
end

print(h, fname, '-dpng', '-r110');
close(h);
