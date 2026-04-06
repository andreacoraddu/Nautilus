%==========================================================================
% NAUTILUS PROJECT - Main entry point for the single-class OOP workflow
%==========================================================================
% This script initializes the execution environment, defines the user-level
% configuration, and launches the hydrostatics/stability pipeline for one
% or more hull geometries.
%==========================================================================

tic
clc; clear; close all;

% Resolve the repository root from the current script location so that the
% project remains portable across machines and folder structures.
repoRoot = fileparts(mfilename('fullpath'));

% Add the main source folders to the MATLAB path:
% - Classes: object-oriented implementation of the processing pipeline
% - Hulls:   geometry files / mesh definitions used as input cases
addpath(fullfile(repoRoot, 'Classes'));
addpath(fullfile(repoRoot, 'Hulls'));

% Force MATLAB to operate from the repository root to ensure that relative
% paths used later in the workflow remain consistent.
cd(repoRoot);

%% ========================================================================
% USER CONFIGURATION
% ========================================================================
% List of hull geometries to be processed.
% Each entry should correspond to a hull mesh or geometry definition
% available in the Hulls folder.
hulls = {'BoxBarge_mesh'};

% Draft sweep definition.
% The bundled Wigley benchmark has:
%   - characteristic depth D = 0.35 m
%   - design draft T = 0.25 m
% The selected range below is therefore kept within a physically meaningful
% interval for this reference hull.
T_start = 0.10;
T_end   = 0.25;
T_step  = 0.01;

% Global configuration structure passed to the Pipeline object.
cfg = struct();

% Stability / righting-arm options
% If computeGZ is enabled, GZ curves are evaluated over the heel-angle
% vector tetaT_vec. Angles are expressed in degrees.
cfg.computeGZ = false;
cfg.tetaT_vec = 0:2.5:60;

% Reference draft index used internally for selecting the baseline draft
% condition in post-processing or comparative analyses.
cfg.refDraftIdx = 4;

% Numerical tolerances
% tolVolume: acceptable mismatch in displaced volume convergence
% tolTrim:   acceptable tolerance in trim equilibrium convergence
cfg.tolVolume = 1e-5;
cfg.tolTrim   = 1e-4;

% Water density [kg/m^3], typically seawater.
cfg.rho = 1025;

% Output control
% If true, figures generated during the analysis are automatically saved.
cfg.savePlots = true;

% Wave-condition parameters
% These are currently set to calm-water conditions:
%   Lwave   = wavelength
%   Hwave   = wave height
%   Xcrest  = crest longitudinal position
%   Wavetype = wave model selector / identifier used internally
cfg.Lwave   = 0;
cfg.Hwave   = 0;
cfg.Xcrest  = 0;
cfg.Wavetype = 2;

% Mesh output options
% autoSaveMesh: if true, derived or processed meshes are exported
% meshOptions.plotSubmerged: controls visualization of the submerged hull
cfg.autoSaveMesh = false;
cfg.meshOptions = struct();
cfg.meshOptions.plotSubmerged = false;

% Geometry sanity-check settings
% areaTolerance: threshold below which face areas are treated as degenerate
% removeDegenerateFaces: remove invalid or near-zero-area elements
% requireWatertight: enforce watertightness if strict mesh integrity is needed
cfg.sanity = struct();
cfg.sanity.areaTolerance = 1e-12;
cfg.sanity.removeDegenerateFaces = true;
cfg.sanity.requireWatertight = false;

%% ========================================================================
% MAIN EXECUTION PIPELINE
% ========================================================================
% Build the draft vector used for the hydrostatic sweep.
T_vec = T_start:T_step:T_end;

% Instantiate the main pipeline object with the user-defined configuration.
pipeline = Pipeline(cfg);

% Run the analysis for each hull listed in the input set.
for k = 1:numel(hulls)
    pipeline.runForHull(hulls{k}, T_vec);
end

% Final execution summary
fprintf('\n+----------------------------------------------+\n');
fprintf('  NAUTILUS pipeline completed in %7.2f s\n', toc);
fprintf('+----------------------------------------------+\n\n');