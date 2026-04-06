%RUN_VALIDATION Quick validation script for NAUTILUS fixes
% Run this in MATLAB to test the fixes

clc; clear; close all;

fprintf('===============================================\n');
fprintf('  NAUTILUS Validation\n');
fprintf('===============================================\n\n');

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'Classes'));
addpath(fullfile(repoRoot, 'Hulls'));
cd(repoRoot);

%% Generate fresh meshes
fprintf('1. Generating Box Barge mesh...\n');
cd(fullfile(repoRoot, 'Hulls'));
generate_box_barge();
cd(repoRoot);

%% Test Box Barge hydrostatics
fprintf('\n2. Testing Box Barge hydrostatics...\n');
loader = GeometryLoader(struct());
[mesh, ~, ~] = loader.loadHullMesh('BoxBarge_mesh');

L = 10; B = 4;
fprintf('   Draft    V_calc    V_theo    Error%%    Af_calc   KB_calc   KB_theo\n');
fprintf('   ---------------------------------------------------------------\n');

for T = [0.5, 1.0, 1.5]
    COND = struct('T', T, 'tetaT', 0, 'tetaL', 0, 'Lwave', 0, 'Hwave', 0, 'Xcrest', 0, 'Wavetype', 2);
    [COND1, ~] = Algorithms.geometriaMesh(COND, mesh, 0);
    
    V_theo = L * B * T;
    V_calc = COND1.V;
    Af_calc = COND1.Af;
    KB_calc = COND1.B(3);
    KB_theo = T/2;
    
    fprintf('   %.1f      %.2f      %.2f      %.1f%%      %.2f      %.3f     %.3f\n', ...
        T, V_calc, V_theo, abs(V_calc-V_theo)/V_theo*100, Af_calc, KB_calc, KB_theo);
end

%% Test full pipeline on Box Barge
fprintf('\n3. Running full pipeline on Box Barge...\n');
cfg = struct('computeGZ', true, 'tetaT_vec', 0:15:45, 'refDraftIdx', 2, ...
    'tolVolume', 1e-5, 'tolTrim', 1e-4, 'rho', 1025, 'savePlots', false, ...
    'Lwave', 0, 'Hwave', 0, 'Xcrest', 0, 'Wavetype', 2, ...
    'autoSaveMesh', false, 'meshOptions', struct(), 'sanity', struct());
cfg.sanity.areaTolerance = 1e-12;
cfg.sanity.removeDegenerateFaces = true;
cfg.sanity.requireWatertight = false;

try
    pipeline = Pipeline(cfg);
    pipeline.runForHull('BoxBarge_mesh', 0.5:0.5:1.5);
    fprintf('   [OK] Pipeline completed\n');
catch ME
    fprintf('   [FAIL] %s\n', ME.message);
end

fprintf('\n===============================================\n');
fprintf('  Validation Complete!\n');
fprintf('===============================================\n');
