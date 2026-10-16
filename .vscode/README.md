# Using NAUTILUS in VS Code — free vs. licensed

Short answer: **the VS Code side is free, the MATLAB engine is not.**
The `MathWorks.language-matlab` extension is open source (MIT) and free to
install, but every useful feature it offers (code analysis, go-to-definition,
formatting, "Run", debugging) is served by a MATLAB installation on your own
machine. That engine needs a licence (campus-wide licence, employer licence,
student licence, or a paid personal licence). With no licence the extension
degrades to plain text editing, so installing it alone buys you nothing.

Free routes that actually run this repository:

| Option | Cost | Runs NAUTILUS? | Notes |
|---|---|---|---|
| **GNU Octave** + free Octave extension in VS Code | Free (GPL) | After patching a few MATLAB-only calls (see below) | Best free option; Octave 9+ recommended, Octave 11 has the best `classdef` support |
| **MATLAB Online** (MathWorks account) | Free tier exists (limited hours/month) | Yes — it is real MATLAB | Runs in the browser, not inside VS Code; good for occasional runs and for checking results against Octave |
| **MATLAB in VS Code** (`MathWorks.language-matlab`) | Extension free, **engine paid** | Yes | The only way to get the real MATLAB experience inside VS Code; needs a local licensed MATLAB R2021b+ |
| Campus / employer MATLAB licence | Usually free for students/staff | Yes | Ask your IT/licensing office first — most universities already have a site licence |

## Free setup with GNU Octave

1. Install GNU Octave (Linux: your distro package or the flatpak; Windows/macOS:
   <https://octave.org/download>). Verify with `octave --version`.
2. In VS Code install the **Octave Language** extension (`leafvmaple.octave`),
   which is already listed in `.vscode/extensions.json`.
3. Uncomment the `files.associations` line in `.vscode/settings.json` so `.m`
   files are highlighted as Octave.
4. Run the pipeline with the task **"Octave: run main.m (free)"**
   (`Terminal ▸ Run Task…`), which is equivalent to:

   ```bash
   octave --no-gui --quiet main.m
   ```

## What needs patching before NAUTILUS runs under Octave

Checked against the current Octave documentation (`wiki.octave.org/Classdef`,
`wiki.octave.org/Differences_between_Octave_and_Matlab`, `octave.org/NEWS-11`)
and a static scan of this repository:

| MATLAB construct | Octave status | Where it is used here |
|---|---|---|
| `tiledlayout` / `nexttile` | **Not implemented** | `Classes/GeometryChecker.m` (3 tiles), `Classes/PostProcessor.m` (6 tiles) → replace with `subplot` |
| `'Interpreter','latex'` in labels/titles/legends | **Not supported** (Octave has `tex`/`none` only) | ~40 calls in `GeometryChecker.m` and `PostProcessor.m` → labels would print raw LaTeX |
| `arguments` blocks with validators (`hullName (1,:) char`) | Argument **validation is not implemented** | `Classes/Pipeline.m:30`, `Classes/GeometryLoader.m:76` → convert to `validateattributes` calls |
| Property validation (`properties; cfg struct; end`) | **Not implemented** | `Classes/PostProcessor.m` property block |
| `classdef` in general | *Limited* support; `Access` attributes partial, `Constant`/`Static`/`Access` methods OK, `enumeration`/`events` absent, `save`/`load` of classdef objects unsupported | all six files in `Classes/` |
| `save(...)` of results (`Output/Hydro_*.mat`) | Octave writes MAT v4/v6/v7 only (**no v7.3**); classdef objects are saved as plain structs | `Classes/PostProcessor.m:34` (structs only, so OK) |
| `onCleanup` | Supported since Octave 3.4 | `Classes/*` (11 uses) — fine |
| `validateattributes`, `strjoin`, `contains`, `startsWith`, `vecnorm`, `interp1`, `hypot` | Available in modern Octave | fine |
| `print(fig, ..., '-dpdf', '-painters')` | Available | fine |

Two practical notes for the free path:

* **Headless machines**: `main.m` saves figures (`cfg.savePlots = true`), so on a
  machine without a display use `graphics_toolkit('gnuplot')` (or run with
  `--no-window-system` only if you set `cfg.savePlots = false`).
* **Numerics**: hydrostatic results from Octave can differ in the last digits
  from MATLAB. Always re-run the benchmarks in the README's validation table
  before trusting Octave output.

## Files in this folder

* `extensions.json` — extension recommendations (MathWorks + free Octave).
* `settings.json` — editor defaults and the optional `.m → octave` association.
* `tasks.json` — ready-made run tasks for MATLAB (`-batch`) and Octave.

Nothing here changes the MATLAB code itself, so committing these files is safe
for MATLAB users; the Octave path is opt-in via the commented association.
