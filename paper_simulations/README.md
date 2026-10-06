# Paper simulations

These scripts are grouped by the output they produce. Run them from the package
root unless a script says otherwise. The existing results below were moved with
their scripts; reorganizing this folder did not rerun any simulation.

| Folder | Contents |
|:--|:--|
| `10sLimit/` | The 12 Gaussian and Poisson scripts used for the 10-second comparison table. Existing logs, CSV summaries, and the table source/PDF are in `10sLimit/results/`. |
| `Figure1_2/` | MeanVar figure script using the `dust` package. |
| `PruningIllustration/` | Four-panel pruning figure, script, PDF, and PNG. |
| `DUSTvsDUSTib/` | Mathematical comparison report and PDF. |
| `MethodComparison10k/` | DUST/DUSTib backend comparison script and CSV results. |
| `IntroductionTimings/` | Earlier introduction timing scripts, separate from the 12 table scripts. |

The `10sLimit` scripts write CSV summaries to `10sLimit/results/` when run
from the package root. The table in that results folder records the existing
measurements; it is not regenerated automatically from the CSV files.

Run `Rscript paper_simulations/Figure1_2/Figure_1_2_meanVar.R` from the package
root to reproduce Figure 1.2. The PNG, PDF, and saved summary are written beside
the script. Use `Rscript paper_simulations/Figure1_2/Figure_1_2_meanVar.R --plot-only`
to redraw the image from that summary without rerunning the simulation.
