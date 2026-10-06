# dust 0.2.0

- Accelerated the Gaussian multivariate exact solver with verified active-face
  searches and exhaustive fallback.
- Renamed the multivariate search budget from nbLoops to nbIterations.
  Added optional epsilon stopping by decision-function gain for coordinate
  descent, projected gradient, and quasi-Newton searches.

# dust 0.1.0

- Added eight one-parameter models and Gaussian mean–variance segmentation,
  with scalar and optional Highway backends and online objects.
- Made `"highway"` the default backend spelling; the former capitalized
  spelling remains an alias.
- Validated model domains and penalties across all 1D methods; empty online
  appends are now no-ops and later penalties cannot silently change.
- Matched mean–variance active-set histories to the original `dust`
  implementation and made degenerate complete segmentations explicit.
- Corrected variance segmentation costs, known-size count normalization,
  and mixed Gaussian decay generation.
