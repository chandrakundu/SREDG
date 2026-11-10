# Robust Euclidean Distance Geometry from Sparsely Corrupted Anchor-Based Measurements

MATLAB code for SREDG and for the phase-transition experiment (synthetic sensor localization) in the paper.

## Abstract

We study robust Euclidean distance geometry from an anchor-distance sketch in which anchor-anchor distances are exact, anchor-target distances contain sparse gross errors, and target-target distances are unavailable. We propose Structured Robust Euclidean Distance Geometry (SREDG), which combines robust principal component analysis with Nyström reconstruction through an anchor-centered Gram representation to recover the point configuration. Under suitable incoherence, sampling, rank, and sparsity conditions, we establish end-to-end recovery guarantees for the Gram matrix and for the point configuration up to an orthogonal transformation. Experiments on synthetic and real-world data show that SREDG achieves lower reconstruction error than robust Euclidean distance matrix optimization and robust matrix completion baselines in most tested settings, often by orders of magnitude. Metallic-glass experiments further show that SREDG typically preserves local atomic structure to numerical precision.

## Quick start

From this folder in MATLAB:

```matlab
run_demo                      % quick run
```

Change any setting in `setup.m` or `run_demo.m`, or in `cfg` directly. For example, to also run PREEEDM and RGNMR (clone them first, see below):

```matlab
cfg = setup();
cfg.methods = {'SREDG', 'PREEEDM', 'RGNMR'};
run_phase_transition(cfg);
```

Each run writes results in `results/run_<date>_<time>/`.

## PREEEDM and RGNMR (optional)

To compare with PREEEDM and RGNMR, clone the authors' repositories into `lib/Baselines/` at the commits we used. 

```bash
git clone https://github.com/ShenglongZhou/PREEEDM.git lib/Baselines/PREEEDM
git -C lib/Baselines/PREEEDM checkout 2dc9624
git clone https://github.com/eilon96/RGNMR.git lib/Baselines/RGNMR
git -C lib/Baselines/RGNMR checkout 66dfb8c
```

## Requirements

- MATLAB R2023a or newer, with the Statistics and Machine Learning Toolbox.
- Optional: the Parallel Computing Toolbox, for parallel workers (`cfg.workers`).
- Optional: the Wavelet Toolbox (a built-in fallback is used without it).

## Citation

Preprint and BibTeX: coming soon.

## License

MIT (`LICENSE`), except the code in `lib/RPCA`, `lib/PROPACK`, and the patch folders, which keep the terms of their authors.

## References

1. C. Kundu, A. Tasissa, and H. Cai, "Structured sampling for robust Euclidean distance geometry," in *Proc. 59th Annual Conference on Information Sciences and Systems (CISS)*, 2025, pp. 1–6. [Paper](https://doi.org/10.1109/CISS64860.2025.10944739)
1. H. Cai, J.-F. Cai, and K. Wei, "Accelerated alternating projections for robust principal component analysis," *Journal of Machine Learning Research*, vol. 20, no. 20, pp. 1–33, 2019. [Paper](https://jmlr.org/papers/v20/18-022.html)
1. S. Zhou, N. Xiu, and H.-D. Qi, "Robust Euclidean embedding via EDM optimization," *Mathematical Programming Computation*, vol. 12, no. 3, pp. 337–387, 2020. [Paper](https://doi.org/10.1007/s12532-019-00168-0)
1. E. Vaknin Laufer and B. Nadler, "RGNMR: A Gauss-Newton method for robust matrix completion with theoretical guarantees," in *Advances in Neural Information Processing Systems (NeurIPS)*, vol. 38, 2025. [Paper](https://papers.nips.cc/paper_files/paper/2025/hash/56bd21259e28ebdc4d7e1503733bf421-Abstract-Conference.html)
1. R. M. Larsen, "Lanczos bidiagonalization with partial reorthogonalization," DAIMI Report PB-357, Aarhus University, 1998. Software: [Code](https://github.com/rmlarsen/propack)
