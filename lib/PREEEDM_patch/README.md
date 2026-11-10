# PREEEDM patch

Our changes to [PREEEDM](https://github.com/ShenglongZhou/PREEEDM), commit `2dc9624`. Clone PREEEDM into `lib/Baselines/PREEEDM` (see the main README). The original files are not changed.

PREEEDM treats a zero distance as missing. In our experiments, negative observed squared distances are set to zero and kept as observed, so PREEEDM uses the same entries as SREDG. 
