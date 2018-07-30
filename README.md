# CTPCP 
## Codes
 
- You can tune parameters for `ctpcp.m` by `expRun.m`.
- Other helping funtions can be found in `methods/`
- Codes listed below for tensor operators are borrowed from [LibADMM](https://github.com/canyilu/LibADMM):

  ```
proximal_operator/
tensor_tools/
methods/trpca_tnn.m
  ```

## Documentation
- `admm_ctpcp/` explains the algorithm for solving CTPCP by ADMM in details.
- `tex/` contains original LaTeX files for the draft.

## References:
- Simulation settings：[Tensor Robust Principal Component Analysis: Exact Recovery of CorruptedLow-Rank Tensors via Convex Optimization](http://www.cis.pku.edu.cn/faculty/vision/zlin/Publications/2016-CVPR-TRPCA.pdf)

- Ituition on solving CTPCP by ADMM: [Dynamic MR Image Reconstruction–Separation FromUndersampled (k, t)-Space via Low-Rank Plus Sparse Prior](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6808502)

- Guassion Sampling: page 6 of [Compressive Principal Component Pursuit](https://arxiv.org/abs/1202.4596)