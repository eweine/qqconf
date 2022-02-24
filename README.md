# qqconf
An R Package to Put Simultaneous Testing Bands on QQ and PP Plots

To install our code from CRAN, run

```
install.packages("qqconf")
```

For an introductory tutorial, please see our [vignette](https://cloud.r-project.org/web/packages/qqconf/vignettes/qqconf_introduction.html).

Full documentation can be found in our [reference manual](https://cloud.r-project.org/web/packages/qqconf/qqconf.pdf) and additional methodological details can be found in [our paper](https://arxiv.org/abs/2111.15082).

## Acknowledgements

Thanks to Amit Moscovich and Boaz Nadler for their help implementing the FFT method for computing testing bands. Our C++ code is copied from [Amit's repository](https://github.com/mosco/crossing-probability), with a few small changes for R compatibility.
