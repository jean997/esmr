Install with vignettes

```
devtools::install_github("jean997/mrScan",  build_vignettes = TRUE)
browseVignettes("mrScan")
```

Note that the vignette requires the following packages which will not be 
installed automatically:

+ DiagrammeR (use `install.packages`)
+ reshape2 (use `install.packages`)
+ GRAPPLE (`devtools::install_github("jingshuw/grapple")`)
