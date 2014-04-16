Motifs found
========================================================

Simple code to see what motifs were found.


```r
motifs = read.csv("../motifsfilesformatted.txt", sep = "\t", h = F)
head(motifs)
```

```
##             V1 V2
## 1 WGD2ANC00031 19
## 2 WGD2ANC00031 16
## 3 WGD2ANC00182 21
## 4 WGD2ANC00182 10
## 5 WGD2ANC00182 11
## 6 WGD2ANC00182  9
```

```r
colnames(motifs) = c("family", "motif.size")
head(motifs)
```

```
##         family motif.size
## 1 WGD2ANC00031         19
## 2 WGD2ANC00031         16
## 3 WGD2ANC00182         21
## 4 WGD2ANC00182         10
## 5 WGD2ANC00182         11
## 6 WGD2ANC00182          9
```


How many motifs are there? What distribution? What families are represented?


```r
library(ggplot2)
library(scales)
g = ggplot(motifs, aes(x = motif.size))
g + geom_histogram(aes(y = ..density..), fill = "darkgreen") + geom_density() + 
    labs(x = "Motif Size", y = "Density") + scale_y_continuous(labels = percent)
```

```
## stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-21.png) 

```r
g + geom_histogram(fill = "darkblue") + labs(x = "Motif Size", y = "Counts")
```

```
## stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-22.png) 



```r
library(reshape)
```

```
## Loading required package: plyr
```

```
## Warning: package 'plyr' was built under R version 3.0.2
```

```
## 
## Attaching package: 'reshape'
## 
## L'objet suivant est masqué from 'package:plyr':
## 
##     rename, round_any
```

```r
cast(motifs, formula = ~family, fun.aggregate = mean)
```

```
## Using motif.size as value column.  Use the value argument to cast to override this choice
```

```
##   value WGD2ANC00031 WGD2ANC00182 WGD2ANC00311 WGD2ANC00387 WGD2ANC00451
## 1 (all)         17.5        12.75            9            7           24
##   WGD2ANC00488 WGD2ANC00492 WGD2ANC00505 WGD2ANC00530 WGD2ANC00581
## 1           11           19           10           12            7
##   WGD2ANC00630 WGD2ANC00651 WGD2ANC00656 WGD2ANC00687 WGD2ANC00726
## 1            9           35           16            8           10
##   WGD2ANC00907 WGD2ANC00928 WGD2ANC00949 WGD2ANC00958 WGD2ANC01001
## 1            7            8           18            8           16
##   WGD2ANC01066 WGD2ANC01102 WGD2ANC01294 WGD2ANC01412 WGD2ANC01495
## 1           23           12          8.5            8            8
##   WGD2ANC01530 WGD2ANC01532 WGD2ANC01563 WGD2ANC01568 WGD2ANC01574
## 1           18           12           24           11            8
##   WGD2ANC01595 WGD2ANC01669 WGD2ANC01715 WGD2ANC01761 WGD2ANC01792
## 1           14            7           19            8            7
##   WGD2ANC01916 WGD2ANC01925 WGD2ANC01930 WGD2ANC01993 WGD2ANC02135
## 1           11            8           15         16.5            8
##   WGD2ANC02164 WGD2ANC02226 WGD2ANC02362 WGD2ANC02383 WGD2ANC02424
## 1           15           12           11           46           16
##   WGD2ANC02466 WGD2ANC02555 WGD2ANC02636 WGD2ANC02697 WGD2ANC02721
## 1           13            9            9           18           12
##   WGD2ANC02833 WGD2ANC02862 WGD2ANC02942 WGD2ANC02960 WGD2ANC03003
## 1           11            8           21           20            7
##   WGD2ANC03049 WGD2ANC03059 WGD2ANC03088 WGD2ANC03116 WGD2ANC03129
## 1            9            8           22            9            8
##   WGD2ANC03139 WGD2ANC03145 WGD2ANC03316 WGD2ANC03407 WGD2ANC03436
## 1           12            7            7            8            8
##   WGD2ANC03453 WGD2ANC03483 WGD2ANC03603 WGD2ANC03636 WGD2ANC03671
## 1           11           13           12           11           14
##   WGD2ANC03690 WGD2ANC03694 WGD2ANC03699 WGD2ANC03869 WGD2ANC03965
## 1            7           15            7            7           11
##   WGD2ANC04128 WGD2ANC04174 WGD2ANC04187 WGD2ANC04245 WGD2ANC04264
## 1            7           12            8           16           10
##   WGD2ANC04356 WGD2ANC04473 WGD2ANC04677 WGD2ANC04783 WGD2ANC04815
## 1           15           15            7           10            8
##   WGD2ANC04917 WGD2ANC04935 WGD2ANC04946 WGD2ANC04989 WGD2ANC05162
## 1           12           13            7           16            7
##   WGD2ANC05188 WGD2ANC05217 WGD2ANC05243 WGD2ANC05283 WGD2ANC05287
## 1           10           12            8            7           28
##   WGD2ANC05311 WGD2ANC05446 WGD2ANC05453 WGD2ANC05475 WGD2ANC05510
## 1            8            7          7.5           12          9.5
##   WGD2ANC05562 WGD2ANC05614 WGD2ANC05654 WGD2ANC05655 WGD2ANC05657
## 1           11            9           11           22           12
##   WGD2ANC05674 WGD2ANC05718 WGD2ANC05721
## 1           11            7          7.5
```

```r
quantile(motifs$motif.size, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
```

```
##   5%  10%  25%  50%  75%  90%  95% 
##  7.0  7.0  8.0 10.0 14.5 21.0 24.0
```

```r
summary(motifs)
```

```
##           family      motif.size  
##  WGD2ANC04473:  5   Min.   : 7.0  
##  WGD2ANC00182:  4   1st Qu.: 8.0  
##  WGD2ANC03088:  3   Median :10.0  
##  WGD2ANC00031:  2   Mean   :12.2  
##  WGD2ANC00505:  2   3rd Qu.:14.5  
##  WGD2ANC01294:  2   Max.   :46.0  
##  (Other)     :113
```

```r

# Show the two families with maximum occurrences

motifs[motifs$family == "WGD2ANC04473", ]
```

```
##           family motif.size
## 96  WGD2ANC04473         22
## 97  WGD2ANC04473          8
## 98  WGD2ANC04473         13
## 99  WGD2ANC04473         25
## 100 WGD2ANC04473          7
```

```r
motifs[motifs$family == "WGD2ANC00182", ]
```

```
##         family motif.size
## 3 WGD2ANC00182         21
## 4 WGD2ANC00182         10
## 5 WGD2ANC00182         11
## 6 WGD2ANC00182          9
```

```r


# Show maximum and minimum
motifs[motifs$motif.size == max(motifs$motif.size), ]
```

```
##          family motif.size
## 52 WGD2ANC02383         46
```

```r
motifs[motifs$motif.size == min(motifs$motif.size), ]
```

```
##           family motif.size
## 8   WGD2ANC00387          7
## 15  WGD2ANC00581          7
## 21  WGD2ANC00907          7
## 38  WGD2ANC01669          7
## 41  WGD2ANC01792          7
## 44  WGD2ANC01925          7
## 56  WGD2ANC02555          7
## 64  WGD2ANC03003          7
## 68  WGD2ANC03088          7
## 74  WGD2ANC03145          7
## 75  WGD2ANC03316          7
## 82  WGD2ANC03636          7
## 85  WGD2ANC03690          7
## 87  WGD2ANC03699          7
## 88  WGD2ANC03869          7
## 90  WGD2ANC04128          7
## 100 WGD2ANC04473          7
## 101 WGD2ANC04677          7
## 106 WGD2ANC04946          7
## 108 WGD2ANC05162          7
## 113 WGD2ANC05283          7
## 116 WGD2ANC05446          7
## 118 WGD2ANC05453          7
## 129 WGD2ANC05718          7
## 130 WGD2ANC05721          7
```


