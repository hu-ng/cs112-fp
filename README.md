# Final Project for CS112

## Description:

This project uses *genetic matching* to extend the research of Lu et al. ([2013](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2354965)) into the low-risk anomaly in the area of stock exchange, which is when low-risk stocks have higher returns than high-risk stocks, contrary to finance fundamentals. The original paper uses *coarsened exact matching* to match the least risky stocks with the most risky ones on the basis of industry group, company size, and trading volume, and found that there is a muted effect of the anomaly. Using genetic matching, we found that the effect of the anomaly is even more muted than the original paper, and the extension results are more reliable due to better balance achieved.

## Running Instructions:
Have RStudio installed, clone the files and run `replication_genmatch_new.R`

```
$ git clone https://github.com/hu-ng/cs112-fp.git
$ cd cs112-fp
$ replication_genmatch_new.R
```
