# Drake plan

Here's a thing I've been trying these days. It's a package called Drake that makes file dependencies and data structure easier. (Here is a nice intro to why it's nice and makes stuff reproducible and easier to debug)[https://books.ropensci.org/drake/]. These mostly duplicate the jm_Allen.R work, but it should save tons of time by cacheing biomart queries, for example.

I followed this guide:
https://milesmcbain.xyz/posts/the-drake-post/

More details provided in:
https://github.com/milesmcbain/dflow

To run the thing, just open the \_drake.R file, run it and then do: `r_make()`. 

Usefult things:
vis_drake_graph(plan) # to see what fails, what's up to date and what not