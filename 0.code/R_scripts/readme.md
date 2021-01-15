# Drake plan

Here's a thing I've been trying these days. It's a package called Drake that makes file dependencies and data structure easier. (Here is a nice intro to why it's nice and makes stuff reproducible and easier to debug)[https://books.ropensci.org/drake/]. These mostly duplicate the jm_Allen.R work, but it saves tons of time by cacheing biomart queries, for example, and makes the code a bit less redundant. I made sure everything is the same, except for a weird bug I have when exporting to xlsx that makes me unable to export data, so I commented that out (if that's ok).

I followed this guide:
https://milesmcbain.xyz/posts/the-drake-post/

Usefult things:
- To run the whole thing, just open the \_drake.R file, run it and then do: `r_make()`. 
- vis_drake_graph(the_plan) # to see what fails, what's up to date and what not
- clean() wipes out the cache in case you want to run each target from the beginning

Folder structure:
- Rscripts holds the drake structure. More details about how the Rscripts folder is organized are provided in:
https://github.com/milesmcbain/dflow
- Output holds the output of the scripts
- file_dependencies ensures the files aren't downloaded from some drive and always available

Overall, apart from the trouble of setting it up, I think it has already saved me hours of biomart-ing and dataset wrangling every time I want to run something. Let me know what you think.