# Circos plot

First, download the program:
```sh
wget http://circos.ca/distribution/circos-0.69-9.tgz
tar xvf circos-0.69-9.tgz
rm circos-0.69-9.tgz
``` 

Second, make sure circos run. This should be easy but tedious. Mostly, you either run `bin/.list.modules` and bin/test.modules, or you run the example (executing `example/run`) and see which modules are missing. Make sure you have Perl installed. To install modules in perl, do:

```sh
perl -MCPAN -e shell
% install Config::General
% install whatever::modules
```

Now for the actual plot: unzip the zip in here (`raul_plot.zip`) into `bin/`. Run it from a terminal there with something like:

```sh
./circos -conf start.conf
```
If you want to change features refer to the documentation (http://circos.ca/documentation/) or ask me.