Between scripts 03 and 05, the BSGenome packages need to be assembled. 

This is to be done on the command line.

Navigate into each folder in ‘output/BSGenome_pkgs’ and call

R CMD build BSGenome.Cxxxx.Wormbase.XXX

When build is complete, install with

R CMD install BSGenome.Cxxxx.Wormbase.XXX_1.0.tar.gz

They should now be able to be loaded within R 