This is very raw instruction on how to install ATTRACT on Linux (does not work on Microsoft Windows nor on Mac):

1. Obtain gfortran, version 13.3.0 
2. Obtain g++, version 11.5.0  
3. Make and activate attract environment with Python 2.7, numpy, scipy and pdb2pqr 
(normally, pyenv/conda/mamba can deal with it)
4. git clone  --recursive https://github.com/sjdv1982/attract.git
5. edit attract/bin/Makefile:
add ` FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer" ` to the end of each line with `$(MAKE)` (but not $(MAKE) clean). 
This is what you want to see:
```

em: as
	cd $(GPU_DIR)/emATTRACT/bin && $(MAKE) TARGET=$(TARGET) FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer"
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/emATTRACT/bin/emATTRACT 

gpuATTRACTcheck: as
	cd $(GPU_DIR)/gpuATTRACTcheck/bin && $(MAKE) TARGET=$(TARGET) FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer"
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/gpuATTRACTcheck/bin/gpuATTRACTcheck

mc: as
	cd $(GPU_DIR)/mcATTRACT/bin && $(MAKE) TARGET=$(TARGET) FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer"
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/mcATTRACT/bin/mcATTRACT
	
sc: as
	cd $(GPU_DIR)/scATTRACT/bin && $(MAKE) TARGET=$(TARGET) FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer"
	cd $(BINDIR) && ln -sf $(GPU_DIR)/scATTRACT/bin/scATTRACT

as: check_cuda_version
	cd $(GPU_DIR)/AttractServer/lib && $(MAKE) TARGET=$(TARGET) FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer"
	cd $(BINDIR) && ln -sf $(GPU_DIR2)/AttractServer/lib/libAttractServer.so

attract:
	cd $(ATTRACTDIR) && $(MAKE) -f _Makefile all FFLAGS="-std=f90 -g -O3 -fno-automatic -ffast-math -fcray-pointer"

```
6. cd attract/bin 
7. make attract

If there are multiple versions of gfortran and g++ on your system, you can point to the required ones by adding theses to Makefile, line 6:
CXX = /path/gcc@11
FC = /path/gfortran@13

If you need to restart installation - run `make clean` first 

8. Export the ATTRACT environment variables. In your .bashrc file, add the following:
 ```
  export ATTRACTDIR=/path/attract/bin (i.e. wherever you installed attract)
  export ATTRACTTOOLS=$ATTRACTDIR/../tools
  export ATTRACTGUI=$ATTRACTDIR/../gui
  export LD_LIBRARY_PATH=$ATTRACTDIR:$LD_LIBRARY_PATH
```
then: 
 ` source ~/.bashrc`
then activate attract environment again

9. Test your installation:
- $ATTRACTDIR/attract" should result in:
  Too few arguments
  usage: $path/attract structures.dat parameterfile receptor.pdb [ligand.pdb] [options]
- "python -c 'import numpy, scipy' " should give no error
- "python -c 'import pdb2pqr' " should give no error