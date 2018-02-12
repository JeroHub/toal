all: yaps3D_Pelagic.so
	#Comment here prevents tab from converting to space
yaps3D_Pelagic.so: yaps3D_Pelagic.cpp
	Rscript --vanilla -e "TMB::compile('yaps3D_Pelagic.cpp','-O0 -g')"

clean:
	rm -rf *o yaps3D_Pelagic
