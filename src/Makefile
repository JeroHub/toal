all: yaps3D_Pelagic.so GroupO2.cpp
	#Comment here prevents tab from converting to space
yaps3D_Pelagic.so: yaps3D_Pelagic.cpp
	Rscript --vanilla -e "TMB::compile('yaps3D_Pelagic.cpp','-O0 -g')"
GroupO2.so: GroupO2.cpp
	Rscript --vanilla -e "TMB::compile('GroupO2.cpp','-O0 -g')"

clean:
	rm -rf *o yaps3D_Pelagic.so GroupO2.so
