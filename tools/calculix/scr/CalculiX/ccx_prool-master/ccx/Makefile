# Prool's makefile
all:
	cd ARPACK;make home=`pwd` lib
	cd SPOOLES.2.2;make lib
	cd CalculiX/ccx_2.12/src;make
	cp CalculiX/ccx_2.12/src/ccx_2.12_MT ccx
clean:
	cd ARPACK;make clean
	cd SPOOLES.2.2;make clean
	cd CalculiX/ccx_2.12/src;make clean
	rm CalculiX/ccx_2.12/*.o
	rm ccx
