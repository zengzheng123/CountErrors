all:
	g++ -o3 -I./bamtools-master/include/ -L./bamtools-master/lib/ CountErrors.cpp -lbamtools -lz -o CountErrors

clean:
	rm -f CountErrors