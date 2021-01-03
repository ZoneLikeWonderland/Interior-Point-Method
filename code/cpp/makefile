ipm: main.cpp
	g++ main.cpp -o ipm -O3 -march=native

test: main.cpp
	g++ main.cpp -o ipm
	./ipm  ../../data/data_matlab/lpi_klein1 --detail

run: ipm
	./ipm  ../../data/data_matlab/capri --detail