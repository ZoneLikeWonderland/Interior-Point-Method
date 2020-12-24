ipm_v2: main_v2.cpp
	g++ main_v2.cpp -o ipm_v2 -O3 -march=native

test: main_v2.cpp
	g++ main_v2.cpp -o ipm_v2
	./ipm_v2  ../../data/data_matlab/brandy --detail

ipm: main.cpp
	g++ main.cpp -o ipm -O3 -march=native

run: ipm_v2
	./ipm_v2  ../../data/data_matlab/grow7 --detail