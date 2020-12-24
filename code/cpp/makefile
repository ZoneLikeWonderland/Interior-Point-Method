VERSION=v3

ipm_$(VERSION): main_$(VERSION).cpp
	g++ main_$(VERSION).cpp -o ipm_$(VERSION) -O3 -march=native

test: main_$(VERSION).cpp
	g++ main_$(VERSION).cpp -o ipm_$(VERSION)
	./ipm_$(VERSION)  ../../data/data_matlab/lpi_klein1 --detail

ipm: main.cpp
	g++ main.cpp -o ipm -O3 -march=native

run: ipm_$(VERSION)
	./ipm_$(VERSION)  ../../data/data_matlab/stocfor1 --detail