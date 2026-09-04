rcpp:
	rm -f src/RcppExports.cpp R/RcppExports.R
	Rscript -e "Rcpp::compileAttributes()"

deps:
	Rscript -e "pak::local_install_dev_deps()"

install:
	Rscript -e "pak::local_install()"