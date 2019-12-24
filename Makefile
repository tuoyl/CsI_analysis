target:
	g++ data_reduction.cxx -I${HEADAS}/include -L${HEADAS}/lib -lcfitsio  -o data_reduction.o
	./data_reduction.o ../test_data/Y201802/20180201-0232/HXMT_20180201T00_HE-Evt_FFFFFF_V1_1K.FITS ../test_data/Y201802/20180201-0232/HXMT_20180201T00_HE-HV_FFFFFF_V1_1K.FITS  ../test_data/Y201802/20180201-0232/HXMT_20180201T00_Orbit_FFFFFF_V1_1K.FITS
