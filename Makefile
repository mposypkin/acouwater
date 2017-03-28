all dep clean indent tests::
	cd alglib && make $@ && cd .. 
	cd sspemdd && make $@ && cd ..	
	cd ahw && make $@ && cd ..	

