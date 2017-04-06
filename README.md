# TCNeighborsFinder

INSTALLATION
------------

- install Shapely (https://pypi.python.org/pypi/Shapely)

then

- compile geos library
	wget http://download.osgeo.org/geos/geos-3.4.2.tar.bz2
	tar xjf geos-3.4.2.tar.bz2
	cd geos-3.4.2
	./configure
	make
	sudo make install
or
- install bunutils package 
	sudo apt-get install binutils

DESCRIPTION OF THE ALGORITHM
----------------------------
STEP 1:

Extraction of configuration codewords
------------------------------------- 
codeword = 7 bits where 1st bit is central wafer type abd the rest are 6 neighbor wafer's types in the clockwise direction
wafertype_small = bit 1
wafertype_large = bit 0

1a) Extracted wafer information on waferIDs and their (row, col) position only on Layer 1 => wafer_info_l1
1b) Extracted waferIDs and their types => wafer_types_map_l1
1c) Extracted "interesting" central wafers with their configurations on Layer 1 
    - the ones that have different configurations, because we don't want all the wafers
    - number of "interesting" central wafers corresponds to number of different configurations on Layer 1
1d) Since some of the configurations are the same in FH as they are in EE (Layer 1), list of interesting wafers of FH is shortened (ejecting the doubles that are already in EE)
1e) Lists of interesting wafers with configurations for EE and FH are divided by type (small or large lists)
    - just for extracting TC neighbors part by part

STEP 2:

TC neighbor Finder Algo
-----------------------
main function => TCNeighborFinder(...)

Foreach interesting wafer 
    1) extract its wafer ring
    2) extract all cellIDs with their position (x,y) from the  current wafer and the ring
    3) convert all cell's (x,y) coordinates to 60 degrees mapping coordinates (referent cell is the one with (minx, miny))
    4) extract all TCs with cellIDs from the current wafer and its ring => waferID (key), TC->cellID map (value)
    5) produce maps cellID->TC and cellID->waferID 
    6) extract neighbor cells for each cellID you extracted (from the current wafer and the ring)
        - if current cell's wafer is the same as the neighbor cell's waferID (you are at the border between same wafer type)
            - extract neighbors using mapping coordinates
        - else (you are at the border between two different wafer types)
            - extract neighbors using polygon intersection
    7) convert current wafer ring waferIDs to relative (this is important for the neighbor tuple generation (TC, relative_waferID))
    8) extract neighbor TCs
        - examine each cellID with neighbor cells 
            - if current cell's TC is different than the neighbor cell's TC then neighbor cell's Tc is a neghbor TC of the current cell's TC
    9) update the key of the Tc neighbors map with the current interesting wafer's configuration          
