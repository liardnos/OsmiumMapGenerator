NAME = mapGen

all:
	g++ osmiumGen.cpp mat/mat.cpp -O4 -lz -lbz2 -lpthread -lexpat -Wall -Wextra -lm -lsfml-graphics -lsfml-window -lsfml-system

run: all
	./a.out ~/openStreetMap/pont-a-marq.osm

launch:
	./a.out ~/openStreetMap/pont-a-marq.osm

simplify:
	#osmium tags-filter france_metro_dom_com_nc.osm.pbf wr/boundary=administrative -o france_filtered.pbf --overwrite