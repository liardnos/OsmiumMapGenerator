NAME = mapGen

SRC = osmiumGen.cpp mat/mat.cpp ByteObject.cpp

FLAGS = -O4 -g -I$(CONDA_PREFIX)/include/python3.10 -L$(CONDA_PREFIX)/lib -lpython3.10 -std=gnu++2a -Wall -Wextra -lcurlpp -lcurl -lz -lbz2 -lpthread -lexpat -lm -lsfml-graphics -lsfml-window -lsfml-system



all:
	g++ $(SRC) $(FLAGS) -o $(NAME)

run: all
	./$(NAME) ~/openStreetMap/pont-a-marq.osm

launch:
	./$(NAME) ~/openStreetMap/pont-a-marq.osm

valgrind:
	clear
	valgrind --leak-check=full ./mapGen Pont-à-Marcq.osm.pbf

callgrindOg:
	rm -f callgrind.*
	clear
	g++ -g -o $(NAME) $(SRC) $(FLAGS)
	-valgrind --tool=callgrind ./$(NAME) Pont-à-Marcq.osm.pbf #&> valgrind_log
	-kcachegrind callgrind.*

simplify:
	./mapGen name=Paris ~/openStreetMap/france_filtered.pbf ~/openStreetMap/france_metro_dom_com_nc.osm.pbf
	#osmium tags-filter france_metro_dom_com_nc.osm.pbf wr/boundary=administrative -o france_filtered.pbf --overwrite


#export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
