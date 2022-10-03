NAME = mapGen

SRC = osmiumGen.cpp mat/mat.cpp ByteObject.cpp

FLAGS = -lz -lbz2 -lpthread -lexpat -Wall -Wextra -lm -lsfml-graphics -lsfml-window -lsfml-system

all:
	g++ $(SRC) $(FLAGS) -g -o $(NAME)

run: all
	./$(NAME) ~/openStreetMap/pont-a-marq.osm

launch:
	./$(NAME) ~/openStreetMap/pont-a-marq.osm

callgrindOg:
	rm -f callgrind.*
	clear
	g++ -g -o $(NAME) $(SRC) $(FLAGS)
	-valgrind --tool=callgrind ./$(NAME) Pont-à-Marcq.osm.pbf #&> valgrind_log
	-kcachegrind callgrind.*

simplify:
	#osmium tags-filter france_metro_dom_com_nc.osm.pbf wr/boundary=administrative -o france_filtered.pbf --overwrite
