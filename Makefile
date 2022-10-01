NAME = mapGen

run:
	g++ osmiumGen.cpp mat/mat.cpp -O4 -lz -lbz2 -lpthread -lexpat -Wall -Wextra -lm -lsfml-graphics -lsfml-window -lsfml-system
	./a.out ~/openStreetMap/pont-a-marq.osm

reSFML:
	rm -rf buildSFML
	mkdir buildSFML; cd buildSFML && cmake .. && make -j6
	cp -r assets buildSFML/
	cp boatcpy.bin buildSFML/

SFML:
	cd buildSFML && rm -f $(NAME)
	cd buildSFML && make -j6 || clear && make

runSFML: SFML


SFMLValgrind: SFML
	##########################################################
	#                       VALGRIND                         #
	##########################################################
	cd buildSFML && valgrind --track-origins=yes ./$(NAME)

	cd buildSFML && ./$(NAME)