


test:
	clear
	g++ -Wall -Wextra -g *.cpp -lsfml-graphics -lsfml-window -lsfml-system -pthread
	./a.out

valgrind:
	clear
	g++ -Wall -Wextra -g *.cpp -lsfml-graphics -lsfml-window -lsfml-system -pthread
	valgrind --track-origins=yes ./a.out