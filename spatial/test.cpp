
#if 1

// #include "UniTreeZone.hpp"
#include "spacialHashZone.hpp"
#include <SFML/Graphics.hpp>
#include <SFML/Graphics.hpp>
#include <unistd.h>

class Point : public Zone<float, 2> {
public:
    Point(Zone<float, 2> zone) : 
        Zone<float, 2>(zone)
    {
        uint32_t *color = (uint32_t *)&_color;
        *color = rand();
        _color[3] = 255;
        static uint count = 0;
        _count = count++;
    }
    unsigned char _color[4];
    uint _count;
};

int main() {
    std::cout << "INIT TREE" << std::endl;
    UniTreeZone<Point, Zone<float, 2>, 2> uniTreeZone = Zone<float, 2>{NDVector<float, 2>{1024/2, 1024/2}, NDVector<float, 2>{1024/2, 1024/2}};
    std::vector<Point> points;
    uint triCount = 256;
    points.reserve(triCount);
    for (uint i = 0; i < triCount; i++) {
        NDVector<float, 2> pos = {(float)80+(rand()%(1024-80*2)), (float)80+(rand()%(1024-80*2))};
        NDVector<float, 2> size = {(float)(rand()%80), (float)(rand()%80)};
        points.emplace_back(Zone<float, 2>{pos, size});
        uniTreeZone.addData(&points[i]);
    }

    NDVector<uint, 2> screenSize = {1024, 1024};
    sf::RenderWindow window(sf::VideoMode(screenSize[0], screenSize[1]), "EVOLVESFML");
    window.setFramerateLimit(60);
    window.setPosition(sf::Vector2i(0, 0));


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed){
                window.close(); 
            }
        }

        if (1) {
            uniTreeZone.forEachReverse(
                [&window](UniTreeZone<Point, Zone<float, 2>, 2> const &treeNode, uint depth)->void {
                    (void)depth;   
                    if (1) {
                        uint outlineThickness = 2;
                        sf::RectangleShape rect({treeNode._zone.size[0]*2-outlineThickness*2, treeNode._zone.size[1]*2-outlineThickness*2});
                        rect.setFillColor({0, 0, 0, 0});
                        rect.setOutlineColor({50, 50, 50, 255});
                        rect.setOutlineThickness(2);
                        rect.setPosition(treeNode._zone.pos[0]-treeNode._zone.size[0]+outlineThickness, treeNode._zone.pos[1]-treeNode._zone.size[1]+outlineThickness);
                        window.draw(rect);
                    }
                }
            );
        }

        uint count = 0;
        if (1) {
            uniTreeZone.forEachReverse(
                [&window, &count](UniTreeZone<Point, Zone<float, 2>, 2> const &treeNode, Point const &data, uint depth)->void {
                    (void)treeNode;
                    (void)depth;
                    if (0) {
                        sf::RectangleShape rect({treeNode._zone.size[0]*2, treeNode._zone.size[1]*2});
                        rect.setFillColor({255/depth, 255/depth, 255/depth, 255});
                        rect.setPosition(treeNode._zone.pos[0]-treeNode._zone.size[0], treeNode._zone.pos[1]-treeNode._zone.size[1]);
                        window.draw(rect);
                    }

                    uint outlineThickness = 2;
                    sf::RectangleShape rect({data.size[0]*2-outlineThickness*2, data.size[1]*2-outlineThickness*2});
                    rect.setFillColor({0, 0, 0, 0});
                    rect.setOutlineColor({data._color[0], data._color[1], data._color[2], 255/2});
                    rect.setOutlineThickness(2);
                    rect.setPosition(data.pos[0]-data.size[0]+outlineThickness, data.pos[1]-data.size[1]+outlineThickness);
                    window.draw(rect);

                    static sf::Text *text = 0;
                    if (!text) {
                        text = new sf::Text();
                        static sf::Font font;
                        if (!font.loadFromFile("../font/font.TTF")){
                            std::cerr << "cannot load font" << std::endl;
                            throw "cannot load font";
                        }
                        text->setFont(font);
                        text->setCharacterSize(20);
                        text->setFillColor({255, 255, 255, 128});
                    }
                    text->setString(std::to_string(data._count));
                    text->setPosition(data.pos[0]-data.size[0]+outlineThickness, data.pos[1]-data.size[1]+outlineThickness);
                    window.draw(*text);

                    count++;
                }
            );
        }
        NDVector<float, 2> rectSize = {30, 50};
        
        Zone<float, 2> rect(NDVector<float, 2>{0, 0}, rectSize);

        auto mousePos = sf::Mouse::getPosition(window);
        rect.pos = {(float)mousePos.x, (float)mousePos.y};

        uint outlineThickness = 2;
        sf::RectangleShape sfRect({rect.size[0]*2-outlineThickness*2, rect.size[1]*2-outlineThickness*2});
        sfRect.setFillColor({0, 0, 0, 0});
        sfRect.setOutlineColor({255, 255, 255, 255});
        sfRect.setOutlineThickness(2);
        sfRect.setPosition(rect.pos[0]-rect.size[0]+outlineThickness, rect.pos[1]-rect.size[1]+outlineThickness);
        window.draw(sfRect);
        std::cout << "loop" << std::endl;

        std::shared_ptr<std::vector<Point*>> res = uniTreeZone.getColides(rect);
        for (Point* &p : *res) {
            Point &rect = *p;
            uint outlineThickness = 2;
            sf::RectangleShape sfRect({rect.size[0]*2-outlineThickness*2, rect.size[1]*2-outlineThickness*2});
            sfRect.setFillColor({0, 0, 0, 0});
            sfRect.setOutlineColor({255, 255, 255, 255});
            sfRect.setOutlineThickness(2);
            sfRect.setPosition(rect.pos[0]-rect.size[0]+outlineThickness, rect.pos[1]-rect.size[1]+outlineThickness);
            window.draw(sfRect);
        }

        window.display();
        window.clear();
    }

    return 0;
}
#endif