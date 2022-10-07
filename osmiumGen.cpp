/*

  EXAMPLE osmium_area_test

  Create multipolygons from OSM data and dump them to stdout in one of two
  formats: WKT or using the built-in Dump format.

  DEMONSTRATES USE OF:
  * file input
  * location indexes and the NodeLocationsForWays handler
  * the MultipolygonManager and Assembler to assemble areas (multipolygons)
  * your own handler that works with areas (multipolygons)
  * the WKTFactory to write geometries in WKT format
  * the Dump handler
  * the DynamicHandler

  SIMPLER EXAMPLES you might want to understand first:
  * osmium_read
  * osmium_count
  * osmium_debug
  * osmium_amenity_list

  LICENSE
  The code in this example file is released into the Public Domain.

*/

#include <cstring>  // for std::strcmp
#include <iostream> // for std::cout, std::cerr

// For assembling multipolygons
#include <osmium/area/assembler.hpp>
#include <osmium/area/multipolygon_manager.hpp>

// For the DynamicHandler class
#include <osmium/dynamic_handler.hpp>

// For the WKT factory
#include <osmium/geom/wkt.hpp>

// For the Dump handler
#include <osmium/handler/dump.hpp>

// For the NodeLocationForWays handler
#include <osmium/handler/node_locations_for_ways.hpp>

// Allow any format of input files (XML, PBF, ...)
#include <osmium/io/any_input.hpp>

// For osmium::apply()
#include <osmium/visitor.hpp>

// For the location index. There are different types of indexes available.
// This will work for all input files keeping the index in memory.
#include <osmium/index/map/flex_mem.hpp>

#include <SFML/Graphics.hpp>

#include <deque>
#include <vector>

#include <float.h>

#include <fstream>

#include "ByteObject.hpp"

#include "utils.hpp"

#include "mat/mat.hpp"

// #include "spatial/UniTreeZone.hpp"
#include "spatial/UniTree.hpp"

// The type of index used. This must match the include file above
using index_type = osmium::index::map::FlexMem<osmium::unsigned_object_id_type, osmium::Location>;

// The location handler always depends on the index type
using location_handler_type = osmium::handler::NodeLocationsForWays<index_type>;

class StreeShape {
public:
    // Segmentd _seg;
    sf::Color _color;
    std::vector<Vector2d> _points;
    std::vector<std::string> _labels;
    bool _drawed = false;

    friend ByteObject &operator<<(ByteObject &obj, StreeShape const &shape) {
        std::vector<Vector2f> vec;
        vec.reserve(shape._points.size());
        for (Vector2d const &p : shape._points)
            vec.push_back(p.cast<float>());

        obj << shape._color << vec << shape._labels;
        return obj;
    }

    friend ByteObject &operator>>(ByteObject &obj, StreeShape &shape) {
        std::vector<Vector2f> vec;
        obj >> shape._color >> vec >> shape._labels;
        shape._points.clear();
        shape._points.reserve(vec.size());
        for (Vector2f const &p : vec)
            shape._points.push_back(p.cast<double>());
        return obj;
    }
};


class Match {
public:
    char const *_key;
    char const *_tag;
    Color const _color;
};

std::vector<Match> g_colorVector = {
    {"wall", "no",  {0.25, 0.25, 0.25, 1.}},
    {"amenity", "parking_space", {0.25, 0.25, 0.25, 1.}},

    {"amenity", "parking", {0.5, 0.25, 0.25, 1.}},
    {"amenity", "parking", {0.5, 0.25, 0.25, 1.}},
    {"highway", "service", {0.5, 0.25, 0.25, 1.}},
    
    {"amenity", "school", {0.5, 0.5, 0.5, 1}},

    {"leisure", "track", {0.2, 0.2, 0.2, 1.}},

    {"landuse", "residential", {1., 0., 0., 1.}},
    
    {"water", "basin", {0., 0., 1., 1.}},
    {"natural", "water", {0., 0., 1., 1.}},
    {"leisure", "swimming_pool", {0., 0., 1., 1.}},

    {"leisure", "garden", {0., 1., 0., 1.}},
    {"leisure", "golf_course", {0., 1., 0., 1.}},
    {"leisure", "park", {0., 1., 0., 1.}},
    {"landuse", "grass", {0., 1., 0., 1.}},
    {"landuse", "village_green", {0., 1., 0., 1.}},
    
    {"landuse", "meadow", {0., 0.75, 0., 1.}},
    {"landuse", "vineyard", {0., 0.75, 0., 1.}},

    {"landuse", "forest", {0.1, 0.5, 0.1, 1.}},

    {"surface", "asphalt", {0.2, 0.2, 0.2, 1.}},

    {"social_facility", "group_home", {0.5, 0.5, 0.5, 1}},
    {"building", 0, {0.5, 0.5, 0.5, 1}},

};

// This handler writes all area geometries out in WKT (Well Known Text) format.
class WKTDump : public osmium::handler::Handler {

    // This factory is used to create a geometry in WKT format from OSM
    // objects. The template parameter is empty here, because we output WGS84
    // coordinates, but could be used for a projection.
    osmium::geom::WKTFactory<> m_factory;

public:

    WKTDump(std::vector<StreeShape> &segs, Segmentd &boundingBox) :
        osmium::handler::Handler(),
        _segs(segs),
        _boundingBox(boundingBox)
    {}

    // This callback is called by osmium::apply for each area in the data.
    void area(const osmium::Area& area) {
        try {
            StreeShape seg;
            for (auto &tag : area.tags()) {
                std::string label = std::string(tag.key()) + "=" + tag.value();
                seg._labels.push_back(label);
            }
            for (const auto& item : area) {
                if (item.type() == osmium::item_type::outer_ring) {
                    bool find = false;
                    Color color = {1, 1, 1, 1};
                    for (auto const &match : g_colorVector) {
                        if (match._tag == 0) {
                            if (area.tags().has_key(match._key)) {
                                color = match._color;
                                find = true;
                                break;
                            }
                        } else if (area.tags().has_tag(match._key, match._tag)) {
                            color = match._color;
                            find = true;
                            break;
                        }
                    }

                    seg._color = {(uint8_t)(color.r*255), (uint8_t)(color.g*255), (uint8_t)(color.b*255), (uint8_t)(color.a*255)};
                    // if (!find) {
                        // std::cout << std::endl;
                        // for (auto &tag : area.tags()) {
                        //     std::cout << tag << std::endl;
                        // }
                    // }

                    auto& ring = static_cast<const osmium::OuterRing&>(item);

                    for (auto const &p : ring) {
                        Vector2d pos = {p.lon(), -p.lat()};

                        _boundingBox._p[0] = std::min(pos[0], _boundingBox._p[0]);
                        _boundingBox._p[1] = std::min(pos[1], _boundingBox._p[1]);
                        _boundingBox._d[0] = std::max(pos[0], _boundingBox._d[0]);
                        _boundingBox._d[1] = std::max(pos[1], _boundingBox._d[1]);

                        seg._points.push_back(pos);
                    }
                    _segs.push_back(seg);

                    if (!(_segs.size() % 1024)) {
                        printf("\r%li ", _segs.size());
                        fflush(stdout);
                    }
                }
            }


        } catch (const osmium::geometry_error& e) {
            std::cout << "GEOMETRY ERROR: " << e.what() << "\n";
        }
    }

    std::vector<StreeShape> &_segs;
    Segmentd &_boundingBox;

};


class WKTAreaFinder : public osmium::handler::Handler {

    osmium::geom::WKTFactory<> m_factory;

public:

    WKTAreaFinder(Segmentd &_boundingBox) :
        osmium::handler::Handler(),
        _boundingBox(_boundingBox)
    {}

    // This callback is called by osmium::apply for each area in the data.
    void area(const osmium::Area& area) {
        try {
            std::cout << std::endl;
            for (auto &tag : area.tags())
                std::cout << std::string(tag.key()) + "=" + tag.value() << std::endl;
            for (const auto& item : area) {
                if (item.type() == osmium::item_type::outer_ring) {
                    auto& ring = static_cast<const osmium::OuterRing&>(item);
                    for (auto const &p : ring) {
                        _boundingBox._p[0] = std::min(p.lon(), _boundingBox._p[0]);
                        _boundingBox._p[1] = std::min(p.lat(), _boundingBox._p[1]);
                        _boundingBox._d[0] = std::max(p.lon(), _boundingBox._d[0]);
                        _boundingBox._d[1] = std::max(p.lat(), _boundingBox._d[1]);
                    }
                }
            }
        } catch (const osmium::geometry_error& e) {
            std::cout << "GEOMETRY ERROR: " << e.what() << "\n";
        }
    }
    Segmentd &_boundingBox;
};

void print_help() {
    std::cout << "osmium_area_test [OPTIONS] OSMFILE\n\n"
              << "Read OSMFILE and build multipolygons from it.\n"
              << "\nOptions:\n"
              << "  -h, --help           This help message\n"
              << "  -w, --dump-wkt       Dump area geometries as WKT\n"
              << "  -o, --dump-objects   Dump area objects\n";
}

void print_usage(const char* prgname) {
    std::cerr << "Usage: " << prgname << " [OPTIONS] OSMFILE\n";
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    float val = (percentage * 100);
    percentage = std::min(percentage, 1.0);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3f2%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

long int fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0) {
        std::cout << filename << " file size = " << st.st_size << std::endl;
        return st.st_size;
    }
    return -1; 
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    char const *fileName = 0;
    if (argc == 2) {
        fileName = argv[1];
    } else if (argc == 4) {
        fileName = argv[2];
        std::string label = argv[1];
        std::string key = label.substr(0, label.rfind('='));
        label = label.substr(label.rfind('=')+1);

        std::cout << DEBUGVAR(key) << std::endl;
        std::cout << DEBUGVAR(label) << std::endl;

        Segmentd boundingBox = {{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}};


        osmium::handler::DynamicHandler handler;
        handler.set<WKTAreaFinder>(boundingBox);
        osmium::io::File input_file{fileName};
        osmium::area::Assembler::config_type assembler_config;

        osmium::TagsFilter filter{false};
        filter.add_rule(true, key, label);
        // filter.add_rule(true, "landuse", "forest");
        // filter.add_rule(true, "natural", "wood");
        osmium::area::MultipolygonManager<osmium::area::Assembler> mp_manager{assembler_config, filter};
        // osmium::area::MultipolygonManager<osmium::area::Assembler> mp_manager{assembler_config};

        std::cerr << "Pass 1...\n";
        osmium::relations::read_relations(input_file, mp_manager);
        std::cerr << "Pass 1 done\n";

        index_type index;
        location_handler_type location_handler{index};

        location_handler.ignore_errors();

        std::cerr << "Pass 2...\n";
        osmium::io::Reader reader{input_file};
        osmium::apply(reader, location_handler, mp_manager.handler([&handler](osmium::memory::Buffer&& buffer) {
            osmium::apply(buffer, handler);
        }));
        reader.close();
        std::cerr << "Pass 2 done\n";

        fileName = strdup((label + std::string(".osm.pbf")).c_str());
        std::cout << fileName << std::endl;
        std::cout << DEBUGVAR(boundingBox) << std::endl;
        std::string cmd = std::string("osmium extract --bbox ")+std::to_string(boundingBox._p[0])+","+std::to_string(boundingBox._p[1])+","+std::to_string(boundingBox._d[0])+","+std::to_string(boundingBox._d[1])+" "+argv[3]+" -o "+fileName+ " --overwrite";
        std::cout << DEBUGVAR(cmd) << std::endl;
        system(cmd.c_str());

    } else if (argc == 6) {
        fileName = "tmp.osm.pbf";

        std::string cmd = std::string("osmium extract --bbox ")+argv[1]+","+argv[2] +","+ argv[3] +","+ argv[4]+" "+argv[5]+" -o "+fileName+ " --overwrite";
        std::cout << DEBUGVAR(cmd) << std::endl;
        system(cmd.c_str());
    }

    try {
        std::vector<StreeShape> segs;
        Segmentd boundingBox = {{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}};

        osmium::handler::DynamicHandler handler;
        handler.set<WKTDump>(segs, boundingBox);
        osmium::io::File input_file{fileName};
        osmium::area::Assembler::config_type assembler_config;

        // osmium::TagsFilter filter{false};
        // filter.add_rule(true, "name", "Pont-à-Marcq");
        // filter.add_rule(true, "landuse", "forest");
        // filter.add_rule(true, "natural", "wood");
        // osmium::area::MultipolygonManager<osmium::area::Assembler> mp_manager{assembler_config, filter};
        osmium::area::MultipolygonManager<osmium::area::Assembler> mp_manager{assembler_config};

        std::cerr << "Pass 1...\n";
        osmium::relations::read_relations(input_file, mp_manager);
        std::cerr << "Pass 1 done\n";

        index_type index;
        location_handler_type location_handler{index};

        location_handler.ignore_errors();

        std::cerr << "Pass 2...\n";
        osmium::io::Reader reader{input_file};
        osmium::apply(reader, location_handler, mp_manager.handler([&handler](osmium::memory::Buffer&& buffer) {
            osmium::apply(buffer, handler);
        }));

        std::cout << "\r" << segs.size() << " areas found" << std::endl;
        reader.close();
        std::cerr << "Pass 2 done\n";





        std::cout << DEBUGVAR(boundingBox) << std::endl;
        boundingBox._d -= boundingBox._p;

        Vector2d center = boundingBox._p+boundingBox._d/2;
        std::cout << DEBUGVAR(boundingBox) << std::endl;

        // do projection && update bounding box
        boundingBox = {{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}};
        for (auto &shape : segs) {
            for (auto &p : shape._points) {
                p[0] *= cos(p[1]/180*M_PI);

                boundingBox._p[0] = std::min(p[0], boundingBox._p[0]);
                boundingBox._p[1] = std::min(p[1], boundingBox._p[1]);
                boundingBox._d[0] = std::max(p[0], boundingBox._d[0]);
                boundingBox._d[1] = std::max(p[1], boundingBox._d[1]);
            }
        }
        boundingBox._d -= boundingBox._p;
        center = boundingBox._p+boundingBox._d/2;

        // center points && scale to meters && update bounding box
        for (auto &shape : segs) {
            for (auto &p : shape._points) {
                p -= center;
                p *= 111194.0;

                boundingBox._p[0] = std::min(p[0], boundingBox._p[0]);
                boundingBox._p[1] = std::min(p[1], boundingBox._p[1]);
                boundingBox._d[0] = std::max(p[0], boundingBox._d[0]);
                boundingBox._d[1] = std::max(p[1], boundingBox._d[1]);
            }
        }
        boundingBox._d -= boundingBox._p;
        center = boundingBox._p+boundingBox._d/2;

        // calculate shapes bounding box
        // for (auto &shape : segs) {
        //     shape._boundingBox = {DBL_MAX, DBL_MAX, DBL_MIN, DBL_MIN};
        //     for (auto &p : shape._points) {
        //         shape._boundingBox._p[0] = std::min(p[0], shape._boundingBox._p[0]);
        //         shape._boundingBox._p[1] = std::min(p[1], shape._boundingBox._p[1]);
        //         shape._boundingBox._d[0] = std::max(p[0], shape._boundingBox._d[0]);
        //         shape._boundingBox._d[1] = std::max(p[1], shape._boundingBox._d[1]);
        //     }
        // }        

        class UniTreeObj : public Vector2d {
            public:
            UniTreeObj(Vector2d const &vec, StreeShape &shape) : 
                Vector2d(vec),
                _shape(shape)
            {
                // std::cout << DEBUGVAR(data[0]) << std::endl;
                // std::cout << DEBUGVAR(data[1]) << std::endl;
            }

            StreeShape &_shape;
        };


        

        // Zone<double, 2> zone(boundingBox._p+boundingBox._d/2, boundingBox._d/2);

        // fill UniTreeZone
        auto treeZone = std::make_unique<UniTree<UniTreeObj, Vector2d, 2>>(boundingBox._p+boundingBox._d/2, boundingBox._d/2);

        std::vector<UniTreeObj> UniTreeObjects;
        uint allocSize = 0;
        for (StreeShape &shape : segs)
            allocSize += shape._points.size();
        std::cout << DEBUGVAR(allocSize) << std::endl;
        UniTreeObjects.reserve(allocSize);
        uint currentSize = 0;
        srand(14);
        std::cout << "building map..." << std::endl;
        for (StreeShape &shape : segs) {
            for (Vector2d const &p : shape._points) {
                UniTreeObjects.emplace_back(p+Vector2d{(currentSize%2+1)*1./10000*currentSize, (currentSize%2+0)*1./10000*currentSize}, shape); // TODO this is bad, I need to make a new spatialisation tool
                treeZone->addData(&UniTreeObjects.back());
                currentSize += 1;
            }
            printProgress(currentSize/(float)allocSize);
        }
        std::cout << "done" << std::endl;



        Vector2f _cameraOffset = {0, 0};
        float _cameraAngle = 0;
        Mat3 _matWorld;
        Mat3 _matWorldInv;
        float _camScale = 1;//1./winSize[0]/std::max(boundingBox._d[0], boundingBox._d[1]);

        Vector2i winSize = {800, 800};

        sf::RenderWindow win(sf::VideoMode(winSize[0], winSize[1]), "Map Gen");
        win.setFramerateLimit(60);

        sf::Font font;
        font.loadFromFile("assets/fonts/ARIBL0.ttf");
        sf::Text text;
        text.setFont(font);
        uint textSize = 12;
        text.setCharacterSize(textSize);
        text.setFillColor(sf::Color::White);

        bool displayLabel = false;

        std::shared_ptr<std::vector<UniTreeObj *>> uniTreeBuffer = std::make_shared<std::vector<UniTreeObj*>>();
        while (win.isOpen()) {
            win.clear();


            float camcos = -cos(-_cameraAngle) * 1.0/_camScale*8;
            float camsin = -sin(-_cameraAngle) * 1.0/_camScale*8;

            sf::Keyboard::isKeyPressed(sf::Keyboard::W) ? _cameraOffset[0] += camsin, _cameraOffset[1] -= camcos : 0;
            sf::Keyboard::isKeyPressed(sf::Keyboard::S) ? _cameraOffset[0] -= camsin, _cameraOffset[1] += camcos : 0;
            sf::Keyboard::isKeyPressed(sf::Keyboard::A) ? _cameraOffset[0] -= camcos, _cameraOffset[1] -= camsin : 0;
            sf::Keyboard::isKeyPressed(sf::Keyboard::D) ? _cameraOffset[0] += camcos, _cameraOffset[1] += camsin : 0;

            sf::Keyboard::isKeyPressed(sf::Keyboard::Z) ? _cameraAngle += 0.025 : 0;
            sf::Keyboard::isKeyPressed(sf::Keyboard::X) ? _cameraAngle -= 0.025 : 0;
            
            sf::Keyboard::isKeyPressed(sf::Keyboard::E) ? _camScale *= 1.02 : 0;
            sf::Keyboard::isKeyPressed(sf::Keyboard::Q) ? _camScale *= 0.98 : 0;

            Vector2f mousePos;
            {
                sf::Vector2i vec = sf::Mouse::getPosition(win);
                mousePos = {(float)vec.x, (float)vec.y};
            }
            Vector2f mousePosCenterRelative = mousePos-winSize.cast<float>()/2;


            sf::Event event;
            while (win.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    win.close();
                } else if (event.type == sf::Event::KeyPressed) {
                    if (event.key.code == sf::Keyboard::L) {
                        displayLabel = !displayLabel;
                    } else if (event.key.code == sf::Keyboard::M) {
                        std::cout << "save..." << std::endl;
                        
                        struct QuickMatch {
                            char const *key;
                            char const *value;
                        };

                        std::vector<QuickMatch> toKeep = {
                            {"building", 0},
                        };

                        std::vector<StreeShape> toSave;
                        toSave.reserve(segs.size());
                        for (StreeShape &shape : segs) {
                            bool keep = false;
                            for (std::string const &s : shape._labels) {
                                std::string key = s.substr(0, s.rfind('='));
                                std::string value = s.substr(s.rfind('=')+1);
                                for (QuickMatch const &ss : toKeep) {
                                    if (key == ss.key && (ss.value == 0 || value == ss.value)) {
                                        keep = true;
                                        goto keepIt;
                                    }
                                }
                            }
                            if (0) {
                                keepIt:
                                toSave.emplace_back(shape);
                            }
                        }

                        ByteObject obj;
                        obj << toSave;

                        std::ofstream file("../BoatFight/map.bmap", std::fstream::binary);
                        file << obj;
                        std::cout << "save...ok" << std::endl;

                    } else if (event.key.code == sf::Keyboard::N) {
                        // segs.clear();
                        UniTreeObjects.clear();

                        std::ifstream file("../BoatFight/map.bmap", std::ifstream::binary);

                        ByteObject obj;
                        file >> obj;

                        obj >> segs;


                        // update bounding box
                        for (auto &shape : segs) {
                            for (auto &p : shape._points) {
                                boundingBox._p[0] = std::min(p[0], boundingBox._p[0]);
                                boundingBox._p[1] = std::min(p[1], boundingBox._p[1]);
                                boundingBox._d[0] = std::max(p[0], boundingBox._d[0]);
                                boundingBox._d[1] = std::max(p[1], boundingBox._d[1]);
                            }
                        }
                        boundingBox._d -= boundingBox._p;
                        center = boundingBox._p+boundingBox._d/2;

                        treeZone = std::make_unique<UniTree<UniTreeObj, Vector2d, 2>>(boundingBox._p+boundingBox._d/2, boundingBox._d/2);


                        uint allocSize = 0;
                        for (StreeShape &shape : segs)
                            allocSize += shape._points.size();
                        std::cout << DEBUGVAR(allocSize) << std::endl;
                        UniTreeObjects.reserve(allocSize);
                        uint currentSize = 0;
                        std::cout << "building map..." << std::endl;
                        srand(14);
                        for (StreeShape &shape : segs) {
                            for (Vector2d const &p : shape._points) {
                                UniTreeObjects.emplace_back(p+Vector2d{(currentSize%2+1)*1./10000*currentSize, (currentSize%2+0)*1./10000*currentSize}, shape); // TODO this is bad, I need to make a new spatialisation tool

                                treeZone->addData(&UniTreeObjects.back());
                                currentSize += 1;
                            }
                            printProgress(currentSize/(float)allocSize);
                        }
                    }
                }
            }



            Mat3 matrix;
            matrix.tx(_cameraOffset[0]);
            matrix.ty(_cameraOffset[1]);
            matrix.rrz(_cameraAngle);
            matrix.scale(_camScale);
            matrix.ttx(-mousePosCenterRelative[0]*0.5);
            matrix.tty(-mousePosCenterRelative[1]*0.5);
            _matWorld = matrix;
            _matWorldInv = matrix.inv();

            sf::CircleShape ship;
            ship.setFillColor(sf::Color::Red);
            ship.setRadius(_camScale*1);
            ship.setPosition(-mousePosCenterRelative[0]*0.8+winSize[0]/2, -mousePosCenterRelative[1]*0.8+winSize[1]/2);

            win.draw(ship);

            Vector2d p1 = (_matWorldInv * (winSize.cast<float>()*Vector2f{-1, -1})/2).cast<double>();  
            Vector2d p2 = (_matWorldInv * (winSize.cast<float>()*Vector2f{-1,  1})/2).cast<double>();
            Vector2d p3 = (_matWorldInv * (winSize.cast<float>()*Vector2f{ 1, -1})/2).cast<double>();
            Vector2d p4 = (_matWorldInv * (winSize.cast<float>()*Vector2f{ 1,  1})/2).cast<double>();
            Vector2d center = (p1+p4);

            Segmentd cameraBoundingBox = {
                std::min(std::min(p1[0], p2[0]), std::min(p3[0], p4[0])),
                std::min(std::min(p1[1], p2[1]), std::min(p3[1], p4[1])),
                std::max(std::max(p1[0], p2[0]), std::max(p3[0], p4[0])),
                std::max(std::max(p1[1], p2[1]), std::max(p3[1], p4[1])),
            };
            cameraBoundingBox._d -= cameraBoundingBox._p;

            Vector2d size = cameraBoundingBox._d/2;

            uniTreeBuffer->clear();
            treeZone->getInArea(center, size, uniTreeBuffer);
            


            // uint allVertexCursor = 0;
            // int shapeDraw = 1024*256;
            uint maxTextDraw = 16;

            Vector2f winSizef = winSize.cast<float>()/2; 
            for (auto &obj : *uniTreeBuffer) {
                auto &shape = obj->_shape;
                if (shape._drawed)  
                        continue;
                shape._drawed = true;

                // display individual
                // shapeDraw -= shape._points.size();
                // if (shapeDraw <= 0) {
                //     break;
                // }
                sf::VertexArray lines(sf::LineStrip, shape._points.size());
                uint i = 0;  
                for (auto &p : shape._points) {
                    Vector2f p1 = _matWorld * p.cast<float>() + winSizef;
                    lines[i] = sf::Vertex(sf::Vector2f(p1[0], p1[1]), shape._color);
                    ++i;
                }

                win.draw(lines);

           
                if (displayLabel) {
                    Vector2f p1 = _matWorld * shape._points[0].cast<float>() + winSize.cast<float>()/2;
                    if (shape._labels.size() && maxTextDraw > 0 && 0 <= p1[0] && p1[0] < winSize[0] && 0 <= p1[1] && p1[1] < winSize[1]) {
                        --maxTextDraw;
                        text.setPosition(p1[0], p1[1]);
                        for (std::string const &string : shape._labels) {
                            text.setString(string);
                            win.draw(text);
                            text.move({0, textSize*1.1});
                        }
                    }
                }
            }
            for (auto &obj : *uniTreeBuffer)
                obj->_shape._drawed = false;
            win.display();
        }





    } catch (const std::exception& e) {
        // All exceptions used by the Osmium library derive from std::exception.
        std::cerr << e.what() << '\n';
        return 1;
    }
}

