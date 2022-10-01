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

#include <float.h>

#include "utils.hpp"

#include "mat/mat.hpp"

// The type of index used. This must match the include file above
using index_type = osmium::index::map::FlexMem<osmium::unsigned_object_id_type, osmium::Location>;

// The location handler always depends on the index type
using location_handler_type = osmium::handler::NodeLocationsForWays<index_type>;

class StreeShape {
public:
    // Segmentd _seg;
    Color _color;
    std::vector<Vector2d> _points;
    std::vector<std::string> _labels;
};


class Match {
public:
    char const *_key;
    char const *_tag;
    Color const _color;
};

std::vector<Match> g_colorVector = {
    {"wall", "no", {0.25, 0.25, 0.25, 1.}},
    {"amenity", "parking_space", {0.25, 0.25, 0.25, 1.}},

    {"amenity", "parking", {0.5, 0.25, 0.25, 1.}},
    {"amenity", "parking", {0.5, 0.25, 0.25, 1.}},
    {"highway", "service", {0.5, 0.25, 0.25, 1.}},
    
    {"amenity", "school", {0.5, 0.5, 0.5, 1}},

    {"leisure", "track", {0.2, 0.2, 0.2, 1.}},

    {"landuse", "residential", {1., 0., 0., 1.}},
    
    {"water", "basin", {0., 0., 1., 1.}},
    {"natural", "water", {0.75, 0.75, 0.75, 1.f}},

    {"leisure", "garden", {0., 1., 0., 1.}},
    {"leisure", "golf_course", {0., 1., 0., 1.}},
    {"landuse", "grass", {0., 1., 0., 1.}},
    
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

    WKTDump(std::deque<StreeShape> &segs, Segmentd &boundingBox) :
        osmium::handler::Handler(),
        _segs(segs),
        _boundingBox(boundingBox)
    {

    }

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

                    seg._color = color;
                    if (!find) {
                        std::cout << std::endl;
                        for (auto &tag : area.tags()) {
                            std::cout << tag << std::endl;
                        }
                    }

                    auto& ring = static_cast<const osmium::OuterRing&>(item);

                    for (auto const &p : ring) {
                        // Vector2d pos = {p.lon()*cos(p.lat()/180*M_PI), -p.lat()};
                        Vector2d pos = {p.lon(), -p.lat()};

                        if (pos[0] < _boundingBox._p[0]) _boundingBox._p[0] = pos[0];
                        if (pos[1] < _boundingBox._p[1]) _boundingBox._p[1] = pos[1];
                        if (pos[0] > _boundingBox._d[0]) _boundingBox._d[0] = pos[0];
                        if (pos[1] > _boundingBox._d[1]) _boundingBox._d[1] = pos[1];

                        seg._points.push_back(pos);
                    }
                    _segs.push_back(seg);
                }
            }


        } catch (const osmium::geometry_error& e) {
            std::cout << "GEOMETRY ERROR: " << e.what() << "\n";
        }
    }

    std::deque<StreeShape> &_segs;
    Segmentd &_boundingBox;

}; // class WKTDump

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

int main(int argc, char* argv[]) {
    if (argc != 2) {
        print_usage(argv[0]);
        return 1;
    }

    // std::make_unique<strategy_smart::Strategy>(m_extracts, m_options);

    try {
        std::deque<StreeShape> segs;
        Segmentd boundingBox = {{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}};

        // Initialize an empty DynamicHandler. Later it will be associated
        // with one of the handlers. You can think of the DynamicHandler as
        // a kind of "variant handler" or a "pointer handler" pointing to the
        // real handler.
        osmium::handler::DynamicHandler handler;

        handler.set<WKTDump>(segs, boundingBox);

        osmium::io::File input_file{argv[1]};

        // Configuration for the multipolygon assembler. Here the default settings
        // are used, but you could change multiple settings.
        osmium::area::Assembler::config_type assembler_config;

        // osmium::TagsFilter filter{false};
        // filter.add_rule(true, "landuse", "forest");
        // filter.add_rule(true, "natural", "wood");
        // osmium::area::MultipolygonManager<osmium::area::Assembler> mp_manager{assembler_config, filter};
        osmium::area::MultipolygonManager<osmium::area::Assembler> mp_manager{assembler_config};

        // We read the input file twice. In the first pass, only relations are
        // read and fed into the multipolygon manager.
        std::cerr << "Pass 1...\n";
        osmium::relations::read_relations(input_file, mp_manager);
        std::cerr << "Pass 1 done\n";

        // Output the amount of main memory used so far. All multipolygon relations
        // are in memory now.
        std::cerr << "Memory:\n";
        osmium::relations::print_used_memory(std::cerr, mp_manager.used_memory());

        // The index storing all node locations.
        index_type index;

        // The handler that stores all node locations in the index and adds them
        // to the ways.
        location_handler_type location_handler{index};

        // If a location is not available in the index, we ignore it. It might
        // not be needed (if it is not part of a multipolygon relation), so why
        // create an error?
        location_handler.ignore_errors();

        // On the second pass we read all objects and run them first through the
        // node location handler and then the multipolygon collector. The collector
        // will put the areas it has created into the "buffer" which are then
        // fed through our "handler".
        std::cerr << "Pass 2...\n";
        osmium::io::Reader reader{input_file};
        osmium::apply(reader, location_handler, mp_manager.handler([&handler](osmium::memory::Buffer&& buffer) {
            osmium::apply(buffer, handler);
        }));
        reader.close();
        std::cerr << "Pass 2 done\n";


        // Output the amount of main memory used so far. All complete multipolygon
        // relations have been cleaned up.
        std::cerr << "Memory:\n";
        osmium::relations::print_used_memory(std::cerr, mp_manager.used_memory());

        // If there were multipolgyon relations in the input, but some of their
        // members are not in the input file (which often happens for extracts)
        // this will write the IDs of the incomplete relations to stderr.
        std::vector<osmium::object_id_type> incomplete_relations_ids;
        mp_manager.for_each_incomplete_relation([&](const osmium::relations::RelationHandle& handle){
            incomplete_relations_ids.push_back(handle->id());
        });
        if (!incomplete_relations_ids.empty()) {
            std::cerr << "Warning! Some member ways missing for these multipolygon relations:";
            for (const auto id : incomplete_relations_ids) {
                std::cerr << " " << id;
            }
            std::cerr << "\n";
        }


        std::cout << DEBUGVAR(boundingBox) << std::endl;
        boundingBox._d -= boundingBox._p;

        Vector2d center = boundingBox._p+boundingBox._d/2;

        // do projection and update bounding box
        boundingBox = {{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}};
        for (auto &shape : segs) {
            for (auto &p : shape._points) {
                p[0] -= center[0];
                p[0] *= cos(p[1]/180*M_PI);

                if (p[0] < boundingBox._p[0]) boundingBox._p[0] = p[0];
                if (p[1] < boundingBox._p[1]) boundingBox._p[1] = p[1];
                if (p[0] > boundingBox._d[0]) boundingBox._d[0] = p[0];
                if (p[1] > boundingBox._d[1]) boundingBox._d[1] = p[1];
            }
        }
        boundingBox._d -= boundingBox._p;




        // center and scale
        double scale = 800/std::max(boundingBox._d[0], boundingBox._d[1]);
        std::cout << DEBUGVAR(boundingBox) << std::endl;
        for (auto &shape : segs) {
            for (auto &p : shape._points) {
                p = (p-boundingBox._p)*scale;
            }
        }

        std::cout << DEBUGVAR(segs.size()) << std::endl;
        Vector2f _cameraOffset = {0, 0};
        float _cameraAngle = 0;
        Mat3 _matWorld;
        Mat3 _matWorldInv;
        float _camScale = 1;

        sf::RenderWindow win(sf::VideoMode(800, 800), "Map Gen");
        win.setFramerateLimit(60);

        sf::Font font;
        font.loadFromFile("assets/fonts/ARIBL0.ttf");
        sf::Text text;
        text.setFont(font);
        uint textSize = 12;
        text.setCharacterSize(textSize);
        text.setFillColor(sf::Color::White);

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

            sf::Event event;
            while (win.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    win.close();
                }
            }


            Mat3 matrix;
            matrix.tx(_cameraOffset[0]);
            matrix.ty(_cameraOffset[1]);
            matrix.rrz(_cameraAngle);
            matrix.scale(_camScale);
            _matWorld = matrix;
            _matWorldInv = matrix.inv();

            

            uint maxTextDraw = 16;
            for (auto &shape : segs) {
                sf::VertexArray lines(sf::LineStrip, shape._points.size());
                uint i = 0;
                for (auto &p : shape._points) {
                    Vector2f p1 = _matWorld * p.cast<float>() + 800/2;
                    lines[i] = sf::Vertex(sf::Vector2f(p1[0], p1[1]), {(uint8_t)(shape._color.r*255), (uint8_t)(shape._color.g*255), (uint8_t)(shape._color.b*255), (uint8_t)(shape._color.a*255)}),
                    ++i;
                }

                win.draw(lines);

                Vector2f p1 = _matWorld * shape._points[0].cast<float>() + 800/2;
                if (shape._labels.size() && maxTextDraw > 0 && 0 <= p1[0] && p1[0] < 800 && 0 <= p1[1] && p1[1] < 800) {
                    --maxTextDraw;
                    text.setPosition(p1[0], p1[1]);
                    for (std::string const &string : shape._labels) {
                        text.setString(string);
                        win.draw(text);
                        text.move({0, textSize*1.1});
                    }
                }
            }
            win.display();
        }





    } catch (const std::exception& e) {
        // All exceptions used by the Osmium library derive from std::exception.
        std::cerr << e.what() << '\n';
        return 1;
    }
}

