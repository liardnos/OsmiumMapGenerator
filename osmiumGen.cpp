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

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>

#include <curlpp/cURLpp.hpp>
#include <curlpp/Easy.hpp>
#include <curlpp/Options.hpp>

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
#include <sstream>

#include <fstream>

#include "ByteObject.hpp"

#include "utils.hpp"

#include "mat/mat.hpp"

#include "json.hpp"

// #include "spatial/UniTreeZone.hpp"
#include "spatial/UniTreeZone.hpp"

// The type of index used. This must match the include file above
using index_type = osmium::index::map::FlexMem<osmium::unsigned_object_id_type, osmium::Location>;

// The location handler always depends on the index type
using location_handler_type = osmium::handler::NodeLocationsForWays<index_type>;

class StreeShape {
public:
    // Segmentd _seg;
    sf::Color _color;
    std::vector<Vector3d> _points;
    std::vector<std::string> _labels;
    float _height = 0;
    Segmentf _boundingBox;

    friend ByteObject &operator<<(ByteObject &obj, StreeShape const &shape) {
        std::vector<Vector3f> vec;
        vec.reserve(shape._points.size());
        for (Vector3d const &p : shape._points)
            vec.push_back(p.cast<float>());

        obj << shape._color << vec << shape._labels;
        return obj;
    }

    friend ByteObject &operator>>(ByteObject &obj, StreeShape &shape) {
        std::vector<Vector3f> vec;
        obj >> shape._color >> vec >> shape._labels;
        shape._points.clear();
        shape._points.reserve(vec.size());
        for (Vector3f const &p : vec)
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
            // search for height
            try {
                if (area.tags().has_key("height")) {
                    seg._height = std::stof(area.tags().get_value_by_key("height"));
                } else if (area.tags().has_key("min_height")) {
                    seg._height = std::stof(area.tags().get_value_by_key("min_height"));
                } else if (area.tags().has_key("level")) {
                    seg._height = std::stof(area.tags().get_value_by_key("level"))*5;
                }
            } catch (std::exception &e) {
                std::cout << "stof failed" << e.what() << std::endl;
                // std::cout << area.tags().get_value_by_key("height") << std::endl;
                // std::cout << area.tags().get_value_by_key("min_height") << std::endl;
                // std::cout << area.tags().get_value_by_key("level") << std::endl;
            }
            
            // height=2.3
            // min_height=2.3
            // level=-2
            for (const auto& item : area) {
                if (item.type() == osmium::item_type::outer_ring) {
                    // bool find = false;
                    Color color = {1, 1, 1, 1};
                    for (auto const &match : g_colorVector) {
                        if (match._tag == 0) {
                            if (area.tags().has_key(match._key)) {
                                color = match._color;
                                // find = true;
                                break;
                            }
                        } else if (area.tags().has_tag(match._key, match._tag)) {
                            color = match._color;
                            // find = true;
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
                        Vector3d pos = {p.lon(), p.lat(), -seg._height};


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
        int r = system(cmd.c_str());
        (void)r;

    } else if (argc == 6) {
        fileName = "tmp.osm.pbf";

        std::string cmd = std::string("osmium extract --bbox ")+argv[1]+","+argv[2] +","+ argv[3] +","+ argv[4]+" "+argv[5]+" -o "+fileName+ " --overwrite";
        std::cout << DEBUGVAR(cmd) << std::endl;
        int r = system(cmd.c_str());
        (void)r;

    }

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
    std::cerr << "Pass 2 done" << std::endl;

    boundingBox._d -= boundingBox._p;

    Vector2d center = boundingBox._p+boundingBox._d/2;
    
    Vector2f _cameraOffset = {0, 0};
    float _cameraAngle = 0;
    Mat3 _matWorld;
    Mat3 _matWorldInv;
    float _camScale = 10000;//1./winSize[0]/std::max(boundingBox._d[0], boundingBox._d[1]);

    std::cout << DEBUGVAR(boundingBox) << std::endl;
    std::mutex segsMut;



    struct ThreadShared {
        ThreadShared(std::vector<StreeShape> &segs) {
            segIt = segs.begin();
            pIt = segIt->_points.begin();
        }
        uint const threadsCount = 2;
        int const requestCount = 128;
        int aliveThreads = threadsCount;

        std::mutex iteratorMut;
        std::vector<StreeShape>::iterator segIt;
        std::vector<Vector3d>::iterator pIt;
        int currentOffset = 0;
        float heightmoy = 0;
        bool finish = false;
    };


    setenv("PYTHONPATH",".",1);
    Py_Initialize();
    PyObject *pFileName, *pModule, *pDict, *presult, *pValue;

    pFileName = PyUnicode_FromString((char *)"elevation");
    pModule = PyImport_Import(pFileName);
    pDict = PyModule_GetDict(pModule);

    PyObject *pFunc = PyDict_GetItemString(pDict, (char*)"lookUpElevation");


    std::unique_ptr<ThreadShared> threadShared = std::make_unique<ThreadShared>(segs);


    if (1) { // get the height map
        int pCount = 0; // number of point in the graph
        for (auto &seg : segs) {
            pCount += seg._points.size();
            seg._color.a = 50; 
        }
        std::cout << DEBUGVAR(pCount) << std::endl;
        std::shared_ptr<std::vector<Vector3f>> requestBuffer = std::make_shared<std::vector<Vector3f>>(pCount);
        { // backup original lat lon for elevaltion requests
            
            int i = 0;
            for (auto seg : segs)
                for (auto p : seg._points)
                    (*requestBuffer)[i++] = p.cast<float>();
        }

        std::cout << "fetching elevation" << std::endl;

        std::thread *threads[threadShared->threadsCount];
        for (uint threadId = 0; threadId < threadShared->threadsCount; threadId++) {
            threads[threadId] = new std::thread(
            [&threads, &threadShared, &segs, pCount, requestBuffer, &segsMut, &pFileName, &pModule, &pDict, &presult, &pValue, pFunc](){
                int const &requestCount = threadShared->requestCount;
                std::mutex &iteratorMut = threadShared->iteratorMut;
                bool &finish = threadShared->finish;
                auto &segIt = threadShared->segIt;
                auto &pIt = threadShared->pIt;
                auto &currentOffset = threadShared->currentOffset;
                auto &heightmoy = threadShared->heightmoy;
                auto &aliveThreads = threadShared->aliveThreads;
                auto const &threadsCount = threadShared->threadsCount;


                Vector3f buf[requestCount];
                auto lsegIt = segIt;
                auto lpIt = pIt;

                while (1) {
                    uint count = 0;
                    {
                        std::lock_guard<std::mutex> const lockable(iteratorMut);
                        lsegIt = segIt;
                        lpIt = pIt;
                        if (finish)
                            break;
                        // advance iterator for next thread
                        for (; count < requestCount;) {
                            count++;
                            pIt++;
                            if (pIt == segIt->_points.end()) {
                                segIt++;
                                if (segIt == segs.end()) {
                                    finish = true;
                                    break;
                                }
                                pIt = segIt->_points.begin();
                            }
                        }
                        // copy data in buffer                
                        for (uint i = 0; i < count; i++)
                            buf[i] = (*requestBuffer)[currentOffset++];
                    }

                    
                    // curlpp::Cleanup myCleanup;
                    // curlpp::Easy myRequest;
                    // std::string requestString = "127.0.0.1:8080/api/v1/lookup?locations=";
                    // for (uint i = 0; i < count; i++) {
                    //     if (i) requestString += "|";
                    //     requestString += std::to_string(buf[i][1]) + "," + std::to_string(buf[i][0]);
                    // }
                    // myRequest.setOpt<curlpp::options::Url>(requestString);
                    // std::stringstream ss;
                    // ss << myRequest;
                    // nlohmann::json data = nlohmann::json::parse(ss.str());

                    { // data in segs
                        // std::lock_guard<std::mutex> const lockable(segsMut);
                        for (uint i = 0; i < count; i++) {
                            std::cout << "build value" << std::endl;
                            pValue = Py_BuildValue("f f", buf[i][1], buf[i][0]);
                            presult = PyObject_CallObject(pFunc, pValue);
                            (*lpIt)[2] -= (double)PyFloat_AS_DOUBLE(presult);
                            heightmoy += (*lpIt)[2];
                            lpIt++;
                            if (lpIt == lsegIt->_points.end()) {
                                lsegIt->_color.a = 255; 
                                lsegIt++;
                                lpIt = lsegIt->_points.begin();
                            }


                        }
                        printProgress((float)currentOffset/pCount);
                    }
                }


                {
                    iteratorMut.lock();
                    aliveThreads--;
                    if (aliveThreads == 0) {
                        printProgress(1);
                        heightmoy /= pCount;
                        std::lock_guard<std::mutex> segsMutLock(segsMut);
                        for (auto &seg : segs)
                            for (auto &p : seg._points)
                                p[2] -= heightmoy;
                        iteratorMut.unlock();
                        goto next;
                    }
                    iteratorMut.unlock();
                    next:;
                }
            });
        }
    }


    // do projection && update bounding box
    boundingBox = {{DBL_MAX, DBL_MAX}, {-DBL_MAX, -DBL_MAX}};
    for (auto &shape : segs) {
        for (auto &p : shape._points) { 
            p[0] *= cos(p[1]/180*M_PI);
            p[1] *= -1;

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
        shape._boundingBox = {FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (auto &p : shape._points) {
            p[0] -= center[0];
            p[1] -= center[1];
            p[0] *= 111194.0;
            p[1] *= 111194.0;

            boundingBox._p[0] = std::min(p[0], boundingBox._p[0]);
            boundingBox._p[1] = std::min(p[1], boundingBox._p[1]);
            boundingBox._d[0] = std::max(p[0], boundingBox._d[0]);
            boundingBox._d[1] = std::max(p[1], boundingBox._d[1]);

            shape._boundingBox._p[0] = std::min((float)p[0], shape._boundingBox._p[0]);
            shape._boundingBox._p[1] = std::min((float)p[1], shape._boundingBox._p[1]);
            shape._boundingBox._d[0] = std::max((float)p[0], shape._boundingBox._d[0]);
            shape._boundingBox._d[1] = std::max((float)p[1], shape._boundingBox._d[1]);
        }
        shape._boundingBox._d -= shape._boundingBox._p;
    }
    boundingBox._d -= boundingBox._p;
    center = boundingBox._p+boundingBox._d/2;


    

    // Zone<double, 2> zone(boundingBox._p+boundingBox._d/2, boundingBox._d/2);

    // fill UniTreeZone
    std::unique_ptr<UniTreeZone<float, StreeShape, 2>> treeZone = std::make_unique<UniTreeZone<float, StreeShape, 2>>(Zone<float, 2>(boundingBox.cast<float>()));

    uint allocSize = segs.size();
    std::cout << DEBUGVAR(allocSize) << std::endl;
    uint currentSize = 0;
    std::cout << "building map..." << std::endl;
    for (StreeShape &shape : segs) {
        treeZone->addData(Zone<float, 2>(shape._boundingBox), &shape);
        if (!(currentSize%1000))
            printProgress(currentSize/(float)allocSize);
        currentSize += 1;
    }
    printProgress(1);
    std::cout << "done" << std::endl;


    Vector2i winSize = {1900, 1000};
    Vector2f winSizef = winSize.cast<float>(); 


    float targetFrameRateMonving = 30;
    float targetFrameRate = 30;
    sf::RenderWindow win(sf::VideoMode(winSize[0], winSize[1]), "Map Gen");
    win.setFramerateLimit(targetFrameRateMonving*1.1);

    sf::Font font;
    font.loadFromFile("assets/fonts/ARIBL0.ttf");
    sf::Text text;
    text.setFont(font);
    uint textSize = 12;
    text.setCharacterSize(textSize);
    text.setFillColor(sf::Color::White);

    bool displayLabel = false;
    float camPitch = 0;

    std::shared_ptr<std::vector<std::shared_ptr<UniTreeZone<float, StreeShape, 2>::Storage>>> uniTreeBuffer = std::make_shared<std::vector<std::shared_ptr<UniTreeZone<float, StreeShape, 2>::Storage>>>();
    // Vector2f mousePosPrev = {0, 0};

    float minSizeGet = 0;
    float minSizeGetFactor = 1;
    float minSizeGetDelta = 0;

    sf::Clock Clock;
    bool hasMovePrev = false;
    bool hasMove = false;

    while (win.isOpen()) {
        hasMovePrev = hasMove;
        hasMove = false;

        float camcos = cos(-_cameraAngle) * _camScale/100;
        float camsin = sin(-_cameraAngle) * _camScale/100;

        sf::Keyboard::isKeyPressed(sf::Keyboard::W) ? _cameraOffset[0] -= camsin, _cameraOffset[1] += camcos, hasMove = true : 0;
        sf::Keyboard::isKeyPressed(sf::Keyboard::S) ? _cameraOffset[0] += camsin, _cameraOffset[1] -= camcos, hasMove = true : 0;
        sf::Keyboard::isKeyPressed(sf::Keyboard::A) ? _cameraOffset[0] += camcos, _cameraOffset[1] += camsin, hasMove = true : 0;
        sf::Keyboard::isKeyPressed(sf::Keyboard::D) ? _cameraOffset[0] -= camcos, _cameraOffset[1] -= camsin, hasMove = true : 0;

        sf::Keyboard::isKeyPressed(sf::Keyboard::R) ? camPitch += 0.01, hasMove = true : 0;
        sf::Keyboard::isKeyPressed(sf::Keyboard::F) ? camPitch -= 0.01, hasMove = true : 0;

        sf::Keyboard::isKeyPressed(sf::Keyboard::Z) ? _cameraAngle += 0.01, hasMove = true : 0;
        sf::Keyboard::isKeyPressed(sf::Keyboard::X) ? _cameraAngle -= 0.01, hasMove = true : 0;
        
        sf::Keyboard::isKeyPressed(sf::Keyboard::Q) ? _camScale *= 1.02, hasMove = true : 0;
        sf::Keyboard::isKeyPressed(sf::Keyboard::E) ? _camScale *= 0.98, hasMove = true : 0;

        Vector2f mousePos;
        {
            sf::Vector2i vec = sf::Mouse::getPosition(win);
            mousePos = {(float)vec.x, (float)vec.y};
        }
        // if (!(mousePosPrev == mousePos))
        //     hasMove = true;
        // Vector2f mousePosCenterRelative = mousePos-winSize.cast<float>()/2;
        // mousePosPrev = mousePos;


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
                        // bool keep = false;
                        for (std::string const &s : shape._labels) {
                            std::string key = s.substr(0, s.rfind('='));
                            std::string value = s.substr(s.rfind('=')+1);
                            goto keepIt;
                            for (QuickMatch const &ss : toKeep) {
                                if (key == ss.key && (ss.value == 0 || value == ss.value)) {
                                    // keep = true;
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
                    std::ifstream file("../BoatFight/map.bmap", std::ifstream::binary);

                    ByteObject obj;
                    file >> obj;

                    obj >> segs;

                    // update bounding box
                    boundingBox = {FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
                    for (auto &shape : segs) {
                        shape._boundingBox = {FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
                        for (auto &p : shape._points) {
                            boundingBox._p[0] = std::min(p[0], boundingBox._p[0]);
                            boundingBox._p[1] = std::min(p[1], boundingBox._p[1]);
                            boundingBox._d[0] = std::max(p[0], boundingBox._d[0]);
                            boundingBox._d[1] = std::max(p[1], boundingBox._d[1]);

                            shape._boundingBox._p[0] = std::min((float)p[0], shape._boundingBox._p[0]);
                            shape._boundingBox._p[1] = std::min((float)p[1], shape._boundingBox._p[1]);
                            shape._boundingBox._d[0] = std::max((float)p[0], shape._boundingBox._d[0]);
                            shape._boundingBox._d[1] = std::max((float)p[1], shape._boundingBox._d[1]);
                        }
                        shape._boundingBox._d -= shape._boundingBox._p;
                    }
                    boundingBox._d -= boundingBox._p;
                    center = boundingBox._p+boundingBox._d/2;

                    treeZone = std::make_unique<UniTreeZone<float, StreeShape, 2>>(Zone<float, 2>(boundingBox.cast<float>()));


                    uint allocSize = segs.size();
                    std::cout << DEBUGVAR(allocSize) << std::endl;
                    uint currentSize = 0;
                    std::cout << "building map..." << std::endl;
                    srand(14);
                    for (StreeShape &shape : segs) {
                        treeZone->addData(Zone<float, 2>(shape._boundingBox), &shape);
                        if (!(currentSize%1000))
                            printProgress(currentSize/(float)allocSize);
                        currentSize += 1;
                    }
                    printProgress(1);
                }
            }
        }



        Mat3 matrix;


        matrix.tz(_camScale);
        matrix.rx(camPitch);
        matrix.rz(_cameraAngle);
        
        matrix.tx(_cameraOffset[0]);
        matrix.ty(_cameraOffset[1]);


        // matrix.scale(_camScale);
        // matrix.ttx(mousePosCenterRelative[0]*0.5);
        // matrix.tty(mousePosCenterRelative[1]*0.5);
        // matrix.rrz(_cameraAngle);




        _matWorld = matrix;
        _matWorldInv = matrix.inv();

        // sf::CircleShape ship;
        // ship.setFillColor(sf::Color::Red);
        // ship.setRadius(_camScale*1);
        // ship.setPosition(-mousePosCenterRelative[0]*0.8+winSize[0]/2, -mousePosCenterRelative[1]*0.8+winSize[1]/2);

        // win.draw(ship);

        Vector3f cameraPos = (_matWorldInv * Vector3f{0, 0, 0});
        Vector2f scale = winSizef/winSizef[0]*0.5;
        Vector3f p1 = (_matWorldInv * Vector3f{-scale[0], -scale[1], 1});
        Vector3f p2 = (_matWorldInv * Vector3f{-scale[0],  scale[1], 1});
        Vector3f p3 = (_matWorldInv * Vector3f{ scale[0], -scale[1], 1});
        Vector3f p4 = (_matWorldInv * Vector3f{ scale[0],  scale[1], 1});
        Vector3f p1d = p1 - cameraPos;
        Vector3f p2d = p2 - cameraPos;
        Vector3f p3d = p3 - cameraPos;
        Vector3f p4d = p4 - cameraPos;
        // std::cout << DEBUGVAR(cameraPos) << std::endl;
        if (cameraPos[2] < 0) {
            p1d[2] = std::max(p1d[2], (float)0.1);
            p2d[2] = std::max(p2d[2], (float)0.1);
            p3d[2] = std::max(p3d[2], (float)0.1);
            p4d[2] = std::max(p4d[2], (float)0.1);
        } else {
            p1d[2] = std::min(p1d[2], (float)-0.1);
            p2d[2] = std::min(p2d[2], (float)-0.1);
            p3d[2] = std::min(p3d[2], (float)-0.1);
            p4d[2] = std::min(p4d[2], (float)-0.1);
        }
        p1 += p1d * -(p1[2] / p1d[2]);
        p2 += p2d * -(p2[2] / p2d[2]);
        p3 += p3d * -(p3[2] / p3d[2]);
        p4 += p4d * -(p4[2] / p4d[2]);

        Segmentd cameraBoundingBox = {
            std::min(std::min(p1[0], p2[0]), std::min(p3[0], p4[0])),
            std::min(std::min(p1[1], p2[1]), std::min(p3[1], p4[1])),
            std::max(std::max(p1[0], p2[0]), std::max(p3[0], p4[0])),
            std::max(std::max(p1[1], p2[1]), std::max(p3[1], p4[1])),
        };
        cameraBoundingBox._d -= cameraBoundingBox._p;


        uniTreeBuffer->clear();
        treeZone->getColides(Zone<float, 2>(cameraBoundingBox.cast<float>()), minSizeGet, uniTreeBuffer);


        if (hasMove != hasMovePrev) {
            minSizeGetFactor = 1;
        }

        // std::cout << DEBUGVAR(uniTreeBuffer->size()) << std::endl;
        // std::cout << DEBUGVAR(_camScale) << std::endl;
        float time = Clock.getElapsedTime().asSeconds();
        Clock.restart();
        float maxDisplayShape = 1.0/ (hasMove ? targetFrameRateMonving : targetFrameRate);
        float delta = maxDisplayShape - time;
        // uint maxDisplayShape = 1024*4;
        // float delta = (float)maxDisplayShape-uniTreeBuffer->size();
        if (std::abs(delta) > maxDisplayShape*0.05) {
            // std::cout << DEBUGVAR(delta) << std::endl;
            float newMinSizeGetDelta = -(delta > 0 ? 1 : -1)*minSizeGetFactor;//*_camScale;
            if (minSizeGetDelta * newMinSizeGetDelta < 0)
                minSizeGetFactor = 0.1;
            else
                minSizeGetFactor *= 1.1;

            minSizeGet += newMinSizeGetDelta;
            minSizeGetDelta = newMinSizeGetDelta;
            // std::cout << DEBUGVAR(uniTreeBuffer->size()) << std::endl;
            // std::cout << DEBUGVAR(newMinSizeGetDelta) << std::endl;
            // std::cout << DEBUGVAR(minSizeGetFactor) << std::endl;
            // std::cout << DEBUGVAR(minSizeGet) << std::endl;
            // std::cout << DEBUGVAR(minSizeGet) << std::endl;
            if (minSizeGet < 0)
                minSizeGet = 0;
        }

        {
            uint maxTextDraw = 32;
            std::lock_guard<std::mutex> segsMutLock(segsMut);
            for (auto &obj : *uniTreeBuffer) {
                auto &shape = *obj->_data;
                sf::VertexArray lines(sf::LineStrip, shape._points.size());
                uint i = 0;  
                for (auto &p : shape._points) {
                    Vector3f vec3 = p.cast<float>();
                    Vector3f p1 = _matWorld * vec3;
                    if (p1[2] < 0)
                        continue;
                    p1[0] /= p1[2];
                    p1[1] /= p1[2];
                    p1 *= Vector3f{winSizef[0], winSizef[0], 0};
                    p1 += Vector3f{winSizef[0]/2, winSizef[1]/2, 0};
                    lines[i++] = sf::Vertex(sf::Vector2f(p1[0], p1[1]), shape._color);
                }
                lines.resize(i);
                win.draw(lines);

                if (displayLabel) {
                    Vector3f p1 = shape._points[0].cast<float>();
                    Vector3f vec3 = p1;
                    Vector3f p1Proj = _matWorld * vec3;
                    if (p1Proj[2] < 0)
                        continue;
                    p1Proj[0] /= p1Proj[2];
                    p1Proj[1] /= p1Proj[2];
                    p1Proj *= Vector3f{winSizef[0], winSizef[0], 0};
                    p1Proj += Vector3f{winSizef[0]/2, winSizef[1]/2, 0};

                    if (shape._labels.size() && maxTextDraw > 0 && 0 <= p1Proj[0] && p1Proj[0] < winSize[0] && 0 <= p1Proj[1] && p1Proj[1] < winSize[1]) {
                        --maxTextDraw;
                        text.setPosition(p1Proj[0], p1Proj[1]);
                        for (std::string const &string : shape._labels) {
                            text.setString(string);
                            win.draw(text);
                            text.move({0, textSize*(float)1.1});
                        }
                    }
                }
            }
        }
        win.display();
        // if (hasMove)
            win.clear();
    }
}

