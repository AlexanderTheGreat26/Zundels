#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <cmath>
#include <iterator>
#include <memory>
#include <algorithm>
#include <string>
#include <array>


typedef std::tuple<double, double, double> data_tuple;
typedef std::vector<data_tuple> frame;
typedef std::vector<frame> frames;

const long atoms_per_frame_count = 216;
const double box_size = 20.915;
const double Zundel_max_length = 2.47;


std::vector<data_tuple> coordinates_read (const std::string & name);

void shoving (std::vector<data_tuple> & coordinates, const double & rib_length);

frames oxygen_frames (std::vector<data_tuple> & coordinates, const long & atoms_per_frame);

void border_layer (frame & atoms);

int Zundel_count (frame & oxygens);

void data_file_creation (const std::string & name, std::vector<std::pair<int, int>> & data);


int main() {
    std::vector<data_tuple> centroid = std::move(coordinates_read("centroid.txt"));

    //shoving(centroid, box_size);



    frames oxygens = oxygen_frames(centroid, atoms_per_frame_count);

    std::vector<std::pair<int, int>> count_vs_frame;


    for (int i = 0; i < oxygens.size(); ++i) {
        border_layer(oxygens[i]);
        count_vs_frame.emplace_back(i, Zundel_count(oxygens[i]));
    }

    data_file_creation("Zundels_per_frame", count_vs_frame);

    return 0;
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


template<typename T, size_t... Is>
std::string tuple_to_string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple_to_string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple_to_string_impl(t, std::make_index_sequence<size>{});
}


void data_file_creation (const std::string & name, std::vector<std::pair<int, int>> & data) {
    std::ofstream fout;
    fout.open(name, std::ios::trunc);
    for (auto & i : data)
        fout << tuple_to_string(i) << '\n';
    fout.close();
}


bool default_area (const double & coordinate, const double & default_area_size) {
    return (coordinate <= default_area_size && coordinate >= 0);
}


template<typename T, size_t... Is>
bool belong_default_area_impl (T const& t, std::index_sequence<Is...>) {
    return (default_area(std::get<Is>(t), box_size) & ...);
}

template <class Tuple>
bool belong_default_area (const Tuple& t) { // Include constants as parameters later;
    constexpr auto size = std::tuple_size<Tuple>{};
    return belong_default_area_impl(t, std::make_index_sequence<size>{});
}


template<typename T, size_t... Is>
double distance_impl (T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return (std::sqrt((std::pow(std::get<Is>(t) - std::get<Is>(t1), 2) + ...)));
}

template <class Tuple>
double distance (const Tuple & t, const Tuple & t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return distance_impl(t, t1, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


int Zundel_count (frame & oxygens) {

    int Z = 0;
    std::vector<int> excludes;

    for (int i = 0; i < oxygens.size(); ++i) {

        if (!belong_default_area(oxygens[i]))
            continue;

        for (int j = 0; j < oxygens.size() && j != i; ++j) {

            if (std::any_of(excludes.begin(), excludes.end(), [&](int k) { return k == j; }))
                continue;

            double dist = distance(oxygens[i], oxygens[j]);

            if (dist > 2.24 && dist <= 2.47) {
                excludes.emplace_back(j);
                ++Z;
            }

        }
    }
    return Z;
}






template<size_t Is = 0, typename... Tp>
void pbc_layer (std::tuple<Tp...>& q, const double & layer_size, const double & l) {

    if (std::get<Is>(q) <= layer_size)
        std::get<Is>(q) += l;
    else if (std::get<Is>(q) > layer_size)
        if (l - std::get<Is>(q) < layer_size)
            std::get<Is>(q) = std::get<Is>(q) - l;

    if constexpr(Is + 1 != sizeof...(Tp))
        pbc_layer<Is + 1>(q, layer_size, l);
}


bool deps (double coordinate, const double & rig_length, const double & section_size) {
    if (coordinate <= section_size) return true;
    if (coordinate > section_size) return (rig_length - coordinate < section_size);
}


template<typename T, size_t... Is>
bool significance_of_pbc_impl (T const& t, std::index_sequence<Is...>) {
    return (deps(std::get<Is>(t), box_size, Zundel_max_length) & ...);
}

template <class Tuple>
bool significance_of_pbc (const Tuple& t) { // Include constants as parameters later;
    constexpr auto size = std::tuple_size<Tuple>{};
    return significance_of_pbc_impl(t, std::make_index_sequence<size>{});
}


//double lammps_box_size (std::vector<double>) {

//}


void Zundel_layer_particle (data_tuple particle, std::vector<data_tuple> & particles) {
        //double extended_box_size;
        pbc_layer(particle, Zundel_max_length, box_size);
        particles.emplace_back(particle);
}


void border_layer (frame & atoms) {
    for (const auto & atom : atoms)
        if (significance_of_pbc(atom))
            Zundel_layer_particle(atom, atoms);
}


frames oxygen_frames (std::vector<data_tuple> & coordinates, const long & atoms_per_frame) {
    long frames_count = coordinates.size() / atoms_per_frame + 1;

    long k = 0;
    frames simulation_steps (frames_count);
    for (long i = 0; i < coordinates.size(); ++i) {
        simulation_steps[k].emplace_back(coordinates[i]);
        if (i % (atoms_per_frame) == 0) ++k;
    }

        return simulation_steps;
}


template<size_t Is = 0, typename... Tp>
void fix_pbc (std::tuple<Tp...>& q, const double & l) {
    if (std::get<Is>(q) < 0) std::get<Is>(q) = l - std::get<Is>(q);
    if (std::get<Is>(q) > l) std::get<Is>(q) -= l;
    if constexpr(Is + 1 != sizeof...(Tp))
        fix_pbc<Is + 1>(q, l);
}


void shoving (std::vector<data_tuple> & coordinates, const double & rib_length) {
    for (auto & coordinate : coordinates)
        fix_pbc(coordinate, rib_length);
}


namespace std {
    istream& operator >> (istream& in, data_tuple & coordinates) {
        double first, second, third;
        in >> first >> second >> third;
        coordinates = {first, second, third};
        return in;
    }

    ostream& operator << (ostream& out, const data_tuple & coordinates) {
        auto [first, second, third] = coordinates;
        out << first << ' ' << second << ' ' << third << ' ';
        return out;
    }
}


// Read data from columns in text-file.
std::vector<data_tuple> coordinates_read (const std::string & name) {
    std::ifstream fin(name);
    if (!fin.is_open()) throw std::runtime_error("Error opening file.");
    std::vector<data_tuple> tuples_vector;
    copy(std::istream_iterator<data_tuple> {fin},
         std::istream_iterator<data_tuple> {},
         back_inserter(tuples_vector));
    //copy(tuples_vector.begin(), tuples_vector.end(), std::ostream_iterator<data>(std::cout, "\n"));
    return tuples_vector;
}

