#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <list>
#include <stack>

using namespace std;

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

const float earthradius = 3963.0;     // [miles]
const float distance_accuracy = 5.0;  // [miles]

const int national_minpop = 1000000;

const float national_dist = 150.0;    // [miles]
const float regional_dist = 100.0;    // [miles]

const float local_maxdist = 50.0;     // [miles]

const float plane_speed = 520.0;      // [mph]
const float truck_speed = 65.0;       // [mph]

enum city_t { LOCAL, METRO, REGIONAL, NATIONAL };
enum color_t { WHITE, BLACK };

class city {
  public:
    //getters
    string get_city_name() const { return city_name; };
    city_t get_type() const { return type; };
    float get_latitude() const { return latitude; };
    float get_longitude() const { return longitude; };
    int get_population() const { return population; };

    //setters
    void set_city_name(string city_name) { this->city_name = city_name; };
    void set_type(city_t type) { this->type = type; };
    void set_latitude(float latitude) { this->latitude = latitude; };
    void set_longitude(float longitude) { this->longitude = longitude;};
    void set_population(int population) { this->population = population; };

    friend ostream & operator<<(ostream &, const city &);
    bool operator<(const city &) const;

    static int w;

  private:
    string city_name;
    city_t type;
    float latitude;
    float longitude;
    int population; 
};

//2d matrix class
//stores distance and time table info
class matrix {
  public:
    //constructor
    matrix (int row) {
      v.resize(row);
      for (size_t i = 0; i < v.size(); i++) {
        v[i].resize(i + 1);
      }
    }
    //get data (distance or time)
    float get_data(int row, int col) const {
      if (col > row) { return v[col][row]; }
      else { return v[row][col]; }
    }
    //set data (dist or time)
    void set_data(int row, int col, float data) {
      if (col > row) { v[col][row] = data; }
      else { v[row][col] = data; }
    }

  private:
    vector<vector<float>> v;
};

int city::w = 0;

istream & operator>>(std::istream &in, city &c) { 
  string city_name, type;
  int population;
  float latitude, longitude;

  in >> city_name >> type >> latitude >> longitude >> population;

  //setting name, lat, long, and population
  c.set_city_name(city_name);
  //setting type
  if (type == "METRO") {
    c.set_type(METRO);
  }
  else {
    c.set_type(LOCAL);
  }
  c.set_latitude(static_cast<float>(latitude * DEG2RAD));
  c.set_longitude(static_cast<float>(longitude * DEG2RAD));
  c.set_population(population);

  return in;
};

string change_type(city_t type) {
  switch (type) {
    case LOCAL: return "LOCAL";
    case METRO: return "METRO";
    case REGIONAL: return "REGIONAL";
    case NATIONAL: return "NATIONAL";
    default: return "unknown";
  }
}

ostream & operator<<(ostream &out, const city &c) { 
  out << c.city_name << "  " << setfill(' ') << setw(8) << change_type(c.type) << "  " << setw(8) 
  << right << c.population << "  " << fixed << setprecision(2) << setw(7) << c.latitude * RAD2DEG 
  << "  " << setw(7) << fixed << setprecision(2) << c.longitude * RAD2DEG;

  return out;
};

//comparing population
bool city::operator<(const city &rhs) const {
  return rhs.population < population;
}

void create_vertex_table(char *fname, vector<city> &vertex_table) {
  ifstream fin;
  fin.open(fname);

  string line;
  //reading in contents from csv file
  while(getline(fin, line)) {
    int count = 0;
    for (char &c : line) {
      //appending city anme
      if (c == ',') {
        count++;
        if (count == 1) {
          c = '_';
        }
        if (count > 1) {
          c = ' ';
        }
      }
    }
    city newcity;
    istringstream sin(line);
    sin >> newcity;
    //created vertex table
    vertex_table.push_back(newcity);
  }
  //sort vertex table form largest population to smallest
  sort(vertex_table.begin(), vertex_table.end());
  fin.close();
}

void update_vertex_table(vector<city> &vertex_table, matrix &dist_table) {
  for (auto &city : vertex_table) {
    //updating metro cities to national or regional
    if (city.get_type() == METRO) {
      if (city.get_population() > national_minpop) {
        city.set_type(NATIONAL);
      }
      else {
        city.set_type(REGIONAL);
      }
    }
  }
  for (size_t i = 0; i < vertex_table.size(); i++) {
    for (size_t j = 0; j < vertex_table.size(); j++) {
      //updating national cities to regional
      if (vertex_table[i].get_type() == NATIONAL && vertex_table[j].get_type() == NATIONAL) {
        if (dist_table.get_data(i, j) < national_dist) {
          if (vertex_table[j].get_population() > vertex_table[i].get_population()) {
            vertex_table[i].set_type(REGIONAL);
          }
        }
      }
    }
  }
  for (size_t i = 0; i < vertex_table.size(); i++) {
    for (size_t j = 0; j < vertex_table.size(); j++) {
      //updating regional cities to local
      if (vertex_table[i].get_type() == REGIONAL && vertex_table[j].get_type() == REGIONAL) {
        if (dist_table.get_data(i, j) < regional_dist) {
          if (vertex_table[j].get_population() > vertex_table[i].get_population()) {
            vertex_table[i].set_type(LOCAL);
          }
        }
      }
    }
  }
}

void create_dist_table(vector<city> &vertex_table, matrix &dist_table) {
  for (size_t i = 0; i < vertex_table.size(); ++i) {
    for (size_t j = 0; j < vertex_table.size(); ++j) {
      //calculating distance using haversine formula
      float lat1 = vertex_table[i].get_latitude();
      float lat2 = vertex_table[j].get_latitude();

      float long1 = vertex_table[i].get_longitude();
      float long2 = vertex_table[j].get_longitude();

      float dlat = (lat2 - lat1);
      float dlong = (long2 - long1);

      float a = pow(sin(dlat / 2), 2) + pow(sin(dlong / 2), 2) * cos(lat1) * cos(lat2);
      float b = 2 * asin(sqrt(a));
      float dist = earthradius * b;
      //rounding to nearest 5 mile distance
      dist = distance_accuracy * round(dist/distance_accuracy);

      //setting distance
      dist_table.set_data(i, j, dist);
      dist_table.set_data(j, i, dist);
    }
  }
}

void create_time_table(vector<city> &vertex_table, matrix &dist_table, matrix &time_table) {
  for (size_t i = 0; i < vertex_table.size(); i++) {
    for (size_t j = 0; j < vertex_table.size(); j++) {
      //getting time using plane and truck speed
      if (vertex_table[i].get_type() == NATIONAL && vertex_table[j].get_type() == NATIONAL) {
        time_table.set_data(i, j, (dist_table.get_data(i, j) / plane_speed));
      }
      else {
        time_table.set_data(i, j, (dist_table.get_data(i, j) / truck_speed));
      }
    }
  }
}

struct compare {
  float distance;
  int index;

//comparing distances
  bool operator<(const compare &rhs) {
    return distance < rhs.distance;
  }
};

void create_edge_table(vector<city> &vertex_table, vector<set<int>> &edge_table, matrix &dist_table) {
  edge_table.resize(vertex_table.size());

  vector<int> national_cities;
  vector<int> regional_cities;
  vector<int> local_cities;

  //pushing cities onto corresponding vectors
  for (size_t i = 0; i < vertex_table.size(); i++) {
    if (vertex_table[i].get_type() == NATIONAL) {
      national_cities.push_back(i);
    }
    else if (vertex_table[i].get_type() == REGIONAL) {
      regional_cities.push_back(i);
    }
    else {
      local_cities.push_back(i);
    }
  }

  //connecting national cities to each other
  for (size_t i = 0; i < national_cities.size(); i++) {
    for (size_t j = 0; j < national_cities.size(); j++) {
      if (i != j) {
        edge_table[national_cities[i]].insert(national_cities[j]);
        edge_table[national_cities[j]].insert(national_cities[i]);
      }
    }
  }

  for (size_t i = 0; i < regional_cities.size(); i++) {
    vector<compare> non_local;
    for (size_t j = 0; j < vertex_table.size(); j++) {
      if (vertex_table[j].get_type() != LOCAL && static_cast<size_t>(regional_cities[i]) != j) {
        compare c;
        c.distance = dist_table.get_data(regional_cities[i], j);
        c.index = j;
        non_local.push_back(c);
      }
    }
    //getting three nearest non local cities
    partial_sort(non_local.begin(), non_local.begin() + 3, non_local.end());
    for (size_t j = 0; j < 3; j++) {
      edge_table[regional_cities[i]].insert(non_local[j].index);
      edge_table[non_local[j].index].insert(regional_cities[i]);
    }
  }

  for (size_t i = 0; i < local_cities.size(); i++) {
    vector<compare> non_local;
    for (size_t j = 0; j < vertex_table.size(); j++) {
      if (vertex_table[j].get_type() != LOCAL) {
        compare c;
        c.distance = dist_table.get_data(local_cities[i], j);
        c.index = j;
        non_local.push_back(c);
      }
      if (vertex_table[j].get_type() == LOCAL && static_cast<size_t>(local_cities[i]) != j) {
        if (dist_table.get_data(local_cities[i], j) < local_maxdist) {
          edge_table[local_cities[i]].insert(j);
          edge_table[j].insert(local_cities[i]);
        }
      }
    }
    //getting five nearest non local cities
    partial_sort(non_local.begin(), non_local.begin() + 5, non_local.end());
    for (size_t j = 0; j < 5; j++) {
      edge_table[local_cities[i]].insert(non_local[j].index);
      edge_table[non_local[j].index].insert(local_cities[i]);
    }
  }
}

//writing out vertex table
void write_vertex_table(vector<city> &vertex_table) {
  ofstream fout;
  fout.open("vertex_table.txt");

  //setting width based on city name
  for (const auto &city_obj : vertex_table) {
    city::w = max(city::w, static_cast<int>(city_obj.get_city_name().length()));
  }

  int index = 0;
  for (size_t i = 0; i < vertex_table.size(); i++) {
    fout << setw(4) << right << index << "  " << left << setw(city::w + 3) << setfill('.') << vertex_table[i] << endl; 
    index++;
  }
  fout.close();
}

//writing out distance table
void write_dist_table(vector<city> &vertex_table, matrix &dist_table) {
  ofstream fout;
  fout.open("dist_table.txt");

  //setting width based on city name
  for (const auto &city_obj : vertex_table) {
    city::w = max(city::w, static_cast<int>(city_obj.get_city_name().length()));
  }

  for (size_t i = 1; i <= vertex_table.size(); i++) {
    for (size_t j = 1; j <= i; j++) {
      if (i != j) {
        fout << setfill(' ') << setw(4) << right << i - 1 << "  " << left << setw(city::w + 3) << setfill('.') << vertex_table[i - 1].get_city_name() 
        << " to " << setw(city::w + 3) << setfill('.') << vertex_table[j - 1].get_city_name() << fixed << setprecision(1) << setw(8) << right << setfill(' ') << dist_table.get_data(i - 1, j - 1)<< " miles" << endl;
      }
      if (j == i - 1) {
        fout << endl;
      }
    }
  }
  fout.close();
}

//writing out time table
void write_time_table(vector<city> &vertex_table, matrix &time_table) {
  ofstream fout;
  fout.open("time_table.txt");

  //setting width based on city name
  for (const auto &city_obj : vertex_table) {
    city::w = max(city::w, static_cast<int>(city_obj.get_city_name().length()));
  }

  for (size_t i = 1; i <= vertex_table.size(); i++) {
    for (size_t j = 1; j <= i; j++) {
      if (i != j) {
        fout << setfill(' ') << setw(4) << right << i - 1 << "  " << left 
        << setw(city::w + 3) << setfill('.') << vertex_table[i - 1].get_city_name() 
        << " to " << setw(city::w + 3) << setfill('.') << vertex_table[j - 1].get_city_name() 
        << fixed << setprecision(1) << setw(8) << right << setfill(' ') << time_table.get_data(i - 1, j - 1) 
        << " hours" << endl;
      }
      if (j == i - 1) {
        fout << endl;
      }
    }
  }
  fout.close();
}

//writing out edge table
void write_edge_table(vector<city> &vertex_table, vector<set<int>> &edge_table, matrix &dist_table, matrix &time_table) {
  ofstream fout;
  fout.open("edge_table.txt");

  //setting width based on city name
  for (const auto &city_obj : vertex_table) {
    city::w = max(city::w, static_cast<int>(city_obj.get_city_name().length()));
  }

  for (size_t i = 0; i < edge_table.size(); i++) {
    fout << setw(4) << right << i << " " << vertex_table[i].get_city_name() << endl;
    for (auto &j : edge_table[i]) {
      fout << "  " << setw(4) << setfill(' ') << j << " " << setw(city::w + 3) 
      << setfill('.') << left << vertex_table[j].get_city_name() << "  "
      << setw(6) << setfill(' ') << fixed << right << setprecision(1) 
      << dist_table.get_data(i, j) << " miles " << fixed << setprecision(1) 
      << right << setw(4) << time_table.get_data(i, j) << " hours" << endl;
    }
    if (i != edge_table.size() -1) {
      fout << endl;
    }
  }
  fout.close();
}

void dijkstra_route(int city_from, int city_to, vector<city> &vertex_table, vector<set<int>> &edge_table, int mode, matrix &dist_table, matrix &time_table) {
  vector<color_t> vcolor;
  vcolor.assign(vertex_table.size(), WHITE);

  vector<float> vdist;
  vdist.assign(vertex_table.size(), numeric_limits<float>::max());
  vdist[city_from] = 0;

  vector<int> vlink;
  vlink.assign(vertex_table.size(), -1);
  vlink[city_from] = city_from;

  while (1) {
    int next_i = -1;
    float mindist = numeric_limits<float>::max();

    for (size_t i = 0; i < vcolor.size(); i++) {
      if (vcolor[i] == WHITE && mindist > vdist[i]) {
        next_i = i;
        mindist = vdist[i];
      }
    }

    int i = next_i;
    if (i == -1) { return; }

    vcolor[i] = BLACK;

    if (i == city_to) { break; }

    for (auto &j : edge_table[i]) {
      float wij;
      //checking mode
      if (mode == 1) {
        wij = dist_table.get_data(i, j);
      }
      else if (mode == 2) {
        wij = time_table.get_data(i, j);
      }

      if (vcolor[j] == BLACK) {
        continue;
      }
      if (vdist[j] > vdist[i] + wij) {
        vdist[j] = vdist[i] + wij;
        vlink[j] = i;
      }
    }
  }
  
  if (vdist[city_to] == numeric_limits<float>::max()) {
    cout << "No path found" << endl;
    return;
  }
  stack<int> S;

  for (int i = city_to; i != city_from; i = vlink[i]) {
    S.push(i);
  }
  S.push(city_from);

  for (const auto &city_obj : vertex_table) {
    city::w = max(city::w, static_cast<int>(city_obj.get_city_name().length()));
  }

  int i = S.top();
  float total_time = 0, total_distance = 0;

//printing out first city
  cout << setw(8) << right << fixed << setprecision(2) << total_distance << " miles " 
    << right << setw(5) << setfill(' ') << fixed << setprecision(2) << total_time << " hours  "
    << setw(3) << i << " " << setw(city::w + 3) << left << setfill('.') 
    << vertex_table[i].get_city_name() << "  " << setfill(' ') << right 
    << setw(8) << vertex_table[i].get_population() << endl;
  S.pop();

//printing out the rest of the cities
  while(!S.empty()) {
    int j = S.top();
    total_time = total_time + time_table.get_data(i, j);
    total_distance = total_distance + dist_table.get_data(i, j);

    cout << setw(8) << right << fixed << setprecision(2) << total_distance << " miles " 
    << right << setw(5) << setfill(' ') << fixed << setprecision(2) << total_time << " hours  "
    << setw(3) << j << " " << setw(city::w + 3) << left << setfill('.') 
    << vertex_table[j].get_city_name() << "  " << setfill(' ') << right 
    << setw(8) << vertex_table[j].get_population() << " " << right << setw(7) << dist_table.get_data(i, j) << " miles  " << time_table.get_data(i, j) << " hours" << endl;
    i = j;
    S.pop();
  }
}

int main(int argc, char *argv[]) {
  //error handling
  if (argc != 3) {
    cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
    return 1;
  }
  
  bool info = false;
  bool dist = false;
  bool time = false;
  // parse commandline arguments
  string a = argv[1];
  if (a == "-info") {
    info = true;
  }
  else if (a == "-dist") {
    dist = true;
  }
  else if (a == "-time") {
    time = true;
  }

//error handling
  if (!info && !dist && !time) {
    cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
    return 1;
  }

  vector<city> vertex_table;
  create_vertex_table(argv[2], vertex_table);

  matrix dist_table(vertex_table.size());
  create_dist_table(vertex_table, dist_table);

  update_vertex_table(vertex_table, dist_table);

  matrix time_table(vertex_table.size());
  create_time_table(vertex_table, dist_table, time_table);

  vector<set<int>> edge_table;
  create_edge_table(vertex_table, edge_table, dist_table);

//if -info, write to files
  if (info) {
    write_vertex_table(vertex_table);
    write_time_table(vertex_table, time_table);
    write_dist_table(vertex_table, dist_table);
    write_edge_table(vertex_table, edge_table, dist_table, time_table);
  }
  else {
    while (true) {
      string city_from_name, city_to_name;
      int city_from, city_to;
      bool found_from = false;
      bool found_to = false;
  
      //ask for input (city from and city to)
      cout << "Enter> ";
      cin >> city_from_name >> city_to_name;

      for (size_t i = 0; i < vertex_table.size(); i++) {
        //set the city from and city to
        if (city_from_name == vertex_table[i].get_city_name()) {
          city_from = i;
          found_from = true;
        }
        if (city_to_name == vertex_table[i].get_city_name()) {
          city_to = i;
          found_to = true;
        }
      }

      if (cin.eof()) {
        break;
      }
      // error checking if city is not found
      if (!found_from) {
        cout << city_from_name << ": prefix not found" << endl;
        cout << endl;
        continue;
      }
      if (!found_to) {
        cout << city_to_name << ": prefix not found" << endl;
        cout << endl;
        continue;
      }

      //call dijkstra function based on command line arguments
      if (dist && found_from && found_to) {
        int d = 1;
        dijkstra_route(city_from, city_to, vertex_table, edge_table, d, dist_table, time_table);
        cout << endl;
      }
      else if (time && found_from && found_to) {
        int t = 2;
        dijkstra_route(city_from, city_to, vertex_table, edge_table, t, dist_table, time_table);
        cout << endl;
      }
    }
    cout << endl;
  }
  return 0;
};