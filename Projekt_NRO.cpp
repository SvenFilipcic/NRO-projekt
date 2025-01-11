
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>

using namespace std;

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    string filename = ("C:\\Users\\filip\\Desktop\\projekt_2024\\projekt_2024\\primer3mreza.txt");

    vector<double> X;
    vector<double> Y;
    vector<vector<int>> celice;
    vector<vector<int>> vozlisca_robnih_pogojev;
    vector<int> tipi_robnih_pogojev;
    vector<double> vrednosti_robnih_pogojev;
    vector<double> vrednosti_prestopa_toplote;


    ifstream file(filename);

    // Preverjanje, ali je datoteka uspešno odprta
    if (!file) { // `ifstream` objekt se preverja z operatorjem `!`
        cerr << "Napaka: Ne morem odpreti datoteke \"" << filename << "\"." << endl;
        return 1; // Zakljuèi program z napako
    }

    string prva_vrstica;
    getline(file, prva_vrstica);
    int n_vozlisc = stoi(prva_vrstica.substr(6));

    for (int line = 0; line < n_vozlisc; line++) {
        string line_str;
        getline(file, line_str);
        line_str = line_str.replace(line_str.find(';'), 1, " ");
        line_str = line_str.replace(line_str.find(','), 1, " ");
        stringstream line_stream(line_str);
        string temp;
        double x, y;
        line_stream >> temp >> x >> y;
        //std::cout << "x: " << x << " .. y: " << y << "\n";
        X.push_back(x);
        Y.push_back(y);
    }

    string prazna_vrstica;
    getline(file, prazna_vrstica);

    string cell_line;
    getline(file, cell_line);
    stringstream cells_string_stream(cell_line);
    string temp;
    int n_celice;
    cells_string_stream >> temp >> n_celice;
    //std::cout << "n_celice: " << n_celice << "\n";

    for (int line = 0; line < n_celice; line++) {
        string line_str;
        getline(file, line_str);
        for (char& c : line_str) {
            if (c == ',') {
                c = ' ';
            }
        }
        line_str = line_str.replace(line_str.find(';'), 1, " ");
        stringstream line_stream(line_str);
        string temp;
        int node1_id, node2_id, node3_id, node4_id;
        line_stream >> temp >> node1_id >> node2_id >> node3_id >> node4_id;
        vector<int> cell = { node1_id, node2_id, node3_id, node4_id };
        celice.push_back(cell);
        //std::cout << "temp: " << temp << ".. node1_id: " << node1_id << ".. node2_id: " << node2_id << ".. node3_id: " << node3_id << ".. node4_id: " << node4_id << "\n";
        //std::cout << "line_str: " << line_str << "\n";
    }

    string prazna_vrstica2;
    getline(file, prazna_vrstica2);

    string robni_pogoji_line;
    getline(file, robni_pogoji_line);
    stringstream robni_pogoji_stream(robni_pogoji_line);
    int n_pogoji;
    robni_pogoji_stream >> temp >> temp >> n_pogoji;

    //std::cout << "robni_pogoji_stream: " << robni_pogoji_line << "\n";

    for (int n = 0; n < n_pogoji; n++) {
        string line;
        getline(file, line);
        stringstream tip_pogoja_stream(line);
        string tip_pogoja;
        tip_pogoja_stream >> temp >> temp >> tip_pogoja;

        if (tip_pogoja == "temperatura") {
            tipi_robnih_pogojev.push_back(0);
            getline(file, line);
            stringstream vrstica_stream(line);
            double temperatura;
            vrstica_stream >> temp >> temperatura;
            vrednosti_robnih_pogojev.push_back(temperatura);
            vrednosti_prestopa_toplote.push_back(-1);
            //std::cout << "temperatura: " << temperatura << "\n";
        }
        else if (tip_pogoja == "toplotni") {
            tipi_robnih_pogojev.push_back(1);
            getline(file, line);
            stringstream vrstica_stream(line);
            double toplotni_tok;
            vrstica_stream >> temp >> temp >> toplotni_tok;
            vrednosti_robnih_pogojev.push_back(toplotni_tok);
            vrednosti_prestopa_toplote.push_back(-1);
            //std::cout << "toplotni tok: " << toplotni_tok << "\n";
        }
        else if (tip_pogoja == "prestop") {
            tipi_robnih_pogojev.push_back(2);
            getline(file, line);
            stringstream vrstica_stream(line);
            double temperatura_prestopa;
            vrstica_stream >> temp >> temperatura_prestopa;
            getline(file, line);
            stringstream vrstica_stream2(line);
            double koeficient_prestopa;
            vrstica_stream2 >> temp >> temp >> koeficient_prestopa;
            vrednosti_robnih_pogojev.push_back(temperatura_prestopa);
            vrednosti_prestopa_toplote.push_back(koeficient_prestopa);
            //std::cout << "temp. prestopa: " << temperatura_prestopa << " .. koef. prestopa: " << koeficient_prestopa << "\n";
        }

        getline(file, line);
        int st_vozlisc_v_robnem_pogoju = stoi(line);
        //std::cout << "line: " << line << "\n";
        vector<int> vozlisca_v_robnem_pogoju;
        for (int vozl = 0; vozl < st_vozlisc_v_robnem_pogoju; vozl++) {
            getline(file, line);
            int id_vozlisce_v_robnem_pogoju = stoi(line);
            //std::cout << "id vozlisce v rp: " << id_vozlisce_v_robnem_pogoju<< "\n";
            vozlisca_v_robnem_pogoju.push_back(id_vozlisce_v_robnem_pogoju);
        }
        vozlisca_robnih_pogojev.push_back(vozlisca_v_robnem_pogoju);

        getline(file, line);
    }

    double deltaX = 1.25;
    double deltaY = 1.25;
    double k = 24;

    vector<vector<int>> sosednja_vozlisca(n_vozlisc, vector<int>(4, -1));


    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        vector<int> node_i_neighbours(4, -1);

        for (int nd = 0; nd < n_celice; nd++) {
            vector<int> trenutna_celica = celice[nd];
            int vozlisce1 = trenutna_celica[0];
            int vozlisce2 = trenutna_celica[1];
            int vozlisce3 = trenutna_celica[2];
            int vozlisce4 = trenutna_celica[3];

            if (node_id == vozlisce1 || node_id == vozlisce2 || node_id == vozlisce3 || node_id == vozlisce4) {
                for (int vozl = 0; vozl < 4; vozl++) {

                    int sosednje_vozlisce = trenutna_celica[vozl];
                    // std::cout << "trenutna_celica[vozl]: " << trenutna_celica[vozl] << "... sosednje vozlisce: " << sosednje_vozlisce << "\n";

                    if (sosednje_vozlisce != node_id) {
                        double x_obravnavano_vozl = X[node_id];
                        double y_obravnavano_vozl = Y[node_id];
                        double x_sosed = X[sosednje_vozlisce];
                        double y_sosed = Y[sosednje_vozlisce];
                        if ((x_obravnavano_vozl - x_sosed) < 1e-9 && (x_obravnavano_vozl - x_sosed) > -1e-9) {
                            if ((y_obravnavano_vozl - y_sosed) > 0) {
                                node_i_neighbours[1] = sosednje_vozlisce;
                            }
                            else {
                                node_i_neighbours[3] = sosednje_vozlisce;
                            }
                        }
                        else if ((y_obravnavano_vozl - y_sosed) < 1e-9 && (y_obravnavano_vozl - y_sosed) > -1e-9) {
                            if ((x_obravnavano_vozl - x_sosed) > 0) {
                                node_i_neighbours[0] = sosednje_vozlisce;
                            }
                            else {
                                node_i_neighbours[2] = sosednje_vozlisce;
                            }
                        }
                    }
                }
            }
        }
        sosednja_vozlisca[node_id] = node_i_neighbours;
    }

    vector<vector<double>> A(n_vozlisc, vector<double>(n_vozlisc, 0));

    vector<double> b(n_vozlisc, 0);

    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        vector<int> sosedi = sosednja_vozlisca[node_id];
        int levi_sosed = sosedi[0];
        int spodnji_sosed = sosedi[1];
        int desni_sosed = sosedi[2];
        int zgornji_sosed = sosedi[3];

        if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1) {
            A[node_id][levi_sosed] = 1;
            A[node_id][spodnji_sosed] = 1;
            A[node_id][desni_sosed] = 1;
            A[node_id][zgornji_sosed] = 1;
            A[node_id][node_id] = -4;
        }
        else {
            for (int robni_pogoj_id = 0; robni_pogoj_id < n_pogoji; robni_pogoj_id++) {
                vector<int> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[robni_pogoj_id];
                for (int id_vozlisce_rp = 0; id_vozlisce_rp < vozlisca_v_trenutnem_rp.size(); id_vozlisce_rp++) {
                    int vozlisce_v_trenutnem_rp = vozlisca_v_trenutnem_rp[id_vozlisce_rp];
                    if (node_id == vozlisce_v_trenutnem_rp) {
                        int tip_robnega_pogoja = tipi_robnih_pogojev[robni_pogoj_id];
                        double vrednost = vrednosti_robnih_pogojev[robni_pogoj_id];
                        double vrednost_prestopa = vrednosti_prestopa_toplote[robni_pogoj_id];

                        if (tip_robnega_pogoja == 0) {
                            A[node_id][node_id] = 1;
                            b[node_id] = vrednost;
                        }
                        else if (tip_robnega_pogoja == 1) {
                            int stevilo_sosedov = 0;
                            for (int st = 0; st < 4; st++) {
                                if (sosedi[st] != -1) {
                                    stevilo_sosedov++;
                                }
                            }

                            if (stevilo_sosedov == 3) {
                                if (levi_sosed == -1) {
                                    A[node_id][node_id] -= 4;
                                    A[node_id][desni_sosed] += 2;
                                    A[node_id][spodnji_sosed] += 1;
                                    A[node_id][zgornji_sosed] += 1;
                                    b[node_id] = -2 * (vrednost * deltaX / k);
                                }
                                if (desni_sosed == -1) {
                                    A[node_id][node_id] -= 4;
                                    A[node_id][levi_sosed] += 2;
                                    A[node_id][spodnji_sosed] += 1;
                                    A[node_id][zgornji_sosed] += 1;
                                    b[node_id] = -2 * (vrednost * deltaX / k);
                                }
                                if (spodnji_sosed == -1) {
                                    A[node_id][node_id] -= 4;
                                    A[node_id][zgornji_sosed] += 2;
                                    A[node_id][levi_sosed] += 1;
                                    A[node_id][desni_sosed] += 1;
                                    b[node_id] = -2 * (vrednost * deltaX / k);
                                }
                                if (zgornji_sosed == -1) {
                                    A[node_id][node_id] -= 4;
                                    A[node_id][spodnji_sosed] += 2;
                                    A[node_id][levi_sosed] += 1;
                                    A[node_id][desni_sosed] += 1;
                                    b[node_id] = -2 * (vrednost * deltaX / k);
                                }
                            }
                        }
                        else if (tip_robnega_pogoja == 2) {
                            int stevilo_sosedov = 0;
                            for (int st = 0; st < 4; st++) {
                                if (sosedi[st] != -1) {
                                    stevilo_sosedov++;
                                }
                            }

                            if (stevilo_sosedov == 3) {
                                if (levi_sosed == -1) {
                                    A[node_id][node_id] -= 2 * (vrednost_prestopa * deltaX / k + 2);
                                    A[node_id][desni_sosed] += 2;
                                    A[node_id][spodnji_sosed] += 1;
                                    A[node_id][zgornji_sosed] += 1;
                                    b[node_id] -= 2 * vrednost_prestopa * deltaX * vrednost / k;
                                }
                                if (desni_sosed == -1) {
                                    A[node_id][node_id] -= 2 * (vrednost_prestopa * deltaX / k + 2);
                                    A[node_id][levi_sosed] += 2;
                                    A[node_id][spodnji_sosed] += 1;
                                    A[node_id][zgornji_sosed] += 1;
                                    b[node_id] -= 2 * vrednost_prestopa * deltaX * vrednost / k;
                                }
                                if (spodnji_sosed == -1) {
                                    A[node_id][node_id] -= 2 * (vrednost_prestopa * deltaX / k + 2);
                                    A[node_id][levi_sosed] += 1;
                                    A[node_id][desni_sosed] += 1;
                                    A[node_id][zgornji_sosed] += 2;
                                    b[node_id] = -2 * vrednost_prestopa * deltaX * vrednost / k;
                                }
                                if (zgornji_sosed == -1) {
                                    A[node_id][node_id] -= 2 * (vrednost_prestopa * deltaX / k + 2);
                                    A[node_id][levi_sosed] += 1;
                                    A[node_id][desni_sosed] += 1;
                                    A[node_id][spodnji_sosed] += 2;
                                    b[node_id] = -2 * vrednost_prestopa * deltaX * vrednost / k;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    vector<double> T(n_vozlisc, 100);

    for (int iitt = 0; iitt < 100; iitt++) {
        //std::cout << "Current loop: " << iitt << "\n";
        for (int jj = 0; jj < n_vozlisc; jj++) {
            double d = b[jj];

            for (int ii = 0; ii < n_vozlisc; ii++) {
                if (jj != ii) {
                    d -= A[jj][ii] * T[ii];
                }
                T[jj] = d / A[jj][jj];
            }
        }
    }

    ofstream output_file("rezultati.vtk");
    output_file << "# vtk DataFile Version 3.0\n";
    output_file << "Mesh_1\n";
    output_file << "ASCII\n";
    output_file << "DATASET UNSTRUCTURED_GRID\n";
    output_file << "POINTS " << n_vozlisc << " float\n";
    for (int koordinata_id = 0; koordinata_id < n_vozlisc; koordinata_id++) {
        output_file << X[koordinata_id] << " " << Y[koordinata_id] << " 0\n";
    }
    output_file << "\n";
    output_file << "CELLS " << n_celice << " " << n_celice * 5 << "\n";
    for (int celica_id = 0; celica_id < n_celice; celica_id++) {
        int vozl_id1 = celice[celica_id][0];
        int vozl_id2 = celice[celica_id][1];
        int vozl_id3 = celice[celica_id][2];
        int vozl_id4 = celice[celica_id][3];
        output_file << "4 " << vozl_id1 << " " << vozl_id2 << " " << vozl_id3 << " " << vozl_id4 << "\n";
    }
    output_file << "\n";
    output_file << "CELL_TYPES " << n_celice << "\n";
    for (int celica_id = 0; celica_id < n_celice; celica_id++) {
        output_file << "9\n";
    }
    output_file << "\n";
    output_file << "POINT_DATA " << n_vozlisc << "\n";
    output_file << "SCALARS Temperature float 1\n";
    output_file << "LOOKUP_TABLE default\n";
    for (int koordinata_id = 0; koordinata_id < n_vozlisc; koordinata_id++) {
        output_file << T[koordinata_id] << "\n";
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Duration: " << duration.count() << "s" << std::endl;
    return 0;
}
