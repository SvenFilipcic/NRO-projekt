#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <omp.h>

using namespace std;
int iteracije = 100;
int num_threads = 1; //st. thredov

int main() {
    // datoteka
    string filename = "C:\\Users\\filip\\Desktop\\projekt_2024\\projekt_2024\\primer3mreza.txt";

    // vektorji
    vector<double> X;
    vector<double> Y;
    vector<vector<int>> celice;
    vector<vector<int>> vozlisca_robnih_pogojev;
    vector<int> tipi_robnih_pogojev; // 0 = temperatura, 1 = toplotni tok
    vector<double> vrednosti_robnih_pogojev; // Vrednosti robnih pogojev
    vector<double> vrednosti_prestopa_toplote; // Vrednosti toplotnega toka

    ifstream file(filename);

    // Branje števila vozlišè iz prve vrstice
    string prva_vrstica;
    getline(file, prva_vrstica);
    int n_vozlisc = stoi(prva_vrstica.substr(6));

    // Branje koordinat vozlišè
    for (int line = 0; line < n_vozlisc; line++) {
        string line_str;
        getline(file, line_str);
        line_str = line_str.replace(line_str.find(';'), 1, " ");
        line_str = line_str.replace(line_str.find(','), 1, " ");
        stringstream line_stream(line_str);
        string temp;
        double x, y;
        line_stream >> temp >> x >> y;
        X.push_back(x);
        Y.push_back(y);
    }

    // Preskoèitev prazne vrstice
    string prazna_vrstica;
    getline(file, prazna_vrstica);

    // Branje števila celic iz datoteke
    string cell_line;
    getline(file, cell_line);
    stringstream cells_string_stream(cell_line);
    string temp;
    int n_celice;
    cells_string_stream >> temp >> n_celice;

    // Branje vozlišè za posamezne celice
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
    }

    // Preskoèitev naslednje prazne vrstice
    string prazna_vrstica2;
    getline(file, prazna_vrstica2);

    // Branje števila in tipov robnih pogojev
    string robni_pogoji_line;
    getline(file, robni_pogoji_line);
    stringstream robni_pogoji_stream(robni_pogoji_line);
    int n_pogoji;
    robni_pogoji_stream >> temp >> temp >> n_pogoji;

    // Procesiranje robnih pogojev
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
        }
        else if (tip_pogoja == "toplotni") {
            tipi_robnih_pogojev.push_back(1);
            getline(file, line);
            stringstream vrstica_stream(line);
            double toplotni_tok;
            vrstica_stream >> temp >> temp >> toplotni_tok;
            vrednosti_robnih_pogojev.push_back(toplotni_tok);
            vrednosti_prestopa_toplote.push_back(-1);
        }

        getline(file, line);
        int st_vozlisc_v_robnem_pogoju = stoi(line);
        vector<int> vozlisca_v_robnem_pogoju;
        for (int vozl = 0; vozl < st_vozlisc_v_robnem_pogoju; vozl++) {
            getline(file, line);
            int id_vozlisce_v_robnem_pogoju = stoi(line);
            vozlisca_v_robnem_pogoju.push_back(id_vozlisce_v_robnem_pogoju);
        }
        vozlisca_robnih_pogojev.push_back(vozlisca_v_robnem_pogoju);

        getline(file, line);
    }

    // Definicija konstanti za izraèun
    double deltaX = 1;
    double deltaY = 1;
    double k = 24;

    // Inicializacija matrike 
    vector<vector<int>> sosednja_vozlisca(n_vozlisc, vector<int>(4, -1));
    double start_time_neighbors = omp_get_wtime();

    // Iskanje sosedov
#pragma omp parallel for num_threads(num_threads)
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

    double end_time_neighbors = omp_get_wtime();
    cout << "èas za sosede: " << (end_time_neighbors - start_time_neighbors) << "s\n";

 
    vector<vector<double>> A(n_vozlisc, vector<double>(n_vozlisc, 0));
    vector<double> b(n_vozlisc, 0);

    // RP in logika za A
    for (int node_id = 0; node_id < n_vozlisc; node_id++) {
        vector<int> sosedi = sosednja_vozlisca[node_id];
        int levi_sosed = sosedi[0];
        int spodnji_sosed = sosedi[1];
        int desni_sosed = sosedi[2];
        int zgornji_sosed = sosedi[3];

        if (levi_sosed != -1 && spodnji_sosed != -1 && desni_sosed != -1 && zgornji_sosed != -1) { //središèna
            A[node_id][levi_sosed] = 1;
            A[node_id][spodnji_sosed] = 1;
            A[node_id][desni_sosed] = 1;
            A[node_id][zgornji_sosed] = 1;
            A[node_id][node_id] = -4;
        }
        else {
            for (int robni_pogoj_id = 0; robni_pogoj_id < n_pogoji; robni_pogoj_id++) {
                vector<int> vozlisca_v_trenutnem_rp = vozlisca_robnih_pogojev[robni_pogoj_id];
                for (int id_vozlisce_rp : vozlisca_v_trenutnem_rp) {
                    if (node_id == id_vozlisce_rp) {
                        int tip_robnega_pogoja = tipi_robnih_pogojev[robni_pogoj_id];
                        double vrednost = vrednosti_robnih_pogojev[robni_pogoj_id];

                        if (tip_robnega_pogoja == 0) {
                            A[node_id][node_id] = 1;
                            b[node_id] = vrednost;
                        }
                        else if (tip_robnega_pogoja == 1) {
                            if (desni_sosed == -1) {
                                A[node_id][node_id] -= 4;
                                A[node_id][levi_sosed] += 2;
                                A[node_id][spodnji_sosed] += 1;
                                A[node_id][zgornji_sosed] += 1;
                                b[node_id] = 0;
                            }
                        }
                    }
                }
            }
        }
    }

    // Shranjevanje matrike A v datoteko
    ofstream matrixAFile("C:\\Users\\filip\\Desktop\\projekt_2024\\projekt_2024\\matrix_A.txt");
    for (const auto& row : A) {
        for (size_t i = 0; i < row.size(); i++) {
            matrixAFile << row[i];
            if (i < row.size() - 1) {
                matrixAFile << " ";
            }
        }
        matrixAFile << "\n";
    }
    matrixAFile.close();

    // Shranjevanje vektorja b v datoteko
    ofstream vectorBFile("C:\\Users\\filip\\Desktop\\projekt_2024\\projekt_2024\\vector_b.txt");
    for (const auto& val : b) {
        vectorBFile << val << "\n";
    }
    vectorBFile.close();


    // Zaèetna T
    vector<double> T(n_vozlisc, 150);
    double start_time_gauss_seidel = omp_get_wtime();

    //Gauss_seidel
#pragma omp parallel num_threads(num_threads)
    {
        for (int iitt = 0; iitt < iteracije; iitt++) {
#pragma omp for schedule(static)
            for (int jj = 0; jj < n_vozlisc; jj++) {
                double d = b[jj];
                for (int ii = 0; ii < n_vozlisc; ii++) {
                    if (jj != ii) {
                        d -= A[jj][ii] * T[ii];
                    }
                }
#pragma omp critical
                {
                    T[jj] = d / A[jj][jj];
                }
            }
#pragma omp barrier
        }
    }

    double end_time_gauss_seidel = omp_get_wtime();
    cout << "Time for Gauss-Seidel: " << (end_time_gauss_seidel - start_time_gauss_seidel) << "s\n";

    // VTK 
    ofstream output_file("rezultat_vtk3.vtk");
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

    return 0;
}
