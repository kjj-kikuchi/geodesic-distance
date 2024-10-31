// 測地距離の計算

#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>
#include <map>
#include <chrono>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/LU>

struct Halfedge
{
    int idx;
    int h_opp;
    std::vector<Eigen::Vector3i> &F;

    Halfedge(std::vector<Eigen::Vector3i> &faces) : F(faces) {}

    int face() const {
        return int(idx/3);
    }
    int h_next() const {
        if(idx % 3 == 2) return idx-2;
        else return idx+1;
    }
    int h_prev() const {
        if(idx % 3 == 0) return idx+2;
        else return idx-1;
    }
    int v_src() const {
        return F[face()][(idx+1)%3];
    }
    int v_tgt() const {
        return F[face()][(idx+2)%3];
    }
};

struct Mesh
{
    std::vector<Eigen::Vector3d> V;
    std::vector<Eigen::Vector3i> F;
    std::vector<Halfedge> hEList;
    std::vector<int> h_out;

    Mesh()
    {
        V.resize(0);
        F.resize(0);
    }

    void make_halfedgelist()
    {
        std::pair<int, int> key;
        std::pair<int, int> keyswap;
        std::map<std::pair<int, int>, int> map;
        for (int i=0; i<F.size(); i++) {
            for (int j=0; j<3; j++) {
                Halfedge h(F);
                h.idx = 3*i+j;

                key = std::make_pair(h.v_src(), h.v_tgt());
                keyswap = std::make_pair(h.v_tgt(),h.v_src());
                if (map.contains(keyswap)) {
                    h.h_opp = map.at(keyswap);
                    hEList[map.at(keyswap)].h_opp = 3*i+j;
                } else {
                    h.h_opp = -1;
                    map.emplace(key, 3*i+j);
                }

                hEList.push_back(h);
            }
        }
        // h_out を計算
        h_out.resize(V.size());
        for (int i=0; i<hEList.size(); i++) {
            //  h_out に境界半辺が保存されていない場合のみ更新
            if (hEList[ h_out[hEList[i].v_src()] ].h_opp != -1) {
                h_out[hEList[i].v_src()] = i;
            }
        }
        // 境界半辺の h_opp に次の境界半辺を保存する
        /*for (int i=0; i<h_out.size(); i++) {
            if(hEList[h_out[i]].h_opp == -1){
                hEList[h_out[i]].h_opp = -h_out[ hEList[h_out[i]].v_tgt() ] - 1;
                // h_i1.h_opp = -i2-1;
            }
        }*/
    }
    int h_ccw(int i) const
    {
        if(i < 0) return i;
        return hEList[ hEList[i].h_prev() ].h_opp;
    }
    int h_cw(int i) const
    {
        if (i < 0) return i;
        else if (hEList[i].h_opp < 0) return hEList[i].h_opp;
        return hEList[ hEList[i].h_opp ].h_next();
    }

    Eigen::Vector3d h_vec(int i) const
    {
        return V[hEList[i].v_tgt()] - V[hEList[i].v_src()];
    }

    double get_triangle_area(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2) const
    {
        return ((v1-v0).cross(v2-v0)).norm() / 2.0;
    }

    double get_triangle_area(Eigen::Vector3d e0, Eigen::Vector3d e1) const
    {
        return (e0.cross(e1)).norm() / 2.0;
    }
};

class GeodesicDistance
{
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using SparseVector = Eigen::SparseVector<double>;
    using SparseSolver = Eigen::SparseLU<SparseMatrix>;
    using Triplet = Eigen::Triplet<double>;

//    std::vector<Eigen::Vector3d> const* V;
//    std::vector<Eigen::Vector3i> const* F;
//    std::vector<Halfedge> hEList;
//    std::vector<int> h_out;
    SparseMatrix L_cotan;
    SparseMatrix A;
    SparseMatrix Laplacian;
    Eigen::VectorXd u;
    std::vector<Eigen::Vector3d> grad_u;
    Eigen::VectorXd divX;

    double time;

public:
    Eigen::VectorXd geodist;

//    GeodesicDistance(Mesh const* mesh) {
//        V = &(mesh->V);
//        F = &(mesh->F);
//    }

    void get_time(Mesh const& mesh)
    {
        double h = 0;
        for (auto& he : mesh.hEList) {
            h += (mesh.V[he.v_tgt()] - mesh.V[he.v_src()]).norm();
        }
        h /= mesh.hEList.size();
        time = h * h;
    }

    void make_cotan_laplacian(Mesh const& mesh)
    {
        L_cotan.resize(mesh.V.size(), mesh.V.size());
        std::vector<Triplet> triplets;

        for (int i = 0; i < mesh.V.size(); i++)
        {
            int h_temp = mesh.h_out[i];
            int h_end = h_temp;
            double sumW = 0;

            do {
                double cotanW;
                Eigen::Vector3d a1 = mesh.V[mesh.hEList[h_temp].v_src()] - mesh.V[mesh.hEList[mesh.h_cw(h_temp)].v_tgt()];
                Eigen::Vector3d a2 = mesh.V[mesh.hEList[h_temp].v_tgt()] - mesh.V[mesh.hEList[mesh.h_cw(h_temp)].v_tgt()];
                Eigen::Vector3d b1 = mesh.V[mesh.hEList[h_temp].v_tgt()] - mesh.V[mesh.hEList[mesh.hEList[h_temp].h_next()].v_tgt()];
                Eigen::Vector3d b2 = mesh.V[mesh.hEList[h_temp].v_src()] - mesh.V[mesh.hEList[mesh.hEList[h_temp].h_next()].v_tgt()];
                double alpha = acos( std::clamp( a1.dot(a2) / (a1.norm()*a2.norm()), -1.0, 1.0) );
                double beta = acos( std::clamp( b1.dot(b2) / (b1.norm()*b2.norm()), -1.0, 1.0) );

                cotanW = (1/tan(alpha) + 1/tan(beta)) / 2.0;
                if(isnan(cotanW)) std::cout << "error : cotan weight is not correct.\n";
                triplets.push_back(Triplet(i, mesh.hEList[h_temp].v_tgt(), cotanW));
                sumW += cotanW;

                h_temp = mesh.h_ccw(h_temp);
            } while (h_temp != h_end && h_temp >= 0);

            triplets.push_back(Triplet(i, i, -sumW));
        }

        L_cotan.setFromTriplets(triplets.begin(), triplets.end());
    }

    void make_mass_matrix(Mesh const& mesh)
    {
        A.resize(mesh.V.size(), mesh.V.size());
        std::vector<Triplet> triplets;

        for (int i = 0; i < mesh.V.size(); i++)
        {
            int h_temp = mesh.h_out[i];
            int h_end = h_temp;
            double area = 0;

            do {
                area += mesh.get_triangle_area(mesh.h_vec(h_temp), mesh.h_vec(mesh.h_ccw(h_temp))) / 3.0;

                h_temp = mesh.h_ccw(h_temp);
            } while (h_temp != h_end && h_temp >= 0);

            triplets.push_back(Triplet(i, i, area));
        }

        A.setFromTriplets(triplets.begin(), triplets.end());
    }

    void get_gradient(Mesh const& mesh)
    {
        grad_u.resize(mesh.F.size());

        for (int i = 0; i < mesh.F.size(); i++) {
            Eigen::Vector3i f = mesh.F[i];
            Eigen::Vector3d n = ( (mesh.V[f(1)] - mesh.V[f(0)]).cross(mesh.V[f(2)] - mesh.V[f(0)]) ).normalized();
            grad_u[i] =
            (u(f(0)) * ( n.cross(mesh.V[f(2)] - mesh.V[f(1)]) ) +
             u(f(1)) * ( n.cross(mesh.V[f(0)] - mesh.V[f(2)]) ) +
             u(f(2)) * ( n.cross(mesh.V[f(1)] - mesh.V[f(0)]) ))
            / (2 * mesh.get_triangle_area(mesh.V[f(0)], mesh.V[f(1)], mesh.V[f(2)]));
        }
    }

    void get_divergence(Mesh const& mesh)
    {
        divX.resize(mesh.V.size());

        for (int i = 0; i < mesh.V.size(); i++)
        {
            int h_temp = mesh.h_out[i];
            int h_end = h_temp;

            do {
                Eigen::Vector3d e1 = mesh.h_vec(h_temp);
                Eigen::Vector3d e2 = mesh.h_vec(mesh.h_ccw(h_temp));
                Eigen::Vector3d Xj = grad_u[mesh.hEList[h_temp].face()].normalized();
                double theta1 = acos( std::clamp((-e2).dot(e1-e2) / ((-e2).norm() * (e1-e2).norm()), -1.0, 1.0) );
                double theta2 = acos( std::clamp((-e1).dot(e2-e1) / ((-e1).norm() * (e2-e1).norm()), -1.0, 1.0) );

                divX(i) += 1/tan(theta1) * (e1.dot(Xj)) + 1/tan(theta2) * (e2.dot(Xj));

                h_temp = mesh.h_ccw(h_temp);
            } while (h_temp != h_end && h_temp >= 0);

            divX(i) /= 2.0;
            //if(isnan(divX(i))) std::cout << "error : ∇・Xi is nan.\n";
        }
    }

    void compute(Mesh const& mesh, int src_point)
    {
        // 1.
        get_time(mesh);
        make_cotan_laplacian(mesh);
        make_mass_matrix(mesh);
        Laplacian = A - time * L_cotan;

        SparseSolver solver;
        solver.compute(Laplacian);
        if (solver.info() != Eigen::Success) {
            std::cout << "decomposition failed\n";
        }

        Eigen::VectorXd delta(mesh.V.size());
        for (int i = 0; i < mesh.V.size(); i++){
            if (i == src_point) delta(i) = 1;
            else delta(i) = 0;
        }
        u = solver.solve(delta);
        if (u.hasNaN()) std::cout << "u has NaN\n";

        // 2.
        get_gradient(mesh);
        get_divergence(mesh);

        // 3.
        SparseSolver poisson_solver;
        poisson_solver.compute(L_cotan);
        if (solver.info() != Eigen::Success) {
            std::cout << "decomposition failed\n";
        }
        geodist = poisson_solver.solve(divX);
        if (u.hasNaN()) std::cout << "φ has NaN\n";
    }

};

void read_file(std::string const& filename, Mesh& mesh) {
    std::ifstream ifs(filename);
    if (ifs.fail()){
        std::cerr << "Failed to open file." << "\n";
        std::exit(1);
    }
    std::string line;
    if (std::filesystem::path(filename).extension() == ".obj") {
        while (std::getline(ifs, line)){
            if (line.empty() != 0 || line[0] == '#') {
                // Do nothing.
            } else if (line[0] == 'v') {
                Eigen::Vector3d v;
                std::istringstream string_in{line.substr(1)};
                string_in >> v(0) >> v(1) >> v(2);
                mesh.V.push_back(v);
            } else if (line[0] == 'f') {
                Eigen::Vector3i f;
                std::istringstream string_in{line.substr(1)};
                string_in >> f(0) >> f(1) >> f(2);
                f = f - Eigen::Vector3i{1, 1, 1};
                mesh.F.push_back(f);
            }
        }
    } else if (std::filesystem::path(filename).extension() == ".off") {
        std::getline(ifs, line);    // skip
        std::getline(ifs, line);
        int vsize, fsize;
        sscanf(line.data(), "%d %d", &vsize, &fsize);

        for (int i=0; i<vsize; i++) {
            std::getline(ifs, line);
            Eigen::Vector3d v;
            std::istringstream string_in{line.substr(1)};
            string_in >> v(0) >> v(1) >> v(2);
            mesh.V.push_back(v);
        }
        for (int i=0; i<fsize; i++) {
            std::getline(ifs, line);
            Eigen::Vector3i f;
            std::istringstream string_in{line.substr(1)};
            string_in >> f(0) >> f(1) >> f(2);
            mesh.F.push_back(f);
        }
    }
}

void export_vtk(std::string name, Mesh const& mesh, Eigen::VectorXd const& func)
{
    std::ofstream of;
    name.resize(name.length() - 4);
    std::string filename = name + "_geodesic.vtk";
    of.open(filename, std::ios::out);
    of << "# vtk DataFile Version 3.0\n" << "vtk output\n" << "ASCII\n" << "DATASET POLYDATA\n";
    of << "POINTS " << mesh.V.size() << " float\n";
    for(auto& v : mesh.V){
        of << v(0) << " " << v(1) << " " << v(2) << std::endl;
    }
    of << "POLYGONS " << mesh.F.size() << " " << 4*mesh.F.size() << std::endl;
    for(auto& f: mesh.F){
        of << "3 " << f(0) << " " << f(1) << " " << f(2) << std::endl;
    }
    of << "POINT_DATA " << mesh.V.size() << std::endl;
    of << "SCALARS geodesic_distance float\nLOOKUP_TABLE default\n";
    for (int i = 0; i < mesh.V.size(); i++) {
        of << func(i) << std::endl;
    }
    of.close();
}

void export_source_point(Eigen::Vector3d src_p)
{
    std::ofstream of;
    std::string filename = "source_point.obj";
    of.open(filename, std::ios::out);
    of << "v " << src_p(0) << " " << src_p(1) << " " << src_p(2) << std::endl;
    of << "p 1" << std::endl;
    of.close();
}

// =====================================================================
// main
// =====================================================================

int main(int argc, char *argv[]){

    // ソースの頂点を入力
    int src_point;
    std::cout << "Enter source vertex index : ";
    std::cin >> src_point;

    auto start = std::chrono::system_clock::now();

    // メッシュファイル読み込み
    std::string filename;
    if (argc != 2) {
        std::cout << "usage : " << argv[0] << "  filename" << std::endl;
        std::exit(1);
    }
    filename = std::string(argv[1]);


    Mesh mesh;
    read_file(filename, mesh);
    mesh.make_halfedgelist();

    // 測地距離計算
    GeodesicDistance geodesic_distance;
    geodesic_distance.compute(mesh, src_point);

    auto end = std::chrono::system_clock::now();

    // 出力
    export_vtk(filename, mesh, geodesic_distance.geodist);
    //export_source_point(mesh.V[src_point]);


    using namespace std::chrono_literals;
    std::cout << "φ(src) = " << geodesic_distance.geodist(src_point) << std::endl;
    std::cerr << "実行時間 : " << (end - start) / 1.0s << " 秒\n";
}

