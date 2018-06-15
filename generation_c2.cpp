#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <utility>
#include <stack> 
#include <stdexcept>
#include <vector>
#include <fstream>
#include <giac/config.h>
#include <giac/gen.h>
#include <giac/unary.h>

#include <giac/giac.h>

typedef std::vector<std::pair<int,int>> _graph;
std::set<int> emptyset;
using namespace giac;
//using p_type = polynomial<integer,monomial<int>>;
//    p_type x{"x"},y{"y"},z{"z"};
//    std::cout<< x*y*z + (x-z)*(x+y) << '\n';
//
void initial_processing (std::vector<std::string>& all_lines, std::vector<std::string>& periods, std::vector<std::string>& unprocessed_edges){

    for(auto graph_info : all_lines){
        std::string delim = ":=";

        auto pos = graph_info.find(delim);
        std::string period = graph_info.substr(0,pos);
        graph_info.erase(0,delim.length() + period.length() + 2);
        //std::cout << period << '\n';

        auto second_pos = graph_info.find(']');
        std::string edges = graph_info.substr(0,second_pos);
        //std::cout << edges << '\n';

        periods.push_back(period);
        unprocessed_edges.push_back(edges);
    }

}

void reduced_mod_p (gen& poly, int p){
    if(poly.type != _POLY){
        std::cerr << "reduced mod p not type poly \n";
        exit(0);
    }

    auto polyptr = poly._POLYptr;
    for (auto& monomial : polyptr -> coord){
       monomial.value = monomial.value % p; 
    }
}

gen find_coeff_p (gen& poly, int p, vecteur exps){
    if(poly.type != _POLY){
        std::cerr << "need to be sparse poly type \n";
        exit(0);
    }

    index_t ind;
    for(auto var : exps){
        ind.push_back(p-1);
    }

    index_m i(ind);

    auto polyptr = poly._POLYptr;
    for (auto& monomial : polyptr -> coord){
        if(monomial.index == i){
            std::cout << monomial.index << '\n';
            return monomial.value;
        }
    }

    return 0;
}

//We construct the graph data structures for all graphs in the graphs file.
void populate_graphs(std::vector<std::string>& unprocessed_edges, std::vector<_graph>& graphs,
        std::map<_graph,std::set<int>>& vertex_set, bool decomplete = true){
    for(auto edges : unprocessed_edges){

        int decomplete_vertex;
        _graph  edges_for_graph;
        std::set<int> vertices;
        int first_vertex;
        int second_vertex;
        bool is_first = true;

        for(int i = 0; i < edges.length(); i++){
            std::string number;
            bool digit_found = false;
            while(isdigit(edges[i])){
                number.append(std::string(1,edges[i]));
                i++;
                digit_found = true;
            }
            if(digit_found){
                if(is_first){
                    first_vertex= stoi(number);
                    is_first = false;
                    vertices.insert(first_vertex);
                }
                else{
                    is_first = true;
                    second_vertex = stoi(number);

                    decomplete_vertex = second_vertex;
                    vertices.insert(second_vertex);
                    auto edge_pair = std::make_pair(first_vertex,second_vertex);
                    edges_for_graph.push_back(edge_pair);

                  //  std::cout << first_edge << std::endl;
                }

                digit_found = false;
            }
        }
        //Decomplete the graph by removing a vertex and all associated edges.
        if(decomplete){
            _graph edges_for_graph_decomplete; 

            for(auto edge : edges_for_graph){
                if(edge.first == decomplete_vertex || edge.second == decomplete_vertex){
                    continue;
                }
                edges_for_graph_decomplete.push_back(edge);
            }

            vertices.erase(decomplete_vertex);
            graphs.push_back(edges_for_graph_decomplete);
            vertex_set.insert(std::pair<_graph, std::set<int>> (edges_for_graph_decomplete, vertices));
        }

        else{
            graphs.push_back(edges_for_graph);
            vertex_set.insert(std::pair<_graph, std::set<int>> (edges_for_graph, vertices));
        }

    }
}

void get_incidence_matrices(std::vector<std::vector<std::vector<int>>>& incidence_matrices,
        std::vector<_graph>& graphs, std::map<_graph,std::set<int>>& vertex_set){


    for(auto graph : graphs){
        std::vector<std::vector<int>> incidence_matrix;
        int column_length = vertex_set[graph].size();
        for(auto edges : graph){
            std::vector<int> edge_column;
           
            //minus 1 because we need to remove a column (or row) or vertices
            for (int i = 1; i < column_length; i++){
            //for (int i = 0; i < column_length; i++){

                int to_pushback = 0;
                //std::cout << edges.first << ' ' << edges.second << '\n';
                
                if(i == (edges.first - 1)){
                    to_pushback = -1;
                }

                else if (i == (edges.second - 1)){
                    to_pushback = 1;
                }

                edge_column.push_back(to_pushback);
            }
            incidence_matrix.push_back(edge_column);
        }
        incidence_matrices.push_back(incidence_matrix);

        //for(auto col : incidence_matrix){
        //    for(auto ent : col){
        //        std::cout<< ent << ",";
        //    }
        //    std::cout<< ';' << std::endl;
        //}
    }
}

//expressions : empty vecteur to be populated with the correct edge variables.
//inc_matrix : the incidence matrices. the rows (first index) are the edges, the columns are the vertices
//I : set I with integer edges. 
//J : set J with integer edges. 
//K : set K with integer edges.
vecteur
compute_kirchoff_matrix(vecteur& expressions, const std::vector<std::vector<int>>& inc_matrix, const std::set<int>& I, 
        const std::set<int>& J, const std::set<int>& K){

    if(expressions.size() != 0){
        std::cout << "Have an empty expression before computing kirchoff matrix " << '\n';
        exit(0);
    }
    if(I.size() != J.size()){
        std::cout << " I != J size, exiting." << '\n';
        std::cout << I.size() << J.size() << '\n';
        exit(0);
    }

    int num_edge = inc_matrix.size();
    int num_vertex = inc_matrix[0].size();
    vecteur kirchoff_matrix;
    vecteur exps;
    vecteur exps_to_return;
    
    //populate the expressions vector for determining the poly-variables.
    for(int i = 0; i < num_edge; i++){

        gen edgesym;
        if(I.find(i+1) != I.end() || J.find(i+1) != J.end() || K.find(i+1) != K.end()){
            //Put a dummy vector so we can discern its type.
            exps.push_back(makevecteur(1));
        }

        else{
            std::string edge_var = "a_" + std::to_string(i);
            edgesym = gen(std::string(edge_var),giac::context0);
            exps.push_back(edgesym);
            exps_to_return.push_back(edgesym);

            //Bernard said to do this after parsing a variable.
            eval(edgesym,1,giac::context0);
            //std::cout << edgesym.type << "TYPE \n";
        }
    }

    expressions = exps_to_return;

    for(int i = 0; i < num_edge; i++){

        vecteur kirchoff_row;

        if(I.find(i+1) != I.end()){

//            for(int z = 0; z < num_edge + num_vertex; z++){
//                gen thisgen = 0;
//                kirchoff_row.push_back(0);
//            }
//
//            kirchoff_matrix.push_back(kirchoff_row);
            
            //std::cout << kirchoff_matrix << "KIRCH KAT \n";
            //std::cout << "I found " << i+1 << '\n';
            continue;
        }

        gen edgesym = exps[i];
        gen poly_edgesym;

        //TODO This is kind of messy.
        if(edgesym.type != _VECT){
            poly_edgesym = sym2r(edgesym,exps_to_return,giac::context0);
        }

        else{
            poly_edgesym = 0;
        }
            
        for(int j = 0; j < num_edge;j++){

            if(J.find(j+1) != J.end()){
//                kirchoff_row.push_back(0);
                //std::cout << "J found " << j+1 << '\n';
                continue;
            }

            if(j == i){
                if(K.find(j+1) != K.end()){
                    //std::cout << "K found " << j+1 << '\n';
                    kirchoff_row.push_back(0);
                }

                else{
                    //kirchoff_row.push_back(poly_edgesym);
                    kirchoff_row.push_back(edgesym);
                    //kirchoff_row.push_back(1);
                }
            }

            else{
               //kirchoff_row[j] = sym2r(0,exps,giac::context0);
               //kirchoff_row.push_back(sym2r(0,exps,giac::context0));
               kirchoff_row.push_back(0);
            }
        }

        for(int k = num_edge; k < num_edge + num_vertex; k++){
            //kirchoff_row[k] = sym2r(inc_matrix[i][k - num_edge],exps,giac::context0);
            kirchoff_row.push_back(inc_matrix[i][k - num_edge]);
        }

//        for(auto thing : kirchoff_row){
//            std::cout << thing << ",";
//        }
//        std::cout << std::endl;

        //kirchoff_matrix[i] = kirchoff_row;
        kirchoff_matrix.push_back(kirchoff_row);

    }

    for(int i = 0; i < num_vertex; i++){
        giac::vecteur kirchoff_row;

        for(int j = 0; j < num_edge; j++){
            if(J.find(j+1) != J.end()){
                //kirchoff_row.push_back(0);
                //std::cout << "J found " << j+1 << '\n';
                continue;
            }
            //kirchoff_row[j] = sym2r(inc_matrix[j][i],exps,giac::context0);
            //kirchoff_row[j] = -inc_matrix[j][i]; //negative because the lower left incidence matrix is negative.
//            kirchoff_row.push_back(sym2r(-inc_matrix[j][i],exps,giac::context0));
            kirchoff_row.push_back(-inc_matrix[j][i]);
        }

        for(int k = num_edge; k < num_vertex + num_edge; k++){
            //kirchoff_row[k] = sym2r(0,exps,giac::context0);
            //kirchoff_row[k] = 0;
//            kirchoff_row.push_back(sym2r(0,exps,giac::context0));
            kirchoff_row.push_back(0);
        }

//        for(auto thing : kirchoff_row){
//            std::cout << thing << ",";
//        }
//        std::cout << std::endl;

        kirchoff_matrix.push_back(kirchoff_row);
    }
    
    return kirchoff_matrix;
}

std::vector<int> detect_edge_sequence(std::set<int>& I, std::set<int>& J, std::vector<std::vector<int>>& inc_matrix,
        std::set<int>& del_cont_edges){
    int num_edge= inc_matrix.size();
    int num_vertex = inc_matrix[0].size();

    //TODO We only use row information to make the edge sequence.
    std::vector<std::pair<int,std::vector<int>>> edge_sequences;
    for(int twice = 0; twice < 1; twice++){
        for(int i  = 0 ; i < num_vertex; i++){
            std::vector<int> edge_subseq;
            int connectivity = 0;
            for(int j = 0; j < num_edge ; j++){

                if(twice == 0){
                    if(I.find(j + 1) != I.end()){
                        continue;
                    }
                }

                else{
                    if(J.find(j + 1) != J.end()){
                        continue;
                    }
                }
                    

                std::cout <<inc_matrix[j][i] << ' ';
                if(inc_matrix[j][i] != 0){

                    if(del_cont_edges.find(j+1) == del_cont_edges.end()){

                        connectivity++;
                        edge_subseq.push_back(j+1);
                    }

                    else{
                        edge_subseq.push_back(-1);
                    }
                }
            }
            auto subseq_pair = std::make_pair(connectivity,edge_subseq);
            edge_sequences.push_back(subseq_pair);
            std:: cout << '\n';
        }
    }

    std::sort(edge_sequences.begin(), edge_sequences.end());
    std::vector<int> edge_sequence_good;

    while(!edge_sequences.empty()){
        auto pair = edge_sequences[0];
        auto vect = pair.second;
        if(pair.first == 0){
            edge_sequences.erase(edge_sequences.begin());
            continue;
        }
        int good_edge = -1;
        for(int edge : vect){

            if(edge == -1){
                continue;
            }

            edge_sequence_good.push_back(edge);
            good_edge = edge;
            break;
        }

        std::cout << pair.first << " " << pair.second << '\n';
        for(auto& vectpair : edge_sequences){
            for(auto it = vectpair.second.begin(); it != vectpair.second.end(); ++it){
                if(*it == good_edge){
                    vectpair.first -= 1;
                    vectpair.second.erase(it);
                    break;
               }
            }
        }

        std::sort(edge_sequences.begin(), edge_sequences.end());
    }

    std::cout << "EDGE SEEQ GOOD " << edge_sequence_good << '\n';
    return edge_sequence_good;
}

int main(int argc, char** argv)
{
    std::cout << argv[0] << '\n';

    if(argc == 1){
        printf("needs an argument\n");
        return 0;
    }

    std::vector<_graph> graphs;
    std::vector<std::vector<std::vector<int>>> incidence_matrices;
    std::vector<std::string> periods;
    std::vector<std::string> unprocessed_edges;
    std::map<_graph, std::set<int>> vertex_set;
    std::vector<std::string> all_lines;
    std::string filename = argv[1];
    std::ifstream graph_file(filename);

    // Get each line of the periods file
    if(graph_file.is_open())
    {
        std::string line;
        while(std::getline(graph_file,line)){
            all_lines.push_back(line);
            //std::cout << line << '\n';
        }
        graph_file.close();
    }

    initial_processing(all_lines,periods,unprocessed_edges); 
    populate_graphs(unprocessed_edges, graphs, vertex_set); 
    get_incidence_matrices(incidence_matrices,graphs,vertex_set);

    //Testing things.
    context ct;
    gen x(std::string("x"),&ct);
//    gen y(std::string("y"),&ct);
//    gen f(std::string("7x + 5x^2 -4"),&ct);
//    gen g(std::string("8 + 15y^2*x^2 -3"),&ct);
//    
//    gen cmd = makevecteur(f,x);
//    gen test2 = sym2r(g,makevecteur(x,y),&ct);


    //vecteur exps;
    //find_coeff_p(test2,3,makevecteur(x,y));

//    std::cout << r2e(test2,makevecteur(x),&ct) << "REDUCED MOD P \n";;

//    std::set<int> I;
//    std::set<int> J;
//    std::set<int> K;

    
    //std::cout<< kirchoff_mat << "\n";
    //auto kirchoff_mat = compute_kirchoff_matrix(exps,incidence_matrices[0],I,J,K);
    //auto pol = _det(kirchoff_mat, giac::context0) ;
    //pol = pol*pol;
    
    //std::cout<< r2e(pol,exps,giac::context0)<< '\n';
    //std::cout<< r2e(_factor(pol,giac::context0),exps,giac::context0)<< '\n';
    //pol * pol * pol * pol;
    //
    
    //CONSTRUCTION OF D^6

    //First we find the edges to construct the 6 invariant with.
    //We want 2 sets of edges that are disjoint.

    int mycount = 0;
    auto inc_matrix = incidence_matrices[0];
    int num_edge = inc_matrix.size();
    int num_vertex = inc_matrix[0].size();

    for(auto column : inc_matrix){
        for( auto element : column){
            std::cout << element << ' ';
        }

        std::cout << '\n';
    }

    //Get the 3- valent vertices and the corresponding edges.
    //TODO maybe use this info to get a better edge sequence.
    std::set<int> del_cont_edges;
    std::vector<int> edges_3_valent;
    for(int j = 0; j < inc_matrix[0].size(); j ++){
        std::vector<int> edges;
        int valency = 0; 

        for(int i= 0; i < inc_matrix.size(); i++){
            if(inc_matrix[i][j] != 0){
                valency++;
                edges.push_back(i+1);
            }
        }

        if(valency == 3){
            for(int edge : edges){
                if(mycount == 0 || mycount == 1){
                    edges_3_valent.push_back(edge);
                    del_cont_edges.insert(edge);
                    std::cout << edge << " 3valent" << '\n';
                }
            }
            mycount++;
        }
    }

    //Doing the math with the 2 3-valent vertices gets this as the 6 invariant... we label 
    //edges 1,2,6 for the first 3-valent vertex and 3,4,5 for the other.
    std::set<int> I_1 = {edges_3_valent[0],edges_3_valent[1]};
    std::set<int> J_1 = {edges_3_valent[3],edges_3_valent[4]};
    std::set<int> K_1 = {edges_3_valent[2],edges_3_valent[5]};

    //std::cout << edges_3_valent[0] << edges_3_valent[1]  << edges_3_valent[3] << edges_3_valent[4];

    std::set<int> I_2 = {edges_3_valent[0],edges_3_valent[3],edges_3_valent[5],edges_3_valent[2]};
    std::set<int> J_2 = {edges_3_valent[1],edges_3_valent[4],edges_3_valent[5],edges_3_valent[2]};
    std::set<int> K_2;

    //TODO This is a bad way to do things.
    vecteur exps_all;
    vecteur exps_topass;
    vecteur exps_test;

    std::vector<int> edge_sequence = detect_edge_sequence(I_2,J_2,inc_matrix,del_cont_edges); 
//    Kill all loops first.
//
//    //8,39 testing
//  std::vector<int> edge_sequence = {13,15,5,6,10,9,3,2,4,1};
//  std::vector<int> edge_sequence = {1,13,5,6,10,15,9,3,2,4};
//  std::vector<int> edge_sequence = {2,4,13,5,6,10,15,9,3,1};
//  std::vector<int> edge_sequence = {4,13,5,6,10,15,9,3,1,2};
//  std::vector<int> edge_sequence = {10,14,4,12,17,5,11,6,13,7,21,22,18,19,1,3};
//    std::vector<int> edge_sequence = {21,18,19,22,7,17,10,14,4,12,5,11,6,13,1,3};
//    std::vector<int> edge_sequence = {18,19,22,17,10,21,14,4,12,5,7,11,6,13,1,3};
//    std::vector<int> edge_sequence = {5,11,18,21,7,22,1,6,12,13,19,17,10,14,4,3};

    vecteur dodgson1 = compute_kirchoff_matrix(exps_topass,inc_matrix,I_1,J_1,K_1);
    auto dodgson2 = compute_kirchoff_matrix(exps_all,inc_matrix,I_2,J_2,K_2);

    std::ofstream edges3valent_file("edges3valent");
    std::ofstream edgeseq_file("edgeseq");
    std::ofstream dodgson1_file("dodgson1");
    std::ofstream dodgson2_file("dodgson2");

    for(auto elt : edges_3_valent){
        edges3valent_file << elt << '\n';
    }

    for(auto elt : edge_sequence){
        edgeseq_file << elt << '\n';
    }

    for(int i = 0; i < dodgson1.size(); i++){
        for (int j = 0; j < dodgson1.size(); j++){
            if(j != dodgson1.size()-1){
                dodgson1_file<< dodgson1[i][j] << ',';
            }

            else{
                dodgson1_file<< dodgson1[i][j];
            }
        }
        if(i != dodgson1.size()){
            dodgson1_file<< '\n';
        }
    }
    
    for(int i = 0; i < dodgson2.size(); i++){
        for (int j = 0; j < dodgson2.size(); j++){
            if(j != dodgson2.size()-1){
                dodgson2_file << dodgson2[i][j] << ',';
            }

            else{
                dodgson2_file << dodgson2[i][j];
            }
        }
        if(i != dodgson2.size()){
            dodgson2_file << '\n';
        }
    }

    edges3valent_file.close();
    edgeseq_file.close();
    dodgson1_file.close();
    dodgson2_file.close();
    

    //I_1.insert(bridge);
    //J_1.insert(bridge);
    //K_2.insert(bridge);
    //
    vecteur exp5;
    vecteur exp6;

    vecteur dodgson_1 = compute_kirchoff_matrix(exp5,inc_matrix,I_1,J_1,K_1);
    auto dodgson_2 = compute_kirchoff_matrix(exp6,inc_matrix,I_2,J_2,K_2);
    gen nth_inv = -1 * _det(dodgson_1, giac::context0) * _det(dodgson_2,giac::context0);

    std::cout << exps_all << '\n';
    //std::cout << r2e(nth_inv, exps_all, giac::context0)<< '\n';

    std::set<int> primes = {2};

    for(int p : primes){

        clock_t begin = clock();

        vecteur old_symbols = makevecteur(gen(std::string("b0_0"),&ct),gen(std::string("b0_1"),&ct));
        vecteur new_symbols;
        std::map<std::string, gen> substitute_map;
        std::map<std::string, std::vector<std::set<int>>> IJK;

        std::vector<std::set<int>> to_insert1 = {I_1,J_1,K_1};
        std::vector<std::set<int>> to_insert2 = {I_2,J_2,K_2};

        IJK.insert(std::make_pair(std::string(old_symbols[0]._IDNTptr->id_name),to_insert1));
        IJK.insert(std::make_pair(std::string(old_symbols[1]._IDNTptr->id_name),to_insert2));

        //std::cout << IJK.size() << '\n';

        gen sixinv(1);
        for(int i = 0; i < p-1; i++){
            sixinv = sixinv * old_symbols[0] * old_symbols[1];
        }

        for(int i = 0; i < edge_sequence.size(); i++){
            new_symbols.clear();
            bool last = false;

            if(i == edge_sequence.size() -1){
                last = true;
            }
            for(int j = 0; j < old_symbols.size(); j++){
                std::string key = old_symbols[j]._IDNTptr -> id_name;
                std::string id = key.substr(key.find('_') + 1);
                int term = std::stoi(id);

                std::set<int> i_1 = IJK[key][0];
                std::set<int> j_1 = IJK[key][1];
                std::set<int> k_1 = IJK[key][2];
                std::set<int> i_2 = IJK[key][0];
                std::set<int> j_2 = IJK[key][1];
                std::set<int> k_2 = IJK[key][2];

                int edge = edge_sequence[i];
                i_2.insert(edge);
                j_2.insert(edge);
                k_1.insert(edge);

                bool term1_zero = true;
                bool term2_zero = true;

                vecteur exps1;
                vecteur exps2;

                auto dodgson1 = _det(compute_kirchoff_matrix(exps1,inc_matrix,i_1,j_1,k_1),&ct);

                std::string var1 = "b" + std::to_string(i+1) + "_" + std::to_string(term * 2+1);
                gen a_var = gen(var1,&ct);
                if(!is_zero(dodgson1,&ct)){
                    term1_zero = false;
                    new_symbols.push_back(a_var);

                    IJK.insert(std::make_pair(std::string(a_var._IDNTptr-> id_name), std::vector<std::set<int>>{i_1,j_1,k_1}));
                }


                auto dodgson2 = _det(compute_kirchoff_matrix(exps2,inc_matrix,i_2,j_2,k_2),&ct);
                std::string var2 = "b" + std::to_string(i+1) + "_" + std::to_string(term * 2);
                gen b_var = gen(var2,&ct);
                if(!is_zero(dodgson2,&ct)){
                    term2_zero = false;
                    new_symbols.push_back(b_var);

                    IJK.insert(std::make_pair(std::string(b_var._IDNTptr->id_name), std::vector<std::set<int>>{i_2,j_2,k_2}));
                }
                
                if(!term1_zero){
                    std::cerr << var1 << ',' << dodgson1 << '\n';
                }

                if(!term2_zero){
                    std::cerr << var2 << ',' << dodgson2 << '\n';
                }

                gen subvalue;
                
                if(!last){
                    if(term2_zero){
                        subvalue = a_var; 
                    }

                    else if (term1_zero){
                        subvalue = x * b_var;
                    }
                    
                    else{
                       subvalue = x*b_var + a_var; 
                    }
                }

                else{
                    if(term2_zero){
                        subvalue = dodgson1;
                    }

                    else{
                        subvalue = dodgson2;
                    }
                }


                substitute_map.insert(std::make_pair(key, subvalue));

                if(term1_zero && term2_zero){
                    std::cerr << "both are zero\n";
                    exit(0);
                }
            }

            vecteur tosub;
            for(int k = 0; k < old_symbols.size(); k++){
                std::string identifier = std::string(old_symbols[k]._IDNTptr -> id_name);
                tosub.push_back(substitute_map[identifier]);

            }

            gen subcmd = makesequence(sixinv,old_symbols,tosub);
            gen poly_in_x = _subst(subcmd,&ct);

            gen coeffcmd = makesequence(poly_in_x,x,p-1);

            if(!last){
                poly_in_x = _coeff(coeffcmd,&ct);
            }

            
            sixinv = poly_in_x;
            std::cerr << "<" << sixinv << "ITERATION  " << i << '\n';
            old_symbols = new_symbols;
        }

        clock_t end = clock();

        std::cout << (end-begin)/(double)CLOCKS_PER_SEC << " seconds taken. coeff for p = " << p << ": " << sixinv << ", c_2 is "
            << sixinv % p << '\n';
    }
}


