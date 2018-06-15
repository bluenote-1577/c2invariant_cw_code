#include <string>
#include <cstring>
#include <sstream>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <unordered_map>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>
#include <ostream>
#include <unordered_set>
#include <set>
#include <mpi.h>

//typedef std::unordered_map<std::string,int> variables;
const unsigned int edgestop = 3;
int p = 3;
std::unordered_map<int,int> zero_dodgson_map;
std::unordered_map<int,std::vector<std::pair<std::vector<int>,bool>>> dodgson_map;
std::map<int,int> tocheck;
std::unordered_set<int> allvars;
int numvars = 0;

struct monomial{

    std::map<int,int> vars;

    bool operator==(const monomial &other) const{
        return (vars == other.vars);
    }

    bool operator < (const monomial &other) const{
        
        if(vars.size() != other.vars.size()){
            return vars.size() < other.vars.size();
        }

        int id1 = 1;
        int id2 = 1;
		for(auto it = vars.begin(); it != vars.end(); ++it){
            id1 += 7 *id1 * it->first * it->second;
        }

		for(auto it = other.vars.begin(); it != other.vars.end(); ++it){
            id2 += 7 * id2 * it->first * it->second;
        }

        return id1 < id2;
    }

    monomial(std::map<int,int> a_vars) : vars(a_vars){}
    monomial(const monomial &m){
        vars = m.vars;
    } 
    monomial(std::vector<std::pair<int,int>> a_vars){
        for(unsigned int i = 0; i < a_vars.size(); i++){
           vars.insert(a_vars[i]); 
        }
    }
    monomial(int i, int j){
        vars[i] = j;
    }
    monomial(){};
};

typedef std::unordered_map<monomial,long long int> polynomial;
typedef std::map<monomial,long long int> map_polynomial;

struct vect_mono_key{

    bool operator==(const vect_mono_key &other) const{
        return (vars == other.vars);
    }

    std::vector<int> vars;

    vect_mono_key(int* mon): vars(mon, mon+ numvars + 1){}
    
    vect_mono_key(){};
    vect_mono_key(int null){
        vars = std::vector<int>(null);
    };

};

struct vect_mono{
    

    bool operator==(const vect_mono &other) const{
        return (vars == other.vars);
    }

    std::vector<std::pair<int,int>> vars;

    vect_mono(const monomial& mon)
    {
        vars.reserve(mon.vars.size() + 1);
		for(auto cfpair = mon.vars.begin(); cfpair != mon.vars.end(); ++cfpair){
            this->vars.push_back(*cfpair);
        }
    }

    vect_mono(int* mon){
        for(int i = 0; i < numvars + 1; i++){
            if(mon[i] != 0){
                vars.push_back(std::make_pair(i,mon[i]));
            }
        }
    }
    
    vect_mono(){};

};

typedef std::unordered_map<vect_mono_key,long long int> polynomial_vect;
//typedef ska::flat_hash_map<vect_mono_key,long long int> polynomial_vect;
//typedef google::dense_hash_map<vect_mono_key,long long int> polynomial_vect;
//typedef tsl::hopscotch_map<vect_mono_key,long long int> polynomial_vect;

struct monomial_x{
    vect_mono mono;
    long int cf;
    int exponent;
    monomial_x(std::vector<std::pair<int,int>> a_mono, int a_cf, int a_exponent) : cf(a_cf), exponent(a_exponent){
        mono.vars = a_mono; 
    }
};

typedef std::vector<std::pair<vect_mono,int>> iter_polyn; 

void code_input(std::string& line){
    
    std::string id = line.substr(0,line.find(','));
    int rest = std::stoi(line.substr(line.find('_') + 1, line.find(',')));
    int iter = std::stoi(line.substr(1, line.find('_')));
    std::string dodgson = line.substr(line.find(',')+1);
    int numdigits_rest = 1;
    int tomult = 10;

    int temp = rest;
    temp /= 10;
    while(temp!= 0){
        temp /= 10;
        numdigits_rest++;
        tomult *= 10;
    }

    int to_insert;

    if(iter > 9){
        to_insert = 9 * tomult * 10 * 10 + tomult * iter + rest;
        to_insert *= -1;
    }

    else{
        to_insert = 9 * tomult * 10 + tomult * iter + rest;
    }
    
    std::stringstream ss;
    ss << dodgson;
    int numb = 0;

    if(ss >> numb){
        zero_dodgson_map.insert(std::make_pair(to_insert,numb));
        return;
    }


    std::vector<std::string> minus_tokens;
    std::vector<std::string> mono_tokens;

    size_t pos = 0;
    std::string token;

    while((pos = dodgson.find('+')) != std::string::npos){
        token = dodgson.substr(0,pos);
        minus_tokens.push_back(token);
        dodgson.erase(0,pos+1);
    }

    minus_tokens.push_back(dodgson);

//    for(auto strn : minus_tokens){
	for(auto strn = minus_tokens.begin(); strn!= minus_tokens.end(); ++strn){
        while((pos = strn->find('-',1)) != std::string::npos){
            token = strn->substr(0, pos);
            mono_tokens.push_back(token);
            strn->erase(0,pos);
        }

        mono_tokens.push_back(*strn);
    }
    
    std::vector<std::pair<std::vector<int>,bool>> dodgson_poly;
//    for(auto s : mono_tokens){
	for(auto s = mono_tokens.begin(); s != mono_tokens.end(); ++s){
        std::vector<int> mono;
        bool is_neg = false;
        bool is_first = true;
        if((*s)[0] == '-'){
                is_neg = true;
        }

        while((pos = s->find('*')) != std::string::npos){

           int val; 
           if(is_neg && is_first){
                token = s->substr(1, pos-1);
                val = std::stoi(token.substr(2));
                is_first = false;
           }

            else{
               token = s->substr(0,pos);
                val = std::stoi(token.substr(2));
            }


            int temp = val / 10;
            int tomult = 10;

            while(temp != 0){
                temp /= 10;
                tomult *= 10;
            }

//            int id = tomult* 8 + val;
            int id = val;
            mono.push_back(id);
            s->erase(0,pos+1);
        }

        int val;
        if((*s)[0] == '-'){
            is_neg = true;
            val = std::stoi(s->substr(3));
        }
        else{
            val = std::stoi(s->substr(2));
        }
        int temp = val / 10;
        int tomult = 10;

        while(temp != 0){
            temp /= 10;
            tomult *= 10;
        }

//        int id = tomult* 8 + val;
        int id = val;

//        std::cout << tomult * 8 + val << '\n';;
//        std::cout << val << '\n';
        if(val > numvars){
            numvars = val;
        }

        mono.push_back(id);
        dodgson_poly.push_back(std::make_pair(mono,is_neg));
    }
    
    zero_dodgson_map[to_insert] = 777;
    dodgson_map[to_insert] = dodgson_poly;
}




namespace std{

    template <>
    struct hash<monomial> {
        std::size_t operator()(const monomial& k) const
        {
            std::size_t toreturn = 17;

            using std::hash;
            using std::string;

//            for(auto& key : k.vars){
			for(auto key = k.vars.begin(); key != k.vars.end(); ++key){
                toreturn = toreturn * 31 + ((hash<int>()(key->first) << 1 ) >> 1) * (key->second * 17);
            }

            return toreturn;
        }
    };

    template <>
    struct hash<vect_mono_key> {
        std::size_t operator()(const vect_mono_key& k) const
        {
            std::size_t toreturn = 37;

            using std::hash;

//             return seed;
            int mycount = 1;
			for(auto key = k.vars.begin(); key != k.vars.end(); ++key){
//            for(auto key : k.vars){
//                toreturn = toreturn * 31 + ((hash<int>()(key) << 1 ) >> 1);
//                toreturn = toreturn * 31 + ((key << 1 ) >> 1);
                toreturn ^= *key + 0x9e3779b9 + (toreturn << 6) + (toreturn >> 2);
                mycount++;
            }

            return toreturn;
        }
    };

}


std::ostream& operator<<(std::ostream& os, const monomial& mo){

//    for(auto vals : mo.vars){
    for(auto vals = mo.vars.begin(); vals != mo.vars.end(); ++vals){
        if(vals->second > 1){
            os << vals->first << '^' << vals->second << '*';
        }

        else{
            os << vals->first << '*';
        }
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const monomial_x& m){
        os << m.cf << "*";

        for (unsigned int i = 0; i < m.mono.vars.size(); i++){
            os << m.mono.vars[i].first << "^" <<  m.mono.vars[i].second << ' ';
        }

        return os;
}

std::ostream& operator<<(std::ostream& os, const vect_mono& m){

        for (unsigned int i = 0; i < m.vars.size(); i++){
            os << m.vars[i].first << "^" <<  m.vars[i].second << '*';
        }

        return os;
}

std::ostream& operator<<(std::ostream& os, const vect_mono_key& m){

        for (unsigned int i = 0; i < m.vars.size(); i++){
            if(m.vars[i] != 0){
                os << i << "^" <<  m.vars[i]<< '*';
            }
        }

        return os;
}

bool sort_function ( const std::vector<monomial_x>& vec1, const std::vector<monomial_x>& vec2){
    return (vec1.size() < vec2.size());
}

void
decode (int var, int num_digits, int* toret){
    int ten_it = 1;
    int temp = var;
    std::vector<int> digits;

    if(var < 0){
        ten_it = 2;
        temp = -1 * var;
    }

    for(int i = 1; i < num_digits; i++){
        digits.insert(digits.begin(),temp % 10);
        temp /= 10;
    }

    int iter;
    if(ten_it == 2){
        iter = digits[0] * 10 + digits[1];
    }
    else{
        iter = digits[0];
    }
    iter++;

    int rest1 = 0;
    int rest2 = 0;

    for(unsigned int i = ten_it; i < digits.size(); i ++){
        rest2 += digits[i];
        rest1 += digits[i]; 

        if(i != digits.size() - 1){
        rest1 *= 10;
        rest2 *= 10;
        }

    }

    rest1 *= 2;
    rest2 *= 2;
    rest2++;

    int numdigits1 = 1;
    int numdigits2 = 1;

    int temp1 = rest1 / 10;
    int temp2 = rest2 / 10;

    while(temp1 != 0){
        temp1 /= 10;
        numdigits1++;
    }

    while(temp2 != 0){
        temp2 /= 10;
        numdigits2++;
    }

    int tomult = 1;
    int tomult2 = 1;

    for(int i = 0; i < numdigits1; i++){
        tomult *= 10;
    }

    for(int i = 0; i < numdigits2; i++){
        tomult2 *= 10;
    }

    int a_term1 = 1;
    int a_term2 = 1;

    if(iter < 10){
        a_term1 = 9 * tomult * 10 + iter * tomult + rest1;
        a_term2 = 9 * tomult2 * 10 + iter * tomult2 + rest2;
    }

    else{
        a_term1 = -1 * (9 * tomult * 10 * 10 + iter * tomult + rest1);
        a_term2 = -1 * (9 * tomult2 * 10 * 10+ iter * tomult2 + rest2);
    }

    toret[0] = a_term1;
    toret[1] = a_term2;

}

int 
tabled_lookup_fin (int var, int exponent){
    int num_digits = 1;
    int tempo = var / 10;

    while(tempo !=0){
        tempo /= 10;
        num_digits++;
    }

    int toret[2] = {0,0};
    decode(var,num_digits,toret);    

    int a_term1 = toret[0];
    int a_term2 = toret[1];
    
    if(zero_dodgson_map[a_term1] == 0){

        if(exponent % 2 == 0){
            return zero_dodgson_map[a_term2] * zero_dodgson_map[a_term2];
        }
        else{
            return zero_dodgson_map[a_term2];
        }
    }
    
    else if(zero_dodgson_map[a_term2] == 0){
        if(exponent % 2 == 0){
            return zero_dodgson_map[a_term1] * zero_dodgson_map[a_term1];
        }
        else{
            return zero_dodgson_map[a_term1];
        }
    }

    else{std::cout <<  var << '\n' << a_term1 << a_term2 << '\n'; exit(0);};

}


//tabled lookup just to p=11;
std::vector<monomial_x>
tabled_lookup (int var, int exponent, bool sub_dodgs = false){ 
    int num_digits = 0;
    int tempo = var;

    while(tempo !=0){
        tempo /= 10;
        num_digits++;
    }

    bool a_zero = false;
    bool b_zero = false;

    int toret[2] = {0,0};
    decode(var,num_digits,toret);    

    int a_term1 = toret[0];
    int a_term2 = toret[1];
    
    if(zero_dodgson_map.find(a_term1) == zero_dodgson_map.end()){
        a_zero = true;
    }
    
    if(zero_dodgson_map.find(a_term2) == zero_dodgson_map.end()){
        b_zero = true;

    }

    if(sub_dodgs){
        std::vector<monomial_x> toret;      
        const auto intcheck = zero_dodgson_map.find(var);

        if(intcheck -> second == 777){

//            for(const auto& boolvect : (dodgson_map.find(var)->second)){
            for(auto boolvect = (dodgson_map.find(var)->second).begin(); boolvect != (dodgson_map.find(var)->second).end(); 
					++boolvect){
                std::vector<std::pair<int,int>>topush;
                int cf = 1;

                if(boolvect->second){
                    cf = -1;
                }

//                for(const auto & id : boolvect->first){
                for(auto id = boolvect->first.begin(); id != boolvect->first.end(); ++id){
                    topush.push_back(std::make_pair(*id,1));
                }

                toret.push_back(monomial_x(topush,cf,exponent));
            }
        }
    
        else{
            toret.push_back(monomial_x(std::vector<std::pair<int,int>>(),intcheck->second,exponent));
        }

        return toret;
    }
    
    if(a_zero){
//        if(intcheck2 != 777){
//            int cf = intcheck2;
//            if(exponent % 2 == 0){
//                cf = 1;
//            }
//           return std::vector<monomial_x>({monomial_x(monomial(0,0),cf,0)});
//        }
//        return std::vector<monomial_x>({monomial_x(monomial(a_term2,exponent),1,0)});
        return std::vector<monomial_x>({monomial_x({std::make_pair(a_term2,exponent)},1,0)});
    }

    if(b_zero){
//        if(intcheck1 != 777){
//            int cf = intcheck1;
//            if(exponent % 2 == 0){
//                cf = 1;
//            }
//           return std::vector<monomial_x>({monomial_x(monomial(0,0),cf,exponent)});
//        }
       return std::vector<monomial_x>({monomial_x({std::make_pair(a_term1,exponent)},1,exponent)});
    }

    if(a_zero && b_zero){
        std::cout << "both zero" << '\n';
        exit(0);
    }

    switch(exponent){
        case 1: {
                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,1)},1,0),
                   monomial_x({std::make_pair(a_term1,1)},1,1)
                           });
                }

        case 2: {

                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,2)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,1)},2,1),
                   monomial_x({std::make_pair(a_term1,2)},1,2)
                           });
                }

        case 3: {

                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,3)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,2)},3,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,1)},3,2),
                   monomial_x({std::make_pair(a_term1,3)},1,3)
                   });
                }
        case 4: {
                   
                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,4)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,3)},4,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,2)},6,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,1)},4,3),
                   monomial_x({std::make_pair(a_term1,4)},1,4)
                           });
                }

        case 5: {
                   
                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,5)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,4)},5,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,3)},10,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,2)},10,3),
                   monomial_x({std::make_pair(a_term1,4),std::make_pair(a_term2,1)},5,4),
                   monomial_x({std::make_pair(a_term1,5)},1,5)
                           });
                }
       case 6: {
                   
                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,6)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,5)},6,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,4)},15,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,3)},20,3),
                   monomial_x({std::make_pair(a_term1,4),std::make_pair(a_term2,2)},15,4),
                   monomial_x({std::make_pair(a_term1,5),std::make_pair(a_term2,1)},6,5),
                   monomial_x({std::make_pair(a_term1,6)},1,6)
                           });
                }

      case 7: {
                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,7)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,6)},7,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,5)},21,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,4)},35,3),
                   monomial_x({std::make_pair(a_term1,4),std::make_pair(a_term2,3)},35,4),
                   monomial_x({std::make_pair(a_term1,5),std::make_pair(a_term2,2)},21,5),
                   monomial_x({std::make_pair(a_term1,6),std::make_pair(a_term2,1)},7,6),
                   monomial_x({std::make_pair(a_term1,7)},1,7)
                           });

              }

      case 8: {
                   return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,8)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,7)},8,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,6)},28,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,5)},56,3),
                   monomial_x({std::make_pair(a_term1,4),std::make_pair(a_term2,4)},70,4),
                   monomial_x({std::make_pair(a_term1,5),std::make_pair(a_term2,3)},56,5),
                   monomial_x({std::make_pair(a_term1,6),std::make_pair(a_term2,2)},28,6),
                   monomial_x({std::make_pair(a_term1,7),std::make_pair(a_term2,1)},8,7),
                   monomial_x({std::make_pair(a_term1,8)},1,8)
                           });
              }

      case 9: {
                  return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,9)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,8)},9,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,7)},36,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,6)},84,3),
                   monomial_x({std::make_pair(a_term1,4),std::make_pair(a_term2,5)},126,4),
                   monomial_x({std::make_pair(a_term1,5),std::make_pair(a_term2,4)},126,5),
                   monomial_x({std::make_pair(a_term1,6),std::make_pair(a_term2,3)},84,6),
                   monomial_x({std::make_pair(a_term1,7),std::make_pair(a_term2,2)},36,7),
                   monomial_x({std::make_pair(a_term1,8),std::make_pair(a_term2,1)},9,8),
                   monomial_x({std::make_pair(a_term1,9)},1,9)
                           });

              }
      case 10:{
                    return std::vector<monomial_x>({
                   monomial_x({std::make_pair(a_term2,10)},1,0),
                   monomial_x({std::make_pair(a_term1,1),std::make_pair(a_term2,9)},10,1),
                   monomial_x({std::make_pair(a_term1,2),std::make_pair(a_term2,8)},45,2),
                   monomial_x({std::make_pair(a_term1,3),std::make_pair(a_term2,7)},120,3),
                   monomial_x({std::make_pair(a_term1,4),std::make_pair(a_term2,6)},210,4),
                   monomial_x({std::make_pair(a_term1,5),std::make_pair(a_term2,5)},252,5),
                   monomial_x({std::make_pair(a_term1,6),std::make_pair(a_term2,4)},210,6),
                   monomial_x({std::make_pair(a_term1,7),std::make_pair(a_term2,3)},120,7),
                   monomial_x({std::make_pair(a_term1,8),std::make_pair(a_term2,2)},45,8),
                   monomial_x({std::make_pair(a_term1,9),std::make_pair(a_term2,1)},10,9),
                   monomial_x({std::make_pair(a_term1,10)},1,10)
                           });

              }


        default:{
                    //TODO
                    std::cout << "shouldn't reach here, exponent > 5\n ";
                    exit(0);
                }

    }
};



void
obtain_multiplicands ( 
        std::vector<std::vector<monomial_x*>>& multiplicands, 
        int current_weight,
        std::vector<std::vector<monomial_x>>& to_mult,
        unsigned int iteration, std::vector<monomial_x*>& own_mults){

    
    if(iteration == to_mult.size()){
        multiplicands.push_back(own_mults);


        return;
    }
    
    int max_weight_add = 0;

    for(unsigned int k = iteration + 1; k < to_mult.size(); k++){
        if(to_mult[k].size() != 1){
            max_weight_add += to_mult[k].size() - 1; 
        }

        else{
            max_weight_add += to_mult[k][0].exponent;
        }
    }

    for(unsigned int i = 0; i < to_mult[iteration].size(); i++){
//        monomial_x mon_x = to_mult[iteration][i];

        if(to_mult[iteration][i].exponent + current_weight + max_weight_add < p-1){
            continue;
        }

        else if (to_mult[iteration][i].exponent + current_weight > p-1){
            continue;
        }


        std::vector<monomial_x*> copy_of_own_mults = own_mults;
        copy_of_own_mults.push_back(&(to_mult[iteration][i]));
        obtain_multiplicands(multiplicands, current_weight + to_mult[iteration][i].exponent, to_mult, iteration + 1, copy_of_own_mults);

    }
    
};

int main (int argc, char*argv[]){

//    omp_set_num_threads(4);

	int myid,numprocs;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    std::ifstream dodgson_file("gen_out.txt");
	int max = 0;
	int temp = 0;
    // Get each line of the periods file
    if(dodgson_file.is_open())
    {
        std::string line;
        while(std::getline(dodgson_file,line)){
            if(line[0] != '<'){
			 code_input(line);
            }
        }
        dodgson_file.close();
    }
	
    ////////////////////////////

	std::set<int> primes = {11};
	for(auto it = primes.begin(); it != primes.end(); ++it){
		p = *it;
		vect_mono start;
		start.vars = {std::make_pair(900,p-1),std::make_pair(901,p-1)};
		iter_polyn actual_poly = {std::make_pair(start,1)};
		const int con_numvars = numvars;
		long c2cf = 0;

		for(unsigned int i = 0; i < edgestop; i++){
			polynomial toadd;
//			std::cout << "size " << actual_poly.size() << '\n';

			double begin = omp_get_wtime();
			
			for(unsigned int k = 0; k < actual_poly.size(); k++){

				auto* mono_it = &actual_poly[k];

//				if(k%100 == 0){
//					std::cout << k << '\n';
//				}
				std::vector<std::vector<monomial_x>> to_mult;

				for(auto elt = mono_it->first.vars.begin(); elt !=  mono_it->first.vars.end(); ++elt){
					to_mult.push_back(tabled_lookup(elt->first,elt->second,false));
				}


				std::sort (to_mult.begin(), to_mult.end(), sort_function);
				std::vector<std::vector<monomial_x*>> multiplicands;
				std::vector<monomial_x*> own_mults;
				obtain_multiplicands(multiplicands, 0, to_mult, 0,own_mults);

				unsigned int vecsize = 0;
				if(!multiplicands.empty()){
					vecsize =  multiplicands[0].size();
				}
				for(auto vec = multiplicands.begin(); vec != multiplicands.end(); ++vec){
				   monomial toins; 
				   int cf = mono_it->second;
				   if(cf != 0){
						//TODO Multiply in chunks here.
					   for(unsigned int i = 0; i < vecsize; i++){
						   for(auto keypair = (*vec)[i]->mono.vars.begin(); keypair != (*vec)[i]->mono.vars.end() ; ++keypair){
							   toins.vars[keypair->first] += keypair->second;
						   }

						   cf = cf * (*vec)[i]->cf;
					   }
					   toadd[toins] += cf;
	//                   toadd[toins] %= p;
				   }
				}
			}
			//Look into std::move ... TODO
			actual_poly.clear();
			for(auto monomial = toadd.begin(); monomial != toadd.end(); ++monomial){
				if(monomial->second != 0){
					actual_poly.push_back(std::make_pair(vect_mono(monomial->first),monomial->second));
				}
			}

			double end = omp_get_wtime();
		}
		

		double begin = omp_get_wtime();
		if(myid != 0){
			unsigned int chunksize = actual_poly.size() / (numprocs-1);
			unsigned int range_a = (myid-1) * chunksize;
			unsigned int range_b = range_a + chunksize;
			if(myid == numprocs - 1){
				range_b = range_a + chunksize + actual_poly.size() % (numprocs-1);	
			}

			for(unsigned int k = range_a; k < range_b ; k++){
				auto* mono_it = &actual_poly[k];
//				if(k%100 == 0){
//					std::cout << k << '\n';
//				}
				std::vector<std::vector<monomial_x>> to_mult;
				for(auto elt = mono_it->first.vars.begin(); elt !=  mono_it->first.vars.end(); ++elt){
					to_mult.push_back(tabled_lookup(elt->first,elt->second,true));
				}
			   polynomial_vect* intermediate_poly = new polynomial_vect(); 
			   (*intermediate_poly)[vect_mono_key()] = mono_it->second;
			   
			   ////////////////////////////////////////////////////////////
			   //CRITICAL BLOCK, MOST OF THE TIME IS SPENT HERE
			   ////////////////////////////////////////////////////////////
			   for(unsigned int i = 0; i < to_mult.size(); i++){
					for(int k = 0; k < to_mult[i][0].exponent; k++){
					   polynomial_vect* temp_poly = new polynomial_vect();
					   std::vector<vect_mono_key> to_erase;
						for(unsigned int j = 0; j < to_mult[i].size(); j++){
							for(auto pr = intermediate_poly->begin(),pre = intermediate_poly->end(); pr != pre ; ++pr){
								int cf = pr->second;
								int mon_ins[con_numvars+1];
								memset(mon_ins,0, sizeof(mon_ins));

								 int m = 0;
								 for(auto pr1 = pr->first.vars.begin(),pr1e = pr->first.vars.end(); pr1 != pr1e ; ++pr1){
									 mon_ins[m] += *pr1;
									 m++;
								 }

								bool power_too_high = false;

								auto it = to_mult[i][j].mono.vars.begin();
								auto end = to_mult[i][j].mono.vars.end();

								while(!power_too_high && it != end){
								   mon_ins[it -> first] += it -> second ; 
								   if(mon_ins[it -> first] > p-1){
									   power_too_high = true;
								   }
								   ++it;
								}
									
								if(!power_too_high){
									cf *= to_mult[i][j].cf;
									auto key = vect_mono_key(mon_ins);
									auto check = temp_poly->find(key);

									if(check != temp_poly->end()){
										check->second += cf;
										check->second %= p;
									}

									
									else{
										temp_poly->insert(std::make_pair(key,cf));
									}
								}
							}
						}
						delete intermediate_poly;

						for(auto pr = temp_poly->begin(),pre = temp_poly->end(); pr != pre ; ++pr){
							if(pr->second == 0){
								to_erase.push_back(pr->first);
							}
						}

						for(auto key = to_erase.begin(); key != to_erase.end(); ++key){
							temp_poly -> erase((*key));
						}

						intermediate_poly = temp_poly;
					}
			   }

				if(intermediate_poly->size() ==  1){ 
					{
						MPI_Send (&(intermediate_poly->begin()->second), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	//				c2cf = c2cf + intermediate_poly->begin()->second;
					}
				}
			}
			int neg = -1;

			MPI_Send(&neg,1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		}

		else{
			int partial_c2 = 0;
			int proc_finished_counter = 0;
			while(true){
				MPI_Recv(&partial_c2,1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				if(status.MPI_TAG == 0){
					c2cf += partial_c2;
				}

				else{
					proc_finished_counter++;
				}

				if(proc_finished_counter == numprocs-1){
					break;
				}
			}

			 std::cout << "broken \n";

			std::cout << p << "PRIME" ;
			std::cout << c2cf << " c2 = " << c2cf % p << '\n';;
			double end = omp_get_wtime();
			std::cout <<(end-begin) << '\n';
		}
	}
	MPI_Finalize();
}

