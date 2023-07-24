#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <random>
#include <ctime>
#include <iomanip>
#include "test_function.hpp"

using namespace std;

typedef vector<double> vec;
default_random_engine generator(time(NULL));
uniform_real_distribution<float> unif(0.0, 1.0);
normal_distribution<double> normal_distribute(0, 1);


// data
vector<vec> position;
vec fitness; 

vec best_sol;
vec worst_sol;
double best_fitness; 
double worst_fitness;
int dimension;
int population;
int object_function;
double maxx;
double minn;

// Other function
void write_record(string, vec);
double calculate_avg(vec);
double check_bound(double);
void print_vec_group(vector<vec>);
void print_vec(vec);
vector<int> arg_sort(vec);


// SSA function
void ranking();
void initialization(int, int, int);
double evaluation(vec);
vector<vec> update_position(double, double, double, int);
void determination(vec);

double levy_flight(double, double);
vector<vec> new_update_position(double, double, double, int);
vector<vec> new_update_position_2(double, double, double, int);
vector<vec> new_update_position_3(double, double, double, int);

/* -------------------------------------------- */


void write_record(string name, vec all){

    string filename = "D.txt";
    string bottom_line = "_";
    filename = name + bottom_line + to_string(dimension) + filename;

    ofstream myfile(filename);
    cout << scientific << setprecision(8);
    for(int i = 0; i < all.size(); i++){

        myfile << all[i] << endl;
    }
    
    myfile.close();

}


double calculate_avg(vec v){

    double sum = 0.0;

    for(int i = 0; i < v.size(); i++){
        sum += v[i];
    }

    return (sum / v.size());
} 


double check_bound(double num){

    if(num < minn){
        num = minn;
    }
    else if(num > maxx){
        num = maxx;
    }

    return num;
}


void print_vec_group(vector<vec> v){

    for (int i = 0; i < v.size(); i++){
        for (int k = 0; k  <v[i].size(); k++){
            cout << v[i][k] << " ";
        }
        cout << endl;
    }
}


void print_vec(vec v){

    for(auto i:v){

        cout << i << endl;
    }
}


vector<int> arg_sort(vec v){

    int pop_size = v.size();
    vector<int> index(pop_size, 0);
    iota(index.begin(), index.end(), 0);

    sort(index.begin(), index.end(), [&v](double a, double b){return (v[a]<v[b]); });

    return index;
}

/* -------------------------------------------- */


void ranking(){

    vec temp_fitness(population, 0.0);
    vector<int> rank_index;
    vector<vec> temp_position(population, vec(dimension, 0.0));
    int index;

    rank_index = arg_sort(fitness);

    for(int i = 0; i < population; i++){
        temp_fitness[i] = fitness[rank_index[i]];
        for(int k = 0; k < dimension; k++){
            temp_position[i][k] = position[rank_index[i]][k];
        }
    }

    fitness = temp_fitness;
    position = temp_position;

    best_sol = position[0];
    best_fitness = fitness[0];
    worst_sol = position[population-1];
    worst_fitness = fitness[population-1];

}


void initialization(int dim, int pop, int test_func){

    double fit;
    vec individual_pos;

    dimension = dim;
    population = pop;
    object_function = test_func;

    position.clear();
    fitness.clear();
    
    switch (test_func){

        case 1: //ackley
            minn = -32.768;
            maxx = 32.768;
            break;

        case 2: //griewank
            minn = -600.0;
            maxx = 600.0;
            break;
        
        case 3: //katsuuras
            minn = -100.0;
            maxx = 100.0;
            break;
        
        case 4: //michalewicz
            minn = 0.0;
            maxx = M_PI; 
            break;

        case 5: //rosenbrock
            minn = -10.0;
            maxx = 10.0;
            break;
        
        case 6: //schwefel
            minn = -500.0;
            maxx = 500.0;
            break;
        
        case 7: //bent_cigar
            minn = -100.0;
            maxx = 100.0;
            break;
        
        case 8: //zakharov
            minn = -10.0;
            maxx = 10.0;
            break;

        case 9: //happy_cat
            minn = -20.0;
            maxx = 20.0;
            break;
        
        case 10: //rastrigin
            minn = -5.12;
            maxx = 5.12;
            break;
        
        case 11: //hgbat
            minn = -15.0;
            maxx = 15.0;
            break;


        default:
            break;
    }


    // random initial position
    for(int i = 0; i < pop; i++){

        individual_pos.clear();

        for(int k = 0; k < dim; k++){
            double rand_num = (maxx - minn) * (rand() / (RAND_MAX + 1.0)) + minn;
            individual_pos.push_back(rand_num);
        }

        position.push_back(individual_pos);
    }

    // calculate initial fitness
    for(int i = 0; i < pop; i++){
        fit = evaluation(position[i]);
        fitness.push_back(fit);
    }

}


double evaluation(vec sol){

    double fit;

    switch (object_function){

        case 1: // ackley
            fit = ackley(dimension, sol);
            break;

        case 2: //griewank
            fit = griewank(dimension, sol);
            break;
        
        case 3: //katsuuras
            fit = katsuuras(dimension, sol);
            break;
        
        case 4: //michalewicz
            fit = michalewicz(dimension, sol);
            break;

        case 5: //rosenbrock
            fit = rosenbrock(dimension, sol);
            break;
        
        case 6: //schwefel
            fit = schwefel(dimension, sol);
            break;
        
        case 7: //bent_cigar
            fit = bent_cigar(dimension, sol);
            break;
        
        case 8: //zakharov
            fit = zakharov(dimension, sol);
            break;

        case 9: //happy_cat
            fit = happy_cat(dimension, sol);
            break;
        
        case 10: //rastrigin
            fit = rastrigin(dimension, sol);
            break;
        
        case 11: //hgbat
            fit = hgbat(dimension, sol);
            break;

        default:
            break;
    }

    return fit;

}


vector<vec> update_position(double pd, double sd, double st, int iter_max){

    vector<vec> new_position(population, vec(dimension, 0.0));

    // ----- for update producer
    int producer_num = population * pd;
    double r1, r2;  // r1: alpha -> random number (0,1]
    double new_pos;
    float q;  // random number from normal distribution
    
    double new_best_f = 0.0;
    vec new_best_sol(dimension, 0.0);

    // ----- for update scounger
    vec A(dimension, 0.0);
    double a_num, move;

    // ----- for update danger sparrow
    int danger_num = population * sd;
    vector<int> random_list(population, 0.0);
    int di; // index after random


    /********************************** update producer ******/
    
    r2 = ((double) rand() / (RAND_MAX));

    if(r2 < st){
        for(int i = 0; i < producer_num; i++){
            r1 = ((double) rand() / (RAND_MAX));
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k]*exp(-i/(r1*iter_max));
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }
    else{
        for(int i = 0; i < producer_num; i++){
            q = unif(generator);
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k] + q;
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }


    /********************************** update scrounger ******/

    for(int i = 0; i < (population-producer_num); i++){
        int scrounger = i + producer_num;

        // create A: random 1 or -1 -> matrix(1*d)
        for(int d = 0; d < dimension; d++){
            a_num = ((double) rand() / (RAND_MAX));
            a_num = floor(a_num*2)*2-1;
            A[d] = a_num;
        }

        // update starving spparow (ranking lower): need to find other place for food
        if(scrounger > (population/2)){
            q = unif(generator);

            for(int d = 0; d < dimension; d++){
                new_pos = q * exp((worst_sol[d]-position[scrounger][d])/pow(scrounger, 2));
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
        // update other scrounger
        else{
            move = 0.0;
            for(int d = 0; d < dimension; d++){
                move += abs(position[scrounger][d] - new_best_sol[d]) * (1.0/dimension) * A[d];
            }

            for(int d = 0; d < dimension; d++){
                new_pos = best_sol[d] + move;
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
    }


    /********************************** update sparrow which are aware of dager ******/
    
    iota(random_list.begin(), random_list.end(), 0);
    random_shuffle(random_list.begin(), random_list.end());

    for(int i = 0; i < danger_num; i++){

        di = random_list[i];

        if(fitness[di] > best_fitness){; 
            
            for(int d = 0; d < dimension; d++){
                q = unif(generator);
                new_pos = best_sol[d] + q*abs((position[di][d]-best_sol[d]));
                new_position[di][d] = check_bound(new_pos);                
            }
        }

        else if(fitness[di] = best_fitness){

            double k = 2 * (rand() / (RAND_MAX + 1.0)) - 1; // k: rnadom[-1, 1]
            double miniiiii = pow(10, -50);

            for(int d = 0; d < dimension; d++){
                new_pos = k * (abs(position[di][d]-worst_sol[d])/((fitness[di]-worst_fitness) + miniiiii));
                new_position[di][d] = check_bound(new_pos);
            }
        }
        
    }

    return new_position;

}


void determination(vector<vec> new_sol){

    double new_fitness;


    for(int i = 0; i < population; i++){

        new_fitness = evaluation(new_sol[i]);

        if(new_fitness < fitness[i]){
            position[i] = new_sol[i];
            fitness[i] = new_fitness;
        }

        if(new_fitness < best_fitness){
            best_sol = new_sol[i];
            best_fitness = new_fitness;
        }

    }


}


/* IMPROVEMENT -------------------------------------------- */


double levy_flight(double alpha, double beta){

    double levy_move;
    double a, b, uu, u, v;

    a = tgamma(1.0+beta)*sin(M_PI*beta/2);
    b = tgamma(((1.0+beta)/2)*beta*pow(2, ((beta-1.0)/2)));
    uu = pow((a/b), (1/beta));
    u = normal_distribute(generator);
    v = normal_distribute(generator);

    levy_move = u/pow(fabs(v), (1/beta));

    return alpha*levy_move;

}

// scrounger + levy flight
vector<vec> new_update_position(double pd, double sd, double st, int iter_max){

    vector<vec> new_position(population, vec(dimension, 0.0));

    // ----- for update producer
    int producer_num = population * pd;
    double r1, r2;  // r1: alpha -> random number (0,1]
    double new_pos;
    float q;  // random number from normal distribution
    
    double new_best_f = 0.0;
    vec new_best_sol(dimension, 0.0);

    // ----- for update scounger
    vec A(dimension, 0.0);
    double a_num, move;

    // ----- for update danger sparrow
    int danger_num = population * sd;
    vector<int> random_list(population, 0.0);
    int di; // index after random


    /********************************** update producer ******/
    
    r2 = ((double) rand() / (RAND_MAX));
  
    if(r2 < st){
        for(int i = 0; i < producer_num; i++){
            r1 = ((double) rand() / (RAND_MAX));
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k]*exp(-i/(r1*iter_max));
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }
    else{
        for(int i = 0; i < producer_num; i++){
            q = unif(generator);
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k] + q;
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }


    /********************************** update scrounger ******/

    double lvstep;
    double my_alpha=1;
    double my_beta=1.5;
    // my_alpha: step control, my_beta: according to paper

    // if(dimension<=10)
    //     my_alpha = 0.1;
    // else
    //     my_alpha = 1.0;    


    for(int i = 0; i < (population-producer_num); i++){
        int scrounger = i + producer_num;
        double r_step = ((double) rand() / (RAND_MAX));

        // create A: random 1 or -1 -> matrix(1*d)
        for(int d = 0; d < dimension; d++){
            a_num = ((double) rand() / (RAND_MAX));
            a_num = floor(a_num*2)*2-1;
            A[d] = a_num;
        }

        // update starving spparow (ranking lower): need to find other place for food
        if(scrounger > (population/2)){
            q = unif(generator);

            for(int d = 0; d < dimension; d++){
                lvstep = levy_flight(my_alpha, my_beta);
                new_pos = position[scrounger][d] + lvstep;
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
        // update other scrounger
        else{
            move = 0.0;
            for(int d = 0; d < dimension; d++){
                move += abs(position[scrounger][d] - new_best_sol[d]) * (1.0/dimension) * A[d];
            }

            for(int d = 0; d < dimension; d++){
                new_pos = best_sol[d] + move;
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
    }


    /********************************** update sparrow which are aware of dager ******/
    
    iota(random_list.begin(), random_list.end(), 0);
    random_shuffle(random_list.begin(), random_list.end());

    for(int i = 0; i < danger_num; i++){

        di = random_list[i];

        if(fitness[di] > best_fitness){; 
            
            for(int d = 0; d < dimension; d++){
                q = unif(generator);
                new_pos = best_sol[d] + q*abs((position[di][d]-best_sol[d]));
                new_position[di][d] = check_bound(new_pos);                
            }
        }

        else if(fitness[di] = best_fitness){

            double k = 2 * (rand() / (RAND_MAX + 1.0)) - 1; // k: rnadom[-1, 1]
            double miniiiii = pow(10, -50);

            for(int d = 0; d < dimension; d++){
                new_pos = k * (abs(position[di][d]-worst_sol[d])/((fitness[di]-worst_fitness) + miniiiii));
                new_position[di][d] = check_bound(new_pos);
            }
        }
        
    }

    return new_position;

}

// danger move + levy flight
vector<vec> new_update_position_2(double pd, double sd, double st, int iter_max){

    vector<vec> new_position(population, vec(dimension, 0.0));

    // ----- for update producer
    int producer_num = population * pd;
    double r1, r2;  // r1: alpha -> random number (0,1]
    double new_pos;
    float q;  // random number from normal distribution
    
    double new_best_f = 0.0;
    vec new_best_sol(dimension, 0.0);

    // ----- for update scounger
    vec A(dimension, 0.0);
    double a_num, move;

    // ----- for update danger sparrow
    int danger_num = population * sd;
    vector<int> random_list(population, 0.0);
    int di; // index after random


    /********************************** update producer ******/
    
    r2 = ((double) rand() / (RAND_MAX));

    if(r2 < st){
        for(int i = 0; i < producer_num; i++){
            r1 = ((double) rand() / (RAND_MAX));
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k]*exp(-i/(r1*iter_max));
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }
    else{
        for(int i = 0; i < producer_num; i++){
            q = unif(generator);
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k] + q;
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }


    /********************************** update scrounger ******/

    for(int i = 0; i < (population-producer_num); i++){
        int scrounger = i + producer_num;

        // create A: random 1 or -1 -> matrix(1*d)
        for(int d = 0; d < dimension; d++){
            a_num = ((double) rand() / (RAND_MAX));
            a_num = floor(a_num*2)*2-1;
            A[d] = a_num;
        }

        // update starving spparow (ranking lower): need to find other place for food
        if(scrounger > (population/2)){
            q = unif(generator);

            for(int d = 0; d < dimension; d++){
                new_pos = q * exp((worst_sol[d]-position[scrounger][d])/pow(scrounger, 2));
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
        // update other scrounger
        else{
            move = 0.0;
            for(int d = 0; d < dimension; d++){
                move += abs(position[scrounger][d] - new_best_sol[d]) * (1.0/dimension) * A[d];
            }

            for(int d = 0; d < dimension; d++){
                new_pos = best_sol[d] + move;
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
    }


    /********************************** update sparrow which are aware of dager ******/
    
    iota(random_list.begin(), random_list.end(), 0);
    random_shuffle(random_list.begin(), random_list.end());

    double lvstep;
    double my_alpha=1.0, my_beta=1.5;

    for(int i = 0; i < danger_num; i++){

        di = random_list[i];

        if(fitness[di] > best_fitness){; 
            
            for(int d = 0; d < dimension; d++){
                q = unif(generator);
                new_pos = best_sol[d] + q*abs((position[di][d]-best_sol[d]));
                new_position[di][d] = check_bound(new_pos);                
            }
        }

        else if(fitness[di] = best_fitness){

            double k = 2 * (rand() / (RAND_MAX + 1.0)) - 1; // k: rnadom[-1, 1]
            double miniiiii = pow(10, -50);

            for(int d = 0; d < dimension; d++){
                lvstep = levy_flight(my_alpha, my_beta);
                new_pos = position[di][d] + lvstep;
                new_position[di][d] = check_bound(new_pos);
            }
        }
        
    }

    return new_position;

}

// [scrounger & danger move] + levy flight
vector<vec> new_update_position_3(double pd, double sd, double st, int iter_max){

    vector<vec> new_position(population, vec(dimension, 0.0));

    // ----- for update producer
    int producer_num = population * pd;
    double r1, r2;  // r1: alpha -> random number (0,1]
    double new_pos;
    float q;  // random number from normal distribution
    
    double new_best_f = 0.0;
    vec new_best_sol(dimension, 0.0);

    // ----- for update scounger
    vec A(dimension, 0.0);
    double a_num, move;

    // ----- for update danger sparrow
    int danger_num = population * sd;
    vector<int> random_list(population, 0.0);
    int di; // index after random


    /********************************** update producer ******/
    
    r2 = ((double) rand() / (RAND_MAX));

    if(r2 < st){
        for(int i = 0; i < producer_num; i++){
            r1 = ((double) rand() / (RAND_MAX));
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k]*exp(-i/(r1*iter_max));
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }
    else{
        for(int i = 0; i < producer_num; i++){
            q = unif(generator);
            for(int k = 0; k < dimension; k++){
                new_pos = position[i][k] + q;
                new_position[i][k] = check_bound(new_pos);
            }
            double eval_temp = evaluation(new_position[i]);
            if(eval_temp < new_best_f){
                new_best_f = eval_temp;
                new_best_sol = new_position[i];
            }
        }
    }


    /********************************** update scrounger ******/

    double lvstep;
    double my_alpha=1;
    double my_beta=1.5;
    // my_alpha: step control, my_beta: according to paper

    for(int i = 0; i < (population-producer_num); i++){
        int scrounger = i + producer_num;

        // create A: random 1 or -1 -> matrix(1*d)
        for(int d = 0; d < dimension; d++){
            a_num = ((double) rand() / (RAND_MAX));
            a_num = floor(a_num*2)*2-1;
            A[d] = a_num;
        }

        // update starving spparow (ranking lower): need to find other place for food
        if(scrounger > (population/2)){
            q = unif(generator);

            for(int d = 0; d < dimension; d++){
                lvstep = levy_flight(my_alpha, my_beta);
                new_pos = position[scrounger][d] + lvstep;
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
        // update other scrounger
        else{
            move = 0.0;
            for(int d = 0; d < dimension; d++){
                move += abs(position[scrounger][d] - new_best_sol[d]) * (1.0/dimension) * A[d];
            }

            for(int d = 0; d < dimension; d++){
                new_pos = best_sol[d] + move;
                new_position[scrounger][d] = check_bound(new_pos);
            }
        }
    }


    /********************************** update sparrow which are aware of dager ******/
    
    iota(random_list.begin(), random_list.end(), 0);
    random_shuffle(random_list.begin(), random_list.end());

    for(int i = 0; i < danger_num; i++){

        di = random_list[i];

        if(fitness[di] > best_fitness){; 
            
            for(int d = 0; d < dimension; d++){
                q = unif(generator);
                new_pos = best_sol[d] + q*abs((position[di][d]-best_sol[d]));
                new_position[di][d] = check_bound(new_pos);                
            }
        }

        else if(fitness[di] = best_fitness){

            double k = 2 * (rand() / (RAND_MAX + 1.0)) - 1; // k: rnadom[-1, 1]
            double miniiiii = pow(10, -50);

            for(int d = 0; d < dimension; d++){
                lvstep = levy_flight(my_alpha, my_beta);
                new_pos = position[di][d] + lvstep;
                new_position[di][d] = check_bound(new_pos);
            }
        }
        
    }

    return new_position;

}