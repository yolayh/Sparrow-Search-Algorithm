#include "SSA.hpp"


// test 30 run, evaluation = dimension*10,000
// test 2, 10, 30 dimension

/*
    執行方式: 
        g++ ./SSA.cpp -o ssa
        ./ssa [test_function_number] [version] [dimension] [run]

        [test_function_number]
            0: all  (1, 2, 7, 4, 5, 6, 8, 9, 10, 11)
            1: Ackley
            2: Griewank
            3: Katsuuras
            4: Michalewicz
            5: Rosenbrock
            6: Schwefel2.26"
            7: bent cigar
            8: zakharov
            9: happy cat
            10: rastrigin
            11: hgbat  
        
        [version]
            0: origin
            1: time-variate st
            2: time-variate sd
            3: scrouger movement adapt levy flight
            4: aware danger movement adapt levy flight
            5: version 3 + version 4
            6: version 1 + version 3
            7: version 1 + version 4

*/

double select_version(int, int, int, int, int, double, double, double);
double ssa_search(int, int, int, int, double, double, double);
double ssa_search_v1(int, int, int, int, double, double, double);
double ssa_search_v2(int, int, int, int, double, double, double);
double ssa_search_v3(int, int, int, int, double, double, double);
double ssa_search_v4(int, int, int, int, double, double, double);
double ssa_search_v5(int, int, int, int, double, double, double);
double ssa_search_v6(int, int, int, int, double, double, double);
double ssa_search_v7(int, int, int, int, double, double, double);


int main(int argc, char *argv[]){

    const int t = atoi(argv[1]);
    const int version_input = atoi(argv[2]);
    const int dim_input = atoi(argv[3]);
    const int run = atoi(argv[4]);

    unsigned int seed;
    seed = time (NULL);
    srand(seed);
    cout << scientific << setprecision(8);
    
    int dim = dim_input;             // dimension
    
    /*parameter*/
    int pop = 100;           // population
    double PD_rate = 0.2;    // proportion of producer
    double SD_rate = 0.1;    // proportion of sparrow which are aware of danger (原本:0.1)
    double ST_rate = 0.8;    // safety threshold (0.5)
    /*--------*/

    vector<string> testing = {"Testing function:", "Ackley", "Griewank", "Bent Cigar", "Michalewicz", \
        "Rosenbrock", "Schwefel2.26", "Zakharov", "HappyCat", "Rastrigin", "HGBat"};
    
    vector<string> name_list = {"Testing function:", "Ackley", "Griewank", "Katsuuras", "Michalewicz", \
        "Rosenbrock", "Schwefel2.26", "Bent Cigar", "Zakharov", "HappyCat", "Rastrigin", "HGBat"};
    
    vector<int> testing_list = {1, 2, 7, 4, 5, 6, 8, 9, 10, 11};

    int evaluation = dim*10000; 
    double output;
    vec all_result(run, 0.0);
    double best_result=100000;

    switch(t){

        case 0:  // test all function

            for(int k = 0; k < 10; k++){
                best_result=100000;
                int function_number = testing_list[k];
                for(int i = 0; i < run; i++){
                    output = select_version(version_input, dim, evaluation, function_number, pop, PD_rate, SD_rate, ST_rate);
                    all_result[i] = output;
                    if(output < best_result)
                        best_result = output;
                }
                write_record(testing[k+1], all_result);
                cout << testing[k+1] << ", Avg: " << calculate_avg(all_result) << ", best: " << best_result << endl;
            }

            break;
        
        
      
        default:

            if(t>11){
                cout << "input must <= 11" << endl;
                break;
            }

            else{
                for(int i = 0; i < run; i++){
                    output = select_version(version_input, dim, evaluation, t, pop, PD_rate, SD_rate, ST_rate);
                    all_result[i] = output;
                    if(output < best_result)
                        best_result = output;
                }
                write_record(name_list[t], all_result);
                cout << name_list[t] << ", Avg: " << calculate_avg(all_result) << endl;

                break;
            }
    }

    return 0;
}


double select_version(int version_num, int dim, int eval, int test_func, int pop, double pd, double sd, double st){
    
    double r;

    switch (version_num){
        case 0:
            r = ssa_search(dim, eval, test_func, pop, pd, sd, st);
            break;

        case 1:
            r = ssa_search_v1(dim, eval, test_func, pop, pd, sd, st);
            break;
        
        case 2:
            r = ssa_search_v2(dim, eval, test_func, pop, pd, sd, st);
            break;
        
        case 3:
            r = ssa_search_v3(dim, eval, test_func, pop, pd, sd, st);
            break;
        
        case 4:
            r = ssa_search_v4(dim, eval, test_func, pop, pd, sd, st);
            break;
        
        case 5:
            r = ssa_search_v5(dim, eval, test_func, pop, pd, sd, st);
            break;

        case 6:
            r = ssa_search_v6(dim, eval, test_func, pop, pd, sd, st);
            break;

        case 7:
            r = ssa_search_v7(dim, eval, test_func, pop, pd, sd, st);
            break;
          
        default:
            break;
    }
    return r;
}

// dimension, evaluation_times, testing_function, pop(number of sparrow), pd(number of producer), sd(number to perceive danger), st(safe threshold)
// st 0.5~1.0

double ssa_search(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        new_solution = update_position(pd, sd, st, iter);
        determination(new_solution);

    }


    return best_fitness;

}


double ssa_search_v1(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        // double my_st = 0.2 + 0.7 * i / iter;    //increasing
        double my_st = 0.9 - 0.7 * i / iter;     //decreasing
        new_solution = update_position(pd, sd, my_st, iter);
        determination(new_solution);


    }


    return best_fitness;

}

double ssa_search_v2(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        // double my_sd = 0.2 + 0.7 * i / iter;    //increasing
        double my_sd = 0.9 - 0.7 * i / iter;     //decreasing
        new_solution = update_position(pd, my_sd, st, iter);
        determination(new_solution);

    }


    return best_fitness;
}


double ssa_search_v3(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        new_solution = new_update_position(pd, sd, st, iter);
        determination(new_solution);


    }


    return best_fitness;

}

double ssa_search_v4(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        new_solution = new_update_position_2(pd, sd, st, iter);
        determination(new_solution);

    }


    return best_fitness;

}


double ssa_search_v5(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        new_solution = new_update_position_3(pd, sd, st, iter);
        determination(new_solution);

    }


    return best_fitness;

}


double ssa_search_v6(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        // double my_st = 0.2 + 0.7 * i / iter;    //increasing
        double my_st = 0.9 - 0.7 * i / iter;     //decreasing
        new_solution = new_update_position(pd, sd, st, iter);
        determination(new_solution);

    }

    return best_fitness;

}


double ssa_search_v7(int dim, int eval, int test_func, int pop, double pd, double sd, double st){

    vector<vec> new_solution;

    initialization(dim, pop, test_func);
    int iter = eval / pop;

    for(int i = 0; i < iter; i++){
        
        ranking();
        // double my_st = 0.2 + 0.7 * i / iter;    //increasing
        double my_st = 0.9 - 0.7 * i / iter;     //decreasing
        new_solution = new_update_position_2(pd, sd, st, iter);
        determination(new_solution);

    }

    return best_fitness;

}

