#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <random>
#include <ctime>
#include <iomanip>

using namespace std;

typedef vector<double> vec;

// Test function
double ackley(int, vec);
double griewank(int, vec);
double katsuuras(int, vec);
double michalewicz(int, vec);
double rosenbrock(int, vec);
double schwefel(int, vec);
double bent_cigar(int, vec);
double zakharov(int, vec);
double happy_cat(int, vec);
double rastrigin(int, vec);
double hgbat(int, vec);


/* -------------------------------------------- */


double ackley(int d, vec a){

    /*
        min -32.768, max 32.768
        minimum = 0 at xi = 0
    */

    double quadratic, cos_sum, t;

    for(int i = 0; i < d; i++){
        t = a[i];
        quadratic += t * t;
        cos_sum += cos(2.0*M_PI*a[i]);
    }

    return -20.0 * exp(-0.2*sqrt(quadratic / d)) - exp(cos_sum / d) + 20.0 + M_E;

}


double griewank(int d, vec a){

    /*
        test domain |xi| <= 600
        min -600, max 600
        minimum = 0 at xi = 0
    */

    double quadratic = 0.0;
    double cos_multi = 1.0;

    for(int i = 0; i < d; i++){
        quadratic += pow(a[i], 2);
        cos_multi *= cos(a[i]/sqrt((double)(i+1)));
    }

    return quadratic/4000.0 - (cos_multi - 1.0);

}


double katsuuras(int d, vec a){

    /*
        test domain -100 <= xi <= 100
        minimum = 0 at xi = 0
    */
    
    double power2;
    double sum;
    double multi = 1.0;

    for(int i = 0; i < d; i++){
        sum = 0.0;
        for(int k = 0; k < 32; k++){
            power2 = pow(2, k+1);
            sum += round(power2*a[i]) / power2;
        }
        multi *= 1.0 + (1+i)*sum;
    }

    return multi;

}


double michalewicz(int d, vec a){

    /*
        domain: 0 <= xi <= pi
        minimum = -0.966n
    */

    double sum = 0.0;
    double sin_pow;

    for(int i = 0; i < d; i++){
        sin_pow = pow(sin((i+1)*pow(a[i], 2)/M_PI), 20);
        sum += sin(a[i])*sin_pow;
    }

    if(d==2)
        return  1.8013-sum;
    else if(d==10)
        return 9.66015-sum;
    else
        return -sum;

}


double rosenbrock(int d, vec a){

    /*
        domain: |xi| <= 10
        min -5, max 10
        minimum = 0 at xi = 1
    */

    double sum = 0.0;

    for(int i = 1; i < d; i++){
        sum += 100.0 * pow((a[i]-pow(a[i-1], 2)), 2);
        sum += pow((1-a[i-1]), 2);
    }

    return sum;
}


double schwefel(int d, vec a){

    /*
        domain: |xi| < 500
        min -500, max 500
        minimum = -12569.5 at xi = 420.9687 For d=3
        optimum: 418.9829 * d - sum = 0;
    */

    double sum = 0.0;
    

    for(int i = 0; i < d; i++){
        sum += a[i] * sin(sqrt(abs(a[i])));
    }

    return  418.9829*d - sum;
}


double bent_cigar(int d, vec a){

    /*  max=100.0, min=-100.0 optimal=0  */

    double sum = 0.0;
    for(int i = 0; i < d; i++){
        sum += a[i]*a[i];
    }

    return a[0]*a[0] + pow(10.0, 6)*sum;
}


double zakharov(int d, vec a){
    
    /*  max=10.0, min=-10.0, optimal=0  */

    double sum_1 = 0.0;
    double sum_2 = 0.0;
    for(int i = 0; i < d; i++){
        sum_1 += a[i] * a[i];
        sum_2 += 0.5 * (i+1) * a[i];
    }

    return sum_1 + pow(sum_2, 2) + pow(sum_2, 4);
}


double happy_cat(int d, vec a){

    /*  max 20.0, min -20.0  */

    double sum_1 = 0.0;
    double sum_2 = 0.0;
    for(int i = 0; i < d; i++){
        sum_1 += a[i] * a[i];
        sum_2 += a[i];
    }

    return pow(fabs(sum_1-double(d)), 0.25) + (0.5 * sum_1 + sum_2) / double(d) + 0.5;
}


double rastrigin(int d, vec a){

    /*  max 5.12, min -5.12  */

    double sum = 0.0;
    for(int i = 0; i < d; i++){
        sum += pow(a[i], 2) - (10.0 * cos(2.0 * M_PI * a[i]));
    }

    return sum + 10.0 * double(d);
}


double hgbat(int d, vec a){

    /*  max 15.0, min -15.0  */

    double sum_1 = 0.0;
    double sum_2 = 0.0;
    for(int i = 0; i < d; i++){
        sum_1 += a[i] * a[i];
        sum_2 += a[i];
    }

    return sqrt(fabs((sum_1 * sum_1) - (sum_2 * sum_2))) + (0.5 * sum_1 + sum_2) / double(d) + 0.5;
}