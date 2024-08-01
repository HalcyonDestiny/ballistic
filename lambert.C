#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const int N = 3;
const double PI = 3.14159265358979323846;
const double EPS = 23.45*PI/180;  // rad
const double MU_SUN = 132712440018.0;  // км^3/с^2
const double MU_EARTH = 398600.0;  // км^3/с^2
const double R_0 = 6571;  // км
const double A = 128000000;
const double ACC = 0.000001;
const double I = 51.6*PI/180;
double modulo_r_0, modulo_r_1, phi, t_m, l, delta_m;
double p, e;

double modulo(double []);
double *vect_mult(double [], double []);
double scalar_mult(double [], double []);
double delta(double );
double epsilon(double );
double g(double, int);
double d_g(double, int);
double newton_root(double, double, int);
double lambert();
int sign(double);

int main(){
    lambert();
    return 0;
}
double modulo(double array[]){
    double result;
    result = sqrt(pow(array[0], 2) + pow(array[1], 2) + pow(array[2], 2));
    return result;
}

double *vect_mult(double array1[], double array2[]){
    double *result = (double*)malloc(sizeof (double) * N);
    result[0] = array1[1]*array2[2] - array1[2]*array2[1];
    result[1] = array1[2]*array2[0] - array1[0]*array2[2];
    result[2] = array1[0]*array2[1] - array1[1]*array2[0];
    return result;
}

double scalar_mult(double array1[], double array2[]){
    double result;
    result = array1[0]*array2[0] + array1[1]*array2[1] + array1[2]*array2[2];
    return result;
}

double delta(double a){
    double result;
    result = 2*asin(sqrt((modulo_r_0 + modulo_r_1 - l)/(4*a)));
    return result;
}

double epsilon(double a){
    double result;
    result = 2*asin(sqrt((modulo_r_0 + modulo_r_1 + l)/(4*a)));
    return result;
}

double g(double a, int T_PER){
    double result;
    result = (PI + sign(t_m-T_PER)*(epsilon(a)-sin(epsilon(a)) - PI) - sign(sin(phi))*(delta(a)-sin(delta(a))))*sqrt((pow(a, 3))/(MU_SUN));
    return result - T_PER;
}

double  d_g(double a, int T_PER){
    double t, z, d_eps, d_sin_eps, d_delta, d_sin_delta, d_sqrt_a, sqrt_a, result;
    t = (modulo_r_0 + modulo_r_1 + l) / (a);
    z = (modulo_r_0 + modulo_r_1 - l) / (a);
    sqrt_a = sqrt((pow(a, 3))/(MU_SUN));
    d_eps = -(t)/(2*a*sqrt(t)*sqrt(1-t/4));
    d_sin_eps = cos(epsilon(a))*d_eps;
    d_delta = -(z)/(2*a*sqrt(z)*sqrt(1-z/4));
    d_sin_delta = cos(delta(a))*d_sin_delta;
    d_sqrt_a = sqrt_a*3/(2*a);
    result = (PI*d_sqrt_a + sign(t_m-T_PER)*(sqrt_a*d_eps-d_sqrt_a*epsilon(a)-(sqrt_a*d_sin_eps-d_sqrt_a*sin(epsilon(a))) - PI*d_sqrt_a) - sign(sin(phi))*(sqrt_a*d_delta-d_sqrt_a*delta(a)-(sqrt_a*d_sin_delta-d_sqrt_a*sin(delta(a)))));
    return result;
}

double newton_root(double x0, double eps, int t_per){
    double x1, f0, f1, df0;
    // double iter;
    // iter = 0;
    do {
        f0 = g(x0, t_per);
        df0 = d_g(x0, t_per);
        x1 = x0 - f0/df0;
        f1 = g(x1, t_per);
        x0 = x1;
        // iter++;
    } while (fabs(f1) > eps);
    return x1;
}

double lambert(){
    double r_0[N], r_1[N], V_0[N], V_1[N], V_0S[N], V_1S[N], V_0_hyp[N], V_1_hyp[N], V_0_inf[N];
    double jd1, jd2, t_per, modulo_V_0, modulo_V_1, t_par, a_res, etta, RAAN1, RAAN2, phi1, phi2, V_pi, sigma, p_inf, e_inf, true_anom, AOP1, AOP2, V_round, delta_V;
    double* ro;

    FILE* file_r_0;
    FILE* file_r_1;
    FILE* file_V_0;
    FILE* file_V_1;
    file_r_0 = fopen("data/r_1.txt", "r");
    file_r_1 = fopen("data/r_2.txt", "r");
    file_V_0 = fopen("data/V_1.txt", "r");
    file_V_1 = fopen("data/V_2.txt", "r");
    for(int i = 0; i < N; i++){
        fscanf(file_r_0, "%lf", &r_0[i]);
        fscanf(file_r_1, "%lf", &r_1[i]);
        fscanf(file_V_0, "%lf", &V_0[i]);
        fscanf(file_V_1, "%lf", &V_1[i]);
    }

    jd1 = 2464932.347509600;
    jd2 = 2465386.479085800;
    t_per = (jd2 - jd1)*(86400);

    modulo_r_0 = modulo(r_0);  // км
    modulo_r_1 = modulo(r_1);  // км
    modulo_V_0 = modulo(V_0);  // км/c
    modulo_V_1 = modulo(V_1);  // км/c

    ro = vect_mult(r_0, r_1);  // км^2

    if(sign(ro[2]) >= 0){
        phi = acos((scalar_mult(r_0, r_1)) / (modulo_r_0 * modulo_r_1));
    }
    else{
        phi = 2*PI - acos((scalar_mult(r_0, r_1)) / (modulo_r_0 * modulo_r_1));
    }

    l = sqrt(pow(modulo_r_0, 2) + pow(modulo_r_1, 2) - 2*modulo_r_0*modulo_r_1*cos(phi));  // км

    t_par = (pow(modulo_r_0 + modulo_r_1 + l, 1.5) - pow(modulo_r_0 + modulo_r_1 - l, 1.5)*sign(sin(phi))) / (6*sqrt(MU_SUN));

    delta_m = 2*asin(sqrt((modulo_r_0 + modulo_r_1 - l)/(modulo_r_0 + modulo_r_1 + l)));

    t_m = (pow(modulo_r_0 + modulo_r_1 + l, 1.5)*(PI - sign(sin(phi)) * (delta_m - sin(delta_m))))/(8*sqrt(MU_SUN));

    a_res = newton_root(A, ACC, t_per);  // км

    p = (modulo_r_0*modulo_r_1*sin(phi/2.0)*sin(phi/2.0))/(a_res*pow(sin((epsilon(a_res)-delta(a_res))/(2)), 2));  // км

    e = sqrt(1-p/a_res);

    for(int i = 0; i < N; i++){
        V_0S[i] = (sqrt(MU_SUN*p))/(modulo_r_0*modulo_r_1*sin(phi))*(r_1[i]-(1-(modulo_r_1*(1-cos(phi)))/(p))*r_0[i]);  // км/с
        V_1S[i] = (-r_0[i]+(1-(modulo_r_0*(1-cos(phi)))/(p))*r_1[i]) * (sqrt(MU_SUN*p))/(modulo_r_0*modulo_r_1*sin(phi));  // км/с
        V_0_hyp[i] = V_0S[i] - V_0[i];
        V_1_hyp[i] = V_1S[i] - V_1[i];
    }

    V_0_inf[0] = V_0_hyp[0];
    V_0_inf[1] = V_0_hyp[1]*cos(EPS) - V_0_hyp[2]*sin(EPS);
    V_0_inf[2] = V_0_hyp[1]*sin(EPS) + V_0_hyp[2]*cos(EPS);

    if(sign(V_0_inf[0]) >= 0){
        etta = acos((V_0_inf[1])/(sqrt(V_0_inf[0]*V_0_inf[0] + V_0_inf[1]*V_0_inf[1])));  // rad
    }
    else{
        etta = -acos((V_0_inf[1])/(sqrt(V_0_inf[0]*V_0_inf[0] + V_0_inf[1]*V_0_inf[1])));  // rad
    }

    RAAN1 = -etta + acos(1/tan(I)*(V_0_inf[2])/(sqrt(V_0_inf[0]*V_0_inf[0] + V_0_inf[1]*V_0_inf[1])));  // rad
    RAAN2 = -etta - acos(1/tan(I)*(V_0_inf[2])/(sqrt(V_0_inf[0]*V_0_inf[0] + V_0_inf[1]*V_0_inf[1])));  // rad

    phi1 = sign(V_0_inf[2])*acos((V_0_inf[0]*cos(RAAN1) + V_0_inf[1]*sin(RAAN1))/(modulo(V_0_inf)));  // rad
    phi2 = sign(V_0_inf[2])*acos((V_0_inf[0]*cos(RAAN2) + V_0_inf[1]*sin(RAAN2))/(modulo(V_0_inf)));  // rad

    V_pi = sqrt(2*MU_EARTH/R_0 + pow(modulo(V_0_inf), 2));  // км/с

    sigma = R_0*V_pi; // км^2/с

    p_inf = sigma*sigma / MU_EARTH;
    e_inf = sqrt(1+sigma*sigma*pow(modulo(V_0_inf), 2)/(MU_EARTH*MU_EARTH));

    true_anom = acos(-1/e_inf);  // rad

    AOP1 = phi1 - true_anom;
    AOP2 = phi2 - true_anom;

    V_round = sqrt(MU_EARTH/R_0);  // км/с

    delta_V = V_pi - V_round;
    // for(int i = 0; i < N; i++){
    //     printf("%lf\n", V_0_hyp[i]);
    // }
    printf("%.15lf", modulo(V_1S));
    return modulo(V_0S);
}

int sign(double num){
    int result;
    if(num == 0){
        return 0;
    }
    else{
        result = num / abs(num);
        return result;
    }
}
