#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

long int nt=180000000000   # 時間ステップ数
double dt=0.0000005                # 時間刻み幅
int nx0=10                   # xの全体のセル数
int ny0=10                   # yのセル数
int dx=0.2               # xの刻み幅
int qx = 5
double dy=0.2                 # yの刻み幅
int qy = 5
int lbound=1                  # 仮想セル数
int nx=nx0*qx+2*lbound           # xの総セル数
int ny=ny0*qy+2*lbound           # yの総セル数
int in_diffusion = 1;          #拡散項を入れるか入れないか.考慮する場合は1,しない場合は0

double dlta=1.5
double du1_x = 0.3
double du1_y = dlta*du1_x
double du2_x = 7.5*du1_x
double du2_y = dlta*du2_x

double u1[nx][ny];
double u2[nx][ny];
//先に関数を宣言

//0-1の乱数を返す.
double uniform(void);
void setup();
void bound(double[nx][ny]);
void u_cal();
void laplacian1(double[nx][ny]);
void laplacian2(double[nx][ny]);
double R1(double,double);
double R2(double,double);
void boundary_left(double[nx][ny]);
void boundary_right(double[nx][ny]);
void boundary_upper(double[nx][ny]);
void boundary_lower(double[nx][ny]);
//setup内でu1,u2を初期化する際に用いる.
void initarray(void);
//配列Aを配列Bにコピーする
void copyArray(double[nx][ny] A,double[nx][ny] B);
 

int main(void){
    srand(time(NULL));
    time_t t = time(NULL);
    setup();
    
    for(int k = 0; k < nt; k++){
        u_cal();
        if(k+1 % every_outnum == 0){
            print("nt__________%011d\n",k+1);
        }
    }
    print("%d\n",time(NULL) - t);
    return 0;
}

double uniform(void){
    return rand() / ((double)RAND_MAX);
}

void initarray(void){
    for(int i=0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            u1[i][j] = uniform() - 0.5;
            u2[i][j] = uniform() - 0.5;
        }
    }
}

void setup(){
    initarray();

    u1 = bound(u1);
    u2 = bound(u2);
}

double R1(double a, double b){
    return a - a*a*a - b;
}

double R2(double a, double b){
    return 3*a - 2*b;
}

void boundary_left(double[nx][ny] u){
    for(int j = 0; j < ny; j++) u[lbound-1][j]=u[lbound][j];
}

void boundary_right(double[nx][ny] u){
    for(int j = 0; j < ny; j++) u[nx-lbound][j]=u[nx-lbound-1][j];
}

void boundary_lower(double[nx][ny] u){
    for(int i = 0; i < nx; i++) u[i][lbound-1]=u[i][lbound];
}

void boundary_upper(double[nx][ny] u){
    for(int i = 0; i < nx; i++) u[i][ny-lbound]=u[i][ny-lbound-1];
}

void boundary(double[nx][ny] u){
    double tmp[nx][ny];
    copyArray(u, tmp);
    tmp = boundary_upper(u);
    tmp = boundary_lower(u);
    tmp = boundary_right(u);
    tmp = boundary_left(u);
    u = tmp;
}

void u_cal(){
    double x_numerator1,x_numerator2,y_numerator1,y_numerator2;
    for(int i = lbound; i < nx - lbound; i++){
        for(int j = lbound; j < ny - lbound; j++){
            x_numerator1=du1_x*((u1[i+1][j]-u1[i][j])-(u1[i][j]-u1[i-1][j]));
            y_numerator1=du1_y*((u1[i][j+1]-u1[i][j])-(u1[i][j]-u1[i][j-1]));
            x_numerator2=du2_x*((u2[i+1][j]-u2[i][j])-(u2[i][j]-u2[i-1][j]));
            y_numerator2=du2_y*((u2[i][j+1]-u2[i][j])-(u2[i][j]-u2[i][j-1]));
            
            next_u1 =u1[i][j] +  dt * (x_numerator1 + y_numerator1) / (dx*dx) + dt *in_diffusion* R1(u1[i][j],u2[i][j]);
            next_u2 =u2[i][j] +  dt * (x_numerator2 + y_numerator2) / (dx*dx) + dt *in_diffusion* R2(u1[i][j],u2[i][j]);

            u1[i][j] = next_u1;
            u2[i][j] = next_u2;
        }
    }
    bound(u1);
    bound(u2);
}

void copyArray(double[nx][ny] A, double[nx][ny] B){
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            B[i][j] = A[i][j];
        }
    }
}