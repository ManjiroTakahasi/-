#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

#define nx 53
#define ny 53

long int nt=180000000000;   // 時間ステップ数
double dt=0.0000005;                // 時間刻み幅
int nx0 = 10;
int ny0 = 10;
int dx=0.2;               // xの刻み幅
int qx = 5;
double dy=0.2;                 // yの刻み幅
const int qy = 5;
const int lbound=1;                  // 仮想セル数
//int nx = nx0*qx+2*lbound;
//int ny = ny0*qy+2*lbound;
int in_diffusion = 1;          //拡散項を入れるか入れないか.考慮する場合は1,しない場合は0
int every_outnum=100000;

double dlta=1.5;
double du1_x = 0.3;
double du1_y = 1.5*0.3;
double du2_x = 7.5*0.3;
double du2_y = 1.5*7.5*0.3;

double u1[nx][ny];
double u2[nx][ny];
//先に関数を宣言


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

    boundary(u1);
    boundary(u2);
}

double R1(double a, double b){
    return a - a*a*a - b;
}

double R2(double a, double b){
    return 3*a - 2*b;
}

void boundary_left(){
    for(int j = 0; j < ny; j++){
	    u1[lbound-1][j]=u1[lbound][j];
	    u2[lbound-1][j]=u2[lbound][j];
    }
}

void boundary_right(){
    for(int j = 0; j < ny; j++){
	    u1[nx-lbound][j]=u1[nx-lbound-1][j];
	    u2[nx-lbound][j]=u2[nx-lbound-1][j];
    }
}

void boundary_lower(){
    for(int i = 0; i < nx; i++){
	    u1[i][lbound-1]=u1[i][lbound];
	    u2[i][lbound-1]=u2[i][lbound];
    }
}

void boundary_upper(){
    for(int i = 0; i < nx; i++){
	    u1[i][ny-lbound]=u1[i][ny-lbound-1];
	    u2[i][ny-lbound]=u2[i][ny-lbound-1];
    }
}

void boundary(){
    boundary_upper();
    boundary_lower();
    boundary_right();
    boundary_left();
}

void u_cal(){
    double x_numerator1,x_numerator2,y_numerator1,y_numerator2;
    double next_u1,next_u2;
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
    boundary(u1);
    boundary(u2);
}

 

int main(void){
    srand(time(NULL));
    time_t t = time(NULL);
    setup();
    
    for(int k = 0; k < nt; k++){
	printf("k = %d\n",k);
        u_cal();
        if(k+1 % every_outnum  == 0){
            printf("nt__________%011d\n",k+1);
        }
    }
    printf("%ld\n",time(NULL) - t);
    return 0;
}
