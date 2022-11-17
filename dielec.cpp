#include <math.h>
#include <omp.h>
#include <vector>

extern "C" {
/*
def dielectric(point1,point2):
    "return e*r*r 1/(e*r^2)"
    #print(point1)
    #print(point2)
    from scipy.spatial import distance as D
    r= D.euclidean(point1,point2)
    A= -8.5525
    B= e0 - A
    k= 7.7839
    l= 0.003627
    e= A + (B/(1+k*math.exp(-1.0*l*B*r)))
    return 1/(e*r*r)
*/


typedef struct {float x;float y;float z;} point;




float dielectric(point A, point B, point l_m, point r_m){
    float r = (A.x  - B.x + l_m.x - r_m.x)*(A.x  - B.x + l_m.x - r_m.x) + (A.y - B.y + l_m.y - r_m.y) *(A.y  -B.y + l_m.y -r_m.y) + (A.z -B.z + l_m.z - r_m.z)*(A.z - B.z + l_m.z -r_m.z);
    r = sqrtf(r);

 float cA = -8.5525;
 float e0 = 78.4;
 float cB = e0 -cA;
 float k= 7.7839;
 float l=0.003627;

    float e = cA + (cB / (1 + k * (exp(-1.0 * l *cB * r))));
    return 1/(e * r * r);
}


void calc(point *ligand, point *receptor, point ligand_mean, point receptor_mean,int ligand_n,int receptor_n, float *result[]){
    
    std::vector<std::pair<int,int>> N;
    for (int i=0;i++;i< ligand_n){
        for (int j=0;j++;j< receptor_n){
            N.emplace_back(std::make_pair(i,j));
        }
    }

    for (int i=0;i++;i<N.size()){
        result[N[i].first][N[i].second] =dielectric(ligand[N[i].first],receptor[N[i].second], ligand_mean, receptor_mean);
    }
}


void calc_wrap(float *l, float *r, float *lm,float *rm, int ln, int rn, float *result){


}

void test(float *l, int N, float *r){
    point *ligand;
    ligand= new point[N];
    
    for (int i =0;i<N;i++){
        ligand[i].x = l[3*i];
        ligand[i].y = l[3*i+1];
        ligand[i].z = l[3*i+2];
        r[3*i] = l[3*i];
        r[3*i+1] = l[3*i+1];
        r[3*i+2] = l[3*i +2];
    }
}

}