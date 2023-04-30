#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/times.h>
#include <vector>
#include <cfloat>
#include <numeric>
#include <fstream>

using std::ofstream;

extern "C"
{
#include "f2c.h"
#include "clapack.h"
}

using namespace std;

static const int np = 100000; // np target points and np source points
static const int np3 = np * 3; // The dimension of the matrix, it is the same as n, but the type of n is integer (used for call Lapack)
static const int kmax = 200; // Max Lanczos iteration
static const double k_BT = 1;
static const double eta = 1;
static const double a = 0.1;
static const double pi = 3.14159265358979323846;
static const double Coeff1 = 1.0 / (6.0 * pi * a);
static const double Coeff5 = 1.0 / (8.0 * pi);
static const double Coeff2 = 2.0 * a * a / 3.0;
static const double Coeff3 = 3.0 / (32.0 * a);
static const double Coeff4 = 9.0 / (32.0 * a);
static const double eps = 1e-5;

static const int P = 8; // Order of Far-field approximation
static const int Pflat = (P + 1) * (P + 1) * (P + 1);
static const int N0 = 1000;
static const double sq_theta = 0.36; // theta = 0.6

//=================== Define some structs and variables =========================
struct vec_3d
{
    double val[3];
};

struct xyz // Particle coordinates (physical)
{
    size_t size;
    double* x;
    double* y;
    double* z;
    size_t* index;
    size_t* old_index;

    xyz(size_t numpars_in) : size(numpars_in)
    {
        x = new double[size];
        y = new double[size];
        z = new double[size];
        index = new size_t[size];
        old_index = new size_t[size];
    }

    ~xyz()
    {
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] index;
        delete[] old_index;
    }
};

struct panel
{
    size_t members[2];
    double xinterval[2];
    double yinterval[2];
    double zinterval[2];
    double xc; // Panel center x coordinate
    double yc; // Panel center y coordinate
    double zc; // Panel center z coordinate
    vector<size_t> children;
    double MAC; // r^2 / theta^2
    double moments[3][Pflat];
    int moment_flag;
    double t1[P + 1]; // Interpolation points in x direction
    double t2[P + 1];
    double t3[P + 1];

    // Initialization
    panel() : xc(0.0), yc(0.0), zc(0.0), MAC(0.0), moment_flag(0)
    {
        memset(members, 0, sizeof(members));
        memset(xinterval, 0, sizeof(xinterval));
        memset(yinterval, 0, sizeof(yinterval));
        memset(zinterval, 0, sizeof(zinterval));
        memset(moments, 0, sizeof(moments));
        memset(t1, 0, sizeof(t1));
        memset(t2, 0, sizeof(t2));
        memset(t3, 0, sizeof(t3));

        children.reserve(8);
    }
};

struct xyz particles(np);

vector<panel> tree;
vector<size_t> leaf;

size_t node_count = 0;

double* lambda[3]; // Updated in each Lanczos iteration

vec_3d* velo;
vec_3d* velo_true;

double zTDz;
double zTDz_a;
double xyzminmax[6];

int torder = 0;

//=================== treecode part =========================
long getTickCount()
{
    tms tm;
    return times(&tm);
}

double minval(double* x)
{
    double MinVal = x[0];
    for (int i = 1; i < np; i++) {
        if (MinVal > x[i])
            MinVal = x[i];
    }

    MinVal = MinVal;

    return MinVal;
}

double maxval(double* x)
{
    double MaxVal = x[0];
    for (int i = 1; i < np; i++) {
        if (MaxVal < x[i])
            MaxVal = x[i];
    }

    MaxVal = MaxVal;

    return MaxVal;
}

void build_tree_init()
{
    panel temp_panel;

    // Indices of particles belonging to the panel
    temp_panel.members[0] = 0;
    temp_panel.members[1] = np - 1;

    // Interval defining the panel
    temp_panel.xinterval[0] = xyzminmax[0];
    temp_panel.xinterval[1] = xyzminmax[1];
    temp_panel.yinterval[0] = xyzminmax[2];
    temp_panel.yinterval[1] = xyzminmax[3];
    temp_panel.zinterval[0] = xyzminmax[4];
    temp_panel.zinterval[1] = xyzminmax[5];

    temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
    temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
    temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);

    double xL = temp_panel.xinterval[1] - temp_panel.xinterval[0];
    double yL = temp_panel.yinterval[1] - temp_panel.yinterval[0];
    double zL = temp_panel.zinterval[1] - temp_panel.zinterval[0];

    double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
    temp_panel.MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    tree.push_back(temp_panel);
    node_count = 1;
}

void Swap(size_t i, size_t j)
{
    if (i == j)
        return;

    double x = particles.x[i];
    double y = particles.y[i];
    double z = particles.z[i];
    size_t index = particles.index[i];
    size_t old_index = particles.old_index[i];

    double temp_lambda[3];
    for (int index = 0; index < 3; index++)
        temp_lambda[index] = lambda[index][i];

    particles.x[i] = particles.x[j];
    particles.y[i] = particles.y[j];
    particles.z[i] = particles.z[j];
    particles.index[i] = particles.index[j];
    particles.old_index[i] = particles.old_index[j];

    for (int index = 0; index < 3; index++)
        lambda[index][i] = lambda[index][j];

    particles.x[j] = x;
    particles.y[j] = y;
    particles.z[j] = z;
    particles.index[j] = index;
    particles.old_index[j] = old_index;

    for (int index = 0; index < 3; index++)
        lambda[index][j] = temp_lambda[index];
}

void split_2(size_t panel_index, int split_code)
{
    panel child[2];

    /*
     -------------------
     |        |        |
     |        |        |
     |Child 0 |Child 1 |
     |        |        |
     |        |        |
     -----------------------------> axis A
     start     mid      end
     */

    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];
    double tp_y0 = tree[panel_index].yinterval[0];
    double tp_y1 = tree[panel_index].yinterval[1];
    double tp_z0 = tree[panel_index].zinterval[0];
    double tp_z1 = tree[panel_index].zinterval[1];

    for (int i = 0; i < 2; i++) {
        child[i].xinterval[0] = tp_x0;
        child[i].xinterval[1] = tp_x1;
        child[i].yinterval[0] = tp_y0;
        child[i].yinterval[1] = tp_y1;
        child[i].zinterval[0] = tp_z0;
        child[i].zinterval[1] = tp_z1;
    }

    double xL = tp_x1 - tp_x0;
    double yL = tp_y1 - tp_y0;
    double zL = tp_z1 - tp_z0;

    double* intervalA[2] = {NULL, NULL};
    double* coordA = NULL;
    double startpointA = 0.0, midpointA = 0.0, endpointA = 0.0;

    if (split_code == 4) { // XYZ = 100, A is X
        xL *= 0.5;
        intervalA[0] = child[0].xinterval;
        intervalA[1] = child[1].xinterval;
        coordA = particles.x;
        startpointA = tp_x0;
        endpointA = tp_x1;
    }
    else if (split_code == 2) { // XYZ = 010, A is Y
        yL *= 0.5;
        intervalA[0] = child[0].yinterval;
        intervalA[1] = child[1].yinterval;
        coordA = particles.y;
        startpointA = tp_y0;
        endpointA = tp_y1;
    }
    else if (split_code == 1) { // XYZ = 001, A is Z
        zL *= 0.5;
        intervalA[0] = child[0].zinterval;
        intervalA[1] = child[1].zinterval;
        coordA = particles.z;
        startpointA = tp_z0;
        endpointA = tp_z1;
    }

    midpointA = 0.5 * (startpointA + endpointA);

    // Child 0 ends with mid point on axis A
    intervalA[0][1] = midpointA;

    // Child 1 begins with mid point on axis A
    intervalA[1][0] = midpointA;

    double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
    double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2
    
    for (int i = 0; i < 2; i++) {
        child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
        child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
        child[i].zc = 0.5 * (child[i].zinterval[0] + child[i].zinterval[1]);
        child[i].MAC = MAC;
    }

    vector<size_t> v[2];
    size_t start = tree[panel_index].members[0];
    size_t end = tree[panel_index].members[1];
    size_t* addr_table = new size_t[end - start + 1];

    size_t index;
    for (index = start; index <= end; index++) {
        particles.index[index] = index;
        addr_table[index - start] = index;
        
        if (coordA[index] <= midpointA)
            v[0].push_back(index);
        else
            v[1].push_back(index);
    }

    size_t seq = start;
    for (size_t j = 0; j < 2; j++) {
        size_t size = v[j].size();

        if (size >= 1) {
            for (size_t k = 0; k < size; k++) {
                if (k == 0)
                    child[j].members[0] = seq;
                if (k == size - 1)
                    child[j].members[1] = seq;

                index = v[j][k];
                /*
                 // This is very slow
                 size_t pos;
                 for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++) {
                 if (particles.index[pos] == index)
                 break;
                 }
                 Swap(pos, seq);
                 */
                // This uses an address table
                size_t pos = addr_table[index - start];
                size_t out = particles.index[seq];
                Swap(pos, seq);
                addr_table[index - start] = seq;
                addr_table[out - start] = pos;
                
                seq++;
            }

            node_count++;
            tree[panel_index].children.push_back(node_count - 1);
            tree.push_back(child[j]);
            v[j].clear();
        }
    }

    delete[] addr_table;
}

void split_4(size_t panel_index, int XYZ_flag)
{
    panel child[4];

    /*
     ^ axis B
     |
     end -------------------
     |        |        |
     |Child 2 |Child 3 |
     mid |------------------
     |        |        |
     start |Child 0 |Child 1 |
     -----------------------------> axis A
     start     mid      end
     */

    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];
    double tp_y0 = tree[panel_index].yinterval[0];
    double tp_y1 = tree[panel_index].yinterval[1];
    double tp_z0 = tree[panel_index].zinterval[0];
    double tp_z1 = tree[panel_index].zinterval[1];

    for (int i = 0; i < 4; i++) {
        child[i].xinterval[0] = tp_x0;
        child[i].xinterval[1] = tp_x1;
        child[i].yinterval[0] = tp_y0;
        child[i].yinterval[1] = tp_y1;
        child[i].zinterval[0] = tp_z0;
        child[i].zinterval[1] = tp_z1;
    }

    double xL = tp_x1 - tp_x0;
    double yL = tp_y1 - tp_y0;
    double zL = tp_z1 - tp_z0;

    double* intervalA[4] = {NULL, NULL, NULL, NULL};
    double* intervalB[4] = {NULL, NULL, NULL, NULL};
    double* coordA = NULL;
    double* coordB = NULL;
    double startpointA = 0.0, endpointA = 0.0, midpointA = 0.0;
    double startpointB = 0.0, endpointB = 0.0, midpointB = 0.0;

    switch (XYZ_flag) {
        case 6: // XYZ = 110, A is X and B is Y
            xL *= 0.5;
            yL *= 0.5;
            intervalA[0] = child[0].xinterval;
            intervalB[0] = child[0].yinterval;
            intervalA[1] = child[1].xinterval;
            intervalB[1] = child[1].yinterval;
            intervalA[2] = child[2].xinterval;
            intervalB[2] = child[2].yinterval;
            intervalA[3] = child[3].xinterval;
            intervalB[3] = child[3].yinterval;
            coordA = particles.x;
            coordB = particles.y;
            startpointA = tp_x0;
            endpointA = tp_x1;
            startpointB = tp_y0;
            endpointB = tp_y1;
            break;

        case 3: // XYZ = 011, A is Y and B is Z
            yL *= 0.5;
            zL *= 0.5;
            intervalA[0] = child[0].yinterval;
            intervalB[0] = child[0].zinterval;
            intervalA[1] = child[1].yinterval;
            intervalB[1] = child[1].zinterval;
            intervalA[2] = child[2].yinterval;
            intervalB[2] = child[2].zinterval;
            intervalA[3] = child[3].yinterval;
            intervalB[3] = child[3].zinterval;
            coordA = particles.y;
            coordB = particles.z;
            startpointA = tp_y0;
            endpointA = tp_y1;
            startpointB = tp_z0;
            endpointB = tp_z1;
            break;

        case 5: // XYZ = 101, A is Z and B is X
            zL *= 0.5;
            xL *= 0.5;
            intervalA[0] = child[0].zinterval;
            intervalB[0] = child[0].xinterval;
            intervalA[1] = child[1].zinterval;
            intervalB[1] = child[1].xinterval;
            intervalA[2] = child[2].zinterval;
            intervalB[2] = child[2].xinterval;
            intervalA[3] = child[3].zinterval;
            intervalB[3] = child[3].xinterval;
            coordA = particles.z;
            coordB = particles.x;
            startpointA = tp_z0;
            endpointA = tp_z1;
            startpointB = tp_x0;
            endpointB = tp_x1;
            break;

        default:
            break;
    }

    midpointA = 0.5 * (startpointA + endpointA);
    midpointB = 0.5 * (startpointB + endpointB);

    // Child 0 ends with mid point on axis A, and ends with mid point on axis B
    intervalA[0][1] = midpointA;
    intervalB[0][1] = midpointB;

    // Child 1 begins with mid point on axis A, and ends with mid point on axis B
    intervalA[1][0] = midpointA;
    intervalB[1][1] = midpointB;

    // Child 2 ends with mid point on axis A, and begins with mid point on axis B
    intervalA[2][1] = midpointA;
    intervalB[2][0] = midpointB;

    // Child 3 begins with mid point on axis A, and begins with mid point on axis B
    intervalA[3][0] = midpointA;
    intervalB[3][0] = midpointB;

    double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
    double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    for (int i = 0; i < 4; i++) {
        child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
        child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
        child[i].zc = 0.5 * (child[i].zinterval[0] + child[i].zinterval[1]);
        child[i].MAC = MAC;
    }

    vector<size_t> v[4];
    size_t start = tree[panel_index].members[0];
    size_t end = tree[panel_index].members[1];
    size_t* addr_table = new size_t[end - start + 1];

    size_t index;
    for (index = start; index <= end; index++) {
        particles.index[index] = index;
        addr_table[index - start] = index;

        if (coordA[index] <= midpointA && coordB[index] <= midpointB)
            v[0].push_back(index);
        else if (coordA[index] > midpointA && coordB[index] <= midpointB)
            v[1].push_back(index);
        else if (coordA[index] <= midpointA && coordB[index] > midpointB)
            v[2].push_back(index);
        else if (coordA[index] > midpointA && coordB[index] > midpointB)
            v[3].push_back(index);
    }

    size_t seq = start;
    for (size_t j = 0; j < 4; j++) {
        size_t size = v[j].size();

        if (size >= 1) {
            for (size_t k = 0; k < size; k++) {
                if (k == 0)
                    child[j].members[0] = seq;
                if (k == size - 1)
                    child[j].members[1] = seq;

                index = v[j][k];
                /*
                 // This is very slow
                 size_t pos;
                 for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++) {
                 if (particles.index[pos] == index)
                 break;
                 }
                 Swap(pos, seq);
                 */
                // This uses an address table
                size_t pos = addr_table[index - start];
                size_t out = particles.index[seq];
                Swap(pos, seq);
                addr_table[index - start] = seq;
                addr_table[out - start] = pos;

                seq++;
            }

            node_count++;
            tree[panel_index].children.push_back(node_count - 1);
            tree.push_back(child[j]);
            v[j].clear();
        }
    }

    delete[] addr_table;
}

void split_8(size_t panel_index)
{
    panel child[8];

    /*
     ^ axis y
     |
     end -------------------
     |        |        |
     |Child 2 |Child 3 |
     mid |------------------
     |        |        |
     start |Child 0 |Child 1 |
     -----------------------------> axis x (lower z level)
     start     mid      end
     
     ^ axis y
     |
     end -------------------
     |        |        |
     |Child 6 |Child 7 |
     mid |------------------
     |        |        |
     start |Child 4 |Child 5 |
     -----------------------------> axis x (upper z level)
     start     mid      end
     */

    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];
    double tp_y0 = tree[panel_index].yinterval[0];
    double tp_y1 = tree[panel_index].yinterval[1];
    double tp_z0 = tree[panel_index].zinterval[0];
    double tp_z1 = tree[panel_index].zinterval[1];

    double xL = 0.5 * (tp_x1 - tp_x0);
    double yL = 0.5 * (tp_y1 - tp_y0);
    double zL = 0.5 * (tp_z1 - tp_z0);

    double midpointx = 0.5 * (tp_x0 + tp_x1);
    double midpointy = 0.5 * (tp_y0 + tp_y1);
    double midpointz = 0.5 * (tp_z0 + tp_z1);

    child[0].xinterval[0] = tp_x0;
    child[0].xinterval[1] = midpointx;
    child[0].yinterval[0] = tp_y0;
    child[0].yinterval[1] = midpointy;
    child[0].zinterval[0] = tp_z0;
    child[0].zinterval[1] = midpointz;

    child[1].xinterval[0] = midpointx;
    child[1].xinterval[1] = tp_x1;
    child[1].yinterval[0] = tp_y0;
    child[1].yinterval[1] = midpointy;
    child[1].zinterval[0] = tp_z0;
    child[1].zinterval[1] = midpointz;

    child[2].xinterval[0] = tp_x0;
    child[2].xinterval[1] = midpointx;
    child[2].yinterval[0] = midpointy;
    child[2].yinterval[1] = tp_y1;
    child[2].zinterval[0] = tp_z0;
    child[2].zinterval[1] = midpointz;

    child[3].xinterval[0] = midpointx;
    child[3].xinterval[1] = tp_x1;
    child[3].yinterval[0] = midpointy;
    child[3].yinterval[1] = tp_y1;
    child[3].zinterval[0] = tp_z0;
    child[3].zinterval[1] = midpointz;

    child[4].xinterval[0] = tp_x0;
    child[4].xinterval[1] = midpointx;
    child[4].yinterval[0] = tp_y0;
    child[4].yinterval[1] = midpointy;
    child[4].zinterval[0] = midpointz;
    child[4].zinterval[1] = tp_z1;

    child[5].xinterval[0] = midpointx;
    child[5].xinterval[1] = tp_x1;
    child[5].yinterval[0] = tp_y0;
    child[5].yinterval[1] = midpointy;
    child[5].zinterval[0] = midpointz;
    child[5].zinterval[1] = tp_z1;

    child[6].xinterval[0] = tp_x0;
    child[6].xinterval[1] = midpointx;
    child[6].yinterval[0] = midpointy;
    child[6].yinterval[1] = tp_y1;
    child[6].zinterval[0] = midpointz;
    child[6].zinterval[1] = tp_z1;

    child[7].xinterval[0] = midpointx;
    child[7].xinterval[1] = tp_x1;
    child[7].yinterval[0] = midpointy;
    child[7].yinterval[1] = tp_y1;
    child[7].zinterval[0] = midpointz;
    child[7].zinterval[1] = tp_z1;

    double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
    double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

    for (int i = 0; i < 8; i++) {
        child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
        child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
        child[i].zc = 0.5 * (child[i].zinterval[0] + child[i].zinterval[1]);
        child[i].MAC = MAC;
    }

    vector<size_t> v[8];
    size_t start = tree[panel_index].members[0];
    size_t end = tree[panel_index].members[1];
    size_t* addr_table = new size_t[end - start + 1];

    size_t index;
    for (index = start; index <= end; index++) {
        particles.index[index] = index;
        addr_table[index - start] = index;

        if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
            particles.z[index] <= midpointz)
            v[0].push_back(index);
        else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
                 particles.z[index] <= midpointz)
            v[1].push_back(index);
        else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
                 particles.z[index] <= midpointz)
            v[2].push_back(index);
        else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
                 particles.z[index] <= midpointz)
            v[3].push_back(index);
        else if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
                 particles.z[index] > midpointz)
            v[4].push_back(index);
        else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
                 particles.z[index] > midpointz)
            v[5].push_back(index);
        else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
                 particles.z[index] > midpointz)
            v[6].push_back(index);
        else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
                 particles.z[index] > midpointz)
            v[7].push_back(index);
    }

    size_t seq = start;
    for (size_t j = 0; j < 8; j++) {
        size_t size = v[j].size();

        if (size >= 1) {
            for (size_t k = 0; k < size; k++) {
                if (k == 0)
                    child[j].members[0] = seq;
                if (k == size - 1)
                    child[j].members[1] = seq;

                index = v[j][k];
                /*
                 // This is very slow
                 size_t pos;
                 for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++) {
                 if (particles.index[pos] == index)
                 break;
                 }
                 Swap(pos, seq);
                 */
                // This uses an address table
                size_t pos = addr_table[index - start];
                size_t out = particles.index[seq];
                Swap(pos, seq);
                addr_table[index - start] = seq;
                addr_table[out - start] = pos;

                seq++;
            }

            node_count++;
            tree[panel_index].children.push_back(node_count - 1);
            tree.push_back(child[j]);
            v[j].clear();
        }
    }

    delete[] addr_table;
}

void split_tree_node(size_t panel_index)
{
    double xL = tree[panel_index].xinterval[1] - tree[panel_index].xinterval[0];
    double yL = tree[panel_index].yinterval[1] - tree[panel_index].yinterval[0];
    double zL = tree[panel_index].zinterval[1] - tree[panel_index].zinterval[0];
    double L_max = xL;

    if (yL > L_max)
        L_max = yL;

    if (zL > L_max)
        L_max = zL;

    int XYZ_flag = 0;
    const double ratio = sqrt(2.0);

    if (xL * ratio > L_max)
        XYZ_flag += 4;

    if (yL * ratio > L_max)
        XYZ_flag += 2;

    if (zL * ratio > L_max)
        XYZ_flag += 1;

    switch (XYZ_flag) {
        case 1: // XYZ = 001, split along Z
        case 2: // XYZ = 010, split along Y
        case 4: // XYZ = 100, split along X
#ifdef TREE_STAT
            split2_times++;
#endif
            split_2(panel_index, XYZ_flag);
            break;

        case 3: // XYZ = 011, split along Y and Z
        case 5: // XYZ = 101, split along Z and X
        case 6: // XYZ = 110, split along X and Y
#ifdef TREE_STAT
            split4_times++;
#endif
            split_4(panel_index, XYZ_flag);
            break;

        case 7: // XYZ = 111, split along X, Y, and Z
#ifdef TREE_STAT
            split8_times++;
#endif
            split_8(panel_index);
            break;

        default:
            break;
    }
}

void build_tree_3D_Recursive(size_t panel_index, int level)
{
#ifdef TREE_STAT
    if (level > max_level)
        max_level = level;

    level_cnt[level]++;
#endif

    size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;
    if (n >= (size_t)N0) {
        split_tree_node(panel_index);
        
        for (size_t i = 0; i < tree[panel_index].children.size(); i++) {
            size_t panel_index_new = tree[panel_index].children[i];
            build_tree_3D_Recursive(panel_index_new, level + 1);
        }
    }
    else {
        leaf.push_back(panel_index);

#ifdef TREE_STAT
        double tp_x0 = tree[panel_index].xinterval[0];
        double tp_x1 = tree[panel_index].xinterval[1];
        double tp_y0 = tree[panel_index].yinterval[0];
        double tp_y1 = tree[panel_index].yinterval[1];
        double tp_z0 = tree[panel_index].zinterval[0];
        double tp_z1 = tree[panel_index].zinterval[1];

        double xL = tp_x1 - tp_x0;
        double yL = tp_y1 - tp_y0;
        double zL = tp_z1 - tp_z0;

        double diameter = sqrt(xL * xL + yL * yL + zL * zL);
        if (diameter > max_leaf_diameter)
            max_leaf_diameter = diameter;
        if (diameter < min_leaf_diameter)
            min_leaf_diameter = diameter;

        int pars = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;
        if (pars > max_leaf_pars)
            max_leaf_pars = pars;
        if (pars < min_leaf_pars)
            min_leaf_pars = pars;

        double Lmax = xL;
        if (yL > Lmax)
            Lmax = yL;
        if (zL > Lmax)
            Lmax = zL;

        double Lmin = xL;
        if (yL < Lmin)
            Lmin = yL;
        if (zL < Lmin)
            Lmin = zL;

        double leaf_ratio = Lmax / Lmin;
        if (leaf_ratio > max_leaf_ratio)
            max_leaf_ratio = leaf_ratio;
        if (leaf_ratio < min_leaf_ratio)
            min_leaf_ratio = leaf_ratio;
#endif
    }
}

void Panel_Moment_B(size_t panel_index, double m[][Pflat])
{
    // Intput : panel_index
    // Output : m: moments for panel_index^th panel

    double t1[P + 1];
    double t2[P + 1];
    double t3[P + 1];

    int i, j, k, kk;
    int a1exactind, a2exactind, a3exactind;

    for (i = 0; i < P + 1; i++) {
        t1[i] = tree[panel_index].t1[i];
        t2[i] = tree[panel_index].t2[i];
        t3[i] = tree[panel_index].t3[i];
    }

    kk = -1;
    for (i = 0; i < P + 1; i++) {
        for (j = 0; j < P + 1; j++) {
            for (k = 0; k < P + 1; k++) {
                kk++;
                m[0][kk] = 0;
                m[1][kk] = 0;
                m[2][kk] = 0;;
            }
        }
    }

    double w1i[P + 1];
    double dj[P + 1];
    dj[0] = 0.5;
    dj[P] = 0.5;
    for (j = 1; j < P; j++)
        dj[j] = 1.0;

    for (j = 0; j < P + 1; j++)
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];

    double a1i[P + 1];
    double a2j[P + 1];
    double a3k[P + 1];

    double x, y, z;
    double dx, dy, dz;
    double SumA1;
    double SumA2;
    double SumA3;
    double D;
    double s;

    size_t tp0 = tree[panel_index].members[0];
    size_t tp1 = tree[panel_index].members[1];
    size_t tp_j;

    for (tp_j = tp0; tp_j <= tp1; tp_j++) {
        x = particles.x[tp_j];
        y = particles.y[tp_j];
        z = particles.z[tp_j];

        a1exactind = -1;
        a2exactind = -1;
        a3exactind = -1;

        SumA1 = 0.0;
        SumA2 = 0.0;
        SumA3 = 0.0;

        for (j = 0; j < P + 1; j++) {
            dx = x - t1[j];
            dy = y - t2[j];
            dz = z - t3[j];

            if (fabs(dx) <= DBL_MIN)
                a1exactind = j;
            else {
                a1i[j] = w1i[j] / dx;
                SumA1 += a1i[j];
            }

            if (fabs(dy) <= DBL_MIN)
                a2exactind = j;
            else {
                a2j[j] = w1i[j] / dy;
                SumA2 += a2j[j];
            }

            if (fabs(dz) <= DBL_MIN)
                a3exactind = j;
            else {
                a3k[j] = w1i[j] / dz;
                SumA3 += a3k[j];
            }
        }

        if (a1exactind > -1) {
            SumA1 = 1.0;
            for (j = 0; j < P + 1; j++) a1i[j] = 0.0;
            a1i[a1exactind] = 1.0;
        }

        if (a2exactind > -1) {
            SumA2 = 1.0;
            for (j = 0; j < P + 1; j++) a2j[j] = 0.0;
            a2j[a2exactind] = 1.0;
        }

        if (a3exactind > -1) {
            SumA3 = 1.0;
            for (j = 0; j < P + 1; j++) a3k[j] = 0.0;
            a3k[a3exactind] = 1.0;
        }

        D = 1.0 / (SumA1 * SumA2 * SumA3);

        kk = -1;
        for (i = 0; i < P + 1; i++) {
            for (j = 0; j < P + 1; j++) {
                for (k = 0; k < P + 1; k++) {
                    kk++;
                    s = a1i[i] * a2j[j] * a3k[k] * D;
                    m[0][kk] += s * lambda[0][tp_j];
                    m[1][kk] += s * lambda[1][tp_j];
                    m[2][kk] += s * lambda[2][tp_j];
                }
            }
        }
    }
}

vec_3d Call_Treecode(double x, double y, double z, int panel_index)
{
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;

    double R2, R, Rinv, Rinv3, Rinv5;

    double* dx = new double[torder + 1];
    double* dy = new double[torder + 1];
    double* dz = new double[torder + 1];

    double temp_moments[3];
    double G;
    double DD, CC;

    for (int i = 0; i < torder + 1; i++) {
        dx[i] = x - tree[panel_index].t1[i];
        dy[i] = y - tree[panel_index].t2[i];
        dz[i] = z - tree[panel_index].t3[i];
    }

    int kk = -1;
    double xx, yy, zz;
    double sq_xx, sq_yy, sq_zz;
    for (int i = 0; i < torder + 1; i++) {
        xx = dx[i];
        sq_xx = xx * xx;
        for (int j = 0; j < torder + 1; j++) {
            yy = dy[j];
            sq_yy = yy * yy;
            for (int k = 0; k < torder + 1; k++) {
                zz = dz[k];
                sq_zz = zz * zz;
                kk++;

                temp_moments[0] = tree[panel_index].moments[0][kk];
                temp_moments[1] = tree[panel_index].moments[1][kk];
                temp_moments[2] = tree[panel_index].moments[2][kk];

                R2 = sq_xx + sq_yy + sq_zz;
                R = sqrt(R2);
                Rinv = 1.0 / R;
                Rinv3 = Rinv * Rinv * Rinv;
                Rinv5 = Rinv * Rinv * Rinv3;

                G = (Rinv + Coeff2 * Rinv3) * Coeff5;
                velocity.val[0] += temp_moments[0] * G;
                velocity.val[1] += temp_moments[1] * G;
                velocity.val[2] += temp_moments[2] * G;

                DD = (temp_moments[0] * xx + temp_moments[1] * yy + temp_moments[2] * zz) * Coeff5;
                CC = DD * Rinv3;
                velocity.val[0] += CC * xx;
                velocity.val[1] += CC * yy;
                velocity.val[2] += CC * zz;

                CC = -2.0 * a * a * DD * Rinv5;
                velocity.val[0] += CC * xx;
                velocity.val[1] += CC * yy;
                velocity.val[2] += CC * zz;
            }
        }
    }

    delete[] dx;
    delete[] dy;
    delete[] dz;

    return velocity;
}

vec_3d Call_Ds(size_t limit_1, size_t limit_2, size_t particle_index, double p_x, double p_y, double p_z)
{
    double xx, yy, zz;
    double sq_xx, sq_yy, sq_zz;
    double ff[3];
    double R2, R, Rinv, Rinv3, Rinv5;
    double DD, CC;
    double G;

    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;

    for (size_t jj = limit_1; jj <= limit_2; jj++) {
        ff[0] = lambda[0][jj];
        ff[1] = lambda[1][jj];
        ff[2] = lambda[2][jj];

        if (jj == particle_index) {
            velocity.val[0] += ff[0] * Coeff1;
            velocity.val[1] += ff[1] * Coeff1;
            velocity.val[2] += ff[2] * Coeff1;
        }
        else {
            xx = p_x - particles.x[jj];
            yy = p_y - particles.y[jj];
            zz = p_z - particles.z[jj];

            sq_xx = xx * xx;
            sq_yy = yy * yy;
            sq_zz = zz * zz;

            R2 = sq_xx + sq_yy + sq_zz;
            R = sqrt(R2);
            Rinv = 1.0 / R;
            Rinv3 = Rinv * Rinv * Rinv;
            Rinv5 = Rinv * Rinv * Rinv3;

            if (R < 2.0 * a) {
                G = (1.0 - Coeff4 * R) * Coeff1;
                velocity.val[0] += ff[0] * G;
                velocity.val[1] += ff[1] * G;
                velocity.val[2] += ff[2] * G;
                
                DD = ff[0] * xx + ff[1] * yy + ff[2] * zz;
                CC = DD * Rinv * Coeff3 * Coeff1;
                velocity.val[0] += CC * xx;
                velocity.val[1] += CC * yy;
                velocity.val[2] += CC * zz;
            }
            else {
                G = (Rinv + Coeff2 * Rinv3) * Coeff5;
                velocity.val[0] += ff[0] * G;
                velocity.val[1] += ff[1] * G;
                velocity.val[2] += ff[2] * G;

                DD = (ff[0] * xx + ff[1] * yy + ff[2] * zz) * Coeff5;
                CC = DD * Rinv3;
                velocity.val[0] += CC * xx;
                velocity.val[1] += CC * yy;
                velocity.val[2] += CC * zz;

                CC = -2.0 * a * a * DD * Rinv5;
                velocity.val[0] += CC * xx;
                velocity.val[1] += CC * yy;
                velocity.val[2] += CC * zz;
            }
        }
    }

    return velocity;
}

vec_3d Compute_SUM(size_t particle_index, size_t panel_index)
{
    // Input :
    //         particle_index
    //         panel_index
    // Output :
    //         velocity in 3D
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;

    double p_x = 0.0;
    double p_y = 0.0;
    double p_z = 0.0;
    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;
    size_t limit_1;
    size_t limit_2;
    double R_sq;

    limit_1 = tree[panel_index].members[0];
    limit_2 = tree[panel_index].members[1];

    p_x = particles.x[particle_index];
    p_y = particles.y[particle_index];
    p_z = particles.z[particle_index];

    xc = tree[panel_index].xc;
    yc = tree[panel_index].yc;
    zc = tree[panel_index].zc;

    double tpx = p_x - xc;
    double tpy = p_y - yc;
    double tpz = p_z - zc;

    R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

    if (tree[panel_index].MAC < R_sq)
        velocity = Call_Treecode(p_x, p_y, p_z, panel_index);
    else {
        if (limit_2 - limit_1 < N0) // Otherwise, if cluster is a leaf, use direct sum
            velocity = Call_Ds(limit_1, limit_2, particle_index, p_x, p_y, p_z);
        else { // Otherwise, if cluster is not a leaf, look at children
            size_t length = tree[panel_index].children.size();
            for (size_t i = 0; i < length; i++) {
                size_t index = tree[panel_index].children[i];
                vec_3d temp_result = Compute_SUM(particle_index, index);
                velocity.val[0] += temp_result.val[0];
                velocity.val[1] += temp_result.val[1];
                velocity.val[2] += temp_result.val[2];
            }
        }
    }

    return velocity;
}

void Cluster_Chev_Points(size_t tree_size)
{
    double h = pi / P;
    double t[P + 1] = {0.0};
    for (int i = 0; i < P + 1; i++)
        t[i] = cos(i * h); // Chebyshev interpolation points [-1, 1]

    double x1, x2, y1, y2, z1, z2;
    size_t tree_index;

    for (tree_index = 0; tree_index < tree_size; tree_index++) {
        x1 = tree[tree_index].xinterval[0];
        x2 = tree[tree_index].xinterval[1];
        y1 = tree[tree_index].yinterval[0];
        y2 = tree[tree_index].yinterval[1];
        z1 = tree[tree_index].zinterval[0];
        z2 = tree[tree_index].zinterval[1];

        for (int i = 0; i < P + 1; i++) { // Map to the cluster
            tree[tree_index].t1[i] = x1 + (t[i] + 1.0) * 0.5 * (x2 - x1);
            tree[tree_index].t2[i] = y1 + (t[i] + 1.0) * 0.5 * (y2 - y1);
            tree[tree_index].t3[i] = z1 + (t[i] + 1.0) * 0.5 * (z2 - z1);
        }
    }
}

//============ SLMD part ==========================================
inline double inner_prod(const vector<double>& v1, const vector<double>& v2)
{
    return inner_product(begin(v1), end(v1), begin(v2), 0.0);
}

//=====================
inline void normalize(vector<double>& vec)
{
    double scalar = 1.0 / sqrt(inner_prod(vec, vec));

    int n = vec.size();
    for (int i = 0; i < n; i++)
        vec[i] *= scalar;
}

//=====================
void schmidt_orth(vector<double>& uorth, const vector<vector<double> >& u)
{
    size_t n = uorth.size();

    for (size_t k = 0; k < u.size(); k++) {
        double innprod = inner_prod(uorth, u[k]);

        for (size_t i = 0; i < n; i++)
            uorth[i] -= innprod * u[k][i];
    }
}

//=====================
double vec_norm(double vec[], size_t n)
{
    double temp_r = 0.0;
    for (size_t i = 0; i < n; i++)
        temp_r += vec[i] * vec[i];

     return sqrt(temp_r);
}

//================== Compute ds using gymbutas idea ==================
int compute_Ds()
{
    double temp_x, temp_y, temp_z;
    double px, py, pz;
    double xx, yy, zz;
    double sq_xx, sq_yy, sq_zz;
    double ff[3], sum[3];
    double R2, R, Rinv, Rinv3, Rinv5;
    double DD, CC;
    double G;

    velo_true = new vec_3d[np];

    for (int i = 0; i < np; i++) {
        temp_x = particles.x[i];
        temp_y = particles.y[i];
        temp_z = particles.z[i];

        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;

        for (int j = 0; j < np; j++) {
            ff[0] = lambda[0][j];
            ff[1] = lambda[1][j];
            ff[2] = lambda[2][j];

            if (i == j) {
                sum[0] += ff[0] * Coeff1;
                sum[1] += ff[1] * Coeff1;
                sum[2] += ff[2] * Coeff1;
            }
            else {
                px = particles.x[j];
                py = particles.y[j];
                pz = particles.z[j];

                xx = temp_x - px;
                yy = temp_y - py;
                zz = temp_z - pz;

                sq_xx = xx * xx;
                sq_yy = yy * yy;
                sq_zz = zz * zz;

                R2 = sq_xx + sq_yy + sq_zz;
                R = sqrt(R2);
                Rinv = 1.0 / R;
                Rinv3 = Rinv * Rinv * Rinv;
                Rinv5 = Rinv * Rinv * Rinv3;

                if (R < 2.0 * a) {
                    G = (1.0 - Coeff4 * R) * Coeff1;
                    sum[0] += ff[0] * G;
                    sum[1] += ff[1] * G;
                    sum[2] += ff[2] * G;
                    
                    DD = ff[0] * xx + ff[1] * yy + ff[2] * zz;
                    CC = DD * Rinv * Coeff3 * Coeff1;
                    sum[0] += CC * xx;
                    sum[1] += CC * yy;
                    sum[2] += CC * zz;
                }
                else {
                    G = (Rinv + Coeff2 * Rinv3) * Coeff5;
                    sum[0] += ff[0] * G;
                    sum[1] += ff[1] * G;
                    sum[2] += ff[2] * G;

                    DD = (ff[0] * xx + ff[1] * yy + ff[2] * zz) * Coeff5;
                    CC = DD * Rinv3;
                    sum[0] += CC * xx;
                    sum[1] += CC * yy;
                    sum[2] += CC * zz;
                    
                    CC = -2.0 * a * a * DD * Rinv5;
                    sum[0] += CC * xx;
                    sum[1] += CC * yy;
                    sum[2] += CC * zz;
                }
            }
        }

        velo_true[i].val[0] = sum[0];
        velo_true[i].val[1] = sum[1];
        velo_true[i].val[2] = sum[2];
    }

    return 0;
}

//============== Compute SLDM with treecode ================
void Compute_SLDM_treecode()
{
    int n = np3; // Dimension of the matrix

    double* d = new double[kmax];
    double* e = new double[kmax];
    double* d_back = new double[kmax];
    double* e_back = new double[kmax];

    vector<vector<double> > u; // Lanczos vectors
    u.push_back(vector<double>(n, 0.0));

    vector<double> vk(n, 0.0);

    double alphak = 0.0;
    double betak = 0.0;

    vector<double> uk(n); // Initial q_1
    double v_norm = 0;

    for (int i = 0; i < np; i++) {
        uk[i*3 + 0] = lambda[0][i]; // Initial q_1 is just the vector, i.e. \sqrt{D}vec
        uk[i*3 + 1] = lambda[1][i];
        uk[i*3 + 2] = lambda[2][i];
        v_norm += lambda[0][i] * lambda[0][i] + lambda[1][i] * lambda[1][i] + lambda[2][i] * lambda[2][i];
    }

    v_norm = sqrt(v_norm); // Compute the norm of the incoming vector
    normalize(uk); // Normalize q_1
    u.push_back(uk);

    double* u_p = new double[n];
    double* diff = new double[n];
    double* Q_v_sq_lambda_vT_e = new double[n];
    double* SLDM_DS_q = new double[n];

    for (int i = 0; i< n; i++) {
        diff[i] = 0;
        Q_v_sq_lambda_vT_e[i] = 0;
        u_p[i] = 1.0; // Previous step vector
    }

    //==================== Start Lanczos iteration ======================
    for (int k = 1; k <= kmax; k++) {
        fill(vk.begin(), vk.end(), 0.0);

        for (int i = 0; i < np; i++) {
            lambda[0][i] = uk[i*3 + 0];
            lambda[1][i] = uk[i*3 + 1];
            lambda[2][i] = uk[i*3 + 2];
        }

        // treecode
        size_t tree_size = tree.size();
        for (size_t i = 1; i < tree_size; i++) // Skip roots
            Panel_Moment_B(i, tree[i].moments);

        for (int particle_index = 0; particle_index < np; particle_index++)
            velo[particle_index] = Compute_SUM(particle_index, 0);

        for (int i = 0; i < np; i++) {
            vk[i*3 + 0] = velo[i].val[0];
            vk[i*3 + 1] = velo[i].val[1];
            vk[i*3 + 2] = velo[i].val[2];
        }

        /*
        compute_Ds(); // Computes velo_true

        for (int i = 0; i < np; i++) {
            vk[i*3 + 0] = velo_true[i].val[0];
            vk[i*3 + 1] = velo_true[i].val[1];
            vk[i*3 + 2] = velo_true[i].val[2];
        }
         */

        alphak = inner_prod(u.back(), vk);
        d[k - 1] = alphak; // Diagonal elements of an approximated tridiagonal matrix
        d_back[k - 1] = alphak;

        for (int i = 0; i < n; i++)
            uk[i] = vk[i] - betak * u[k - 1][i] - alphak * u[k][i];

        schmidt_orth(uk, u);

        betak = vec_norm(uk.data(), uk.size());
        e[k - 1] = betak; // Subdiagonal elements of an approximated tridiagonal matrix
        e_back[k - 1] = betak;

        normalize(uk);
        u.push_back(uk);

        if (k > 1) {
            integer kk;
            integer ldz = k;
            double z[kmax * kmax] = {0.0}; // ldz * n
            double work[2 * kmax - 2] = {0.0}; // Dimension (max(1, 2 * N - 2))
            integer info;

            kk = k;

            // Compute e-value and e-vector of T
            dsteqr_("I", &kk, d, e, z, &ldz, work, &info);

            // Compute SLDM: Q_v_sq_lambda_vT_e
            double v_sq_lambda_vT_e[kmax] = {0.0};
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++)
                    v_sq_lambda_vT_e[i] += z[j * k + i] * sqrt(d[j]) * z[k * j];
            }

            for (int i = 0; i < n; i++) {
                Q_v_sq_lambda_vT_e[i] = 0.0;
                for (int j = 0; j < k; j++)
                    Q_v_sq_lambda_vT_e[i] += u[j + 1][i] * v_sq_lambda_vT_e[j];

                Q_v_sq_lambda_vT_e[i] = Q_v_sq_lambda_vT_e[i] * v_norm;
            }

            // ======== Compute error y^Ty with z^TDz ==================
            double yTy = 0.0;

            //========================================================== 
            for (int i = 0; i < n; i++)
                diff[i] = Q_v_sq_lambda_vT_e[i] - u_p[i];

            if (vec_norm(diff, n) / vec_norm(u_p, n) < eps || k == kmax) {
                yTy = 0;
                for (int i = 0; i < n; i++)
                    yTy += Q_v_sq_lambda_vT_e[i] * Q_v_sq_lambda_vT_e[i];

                double Err = abs(yTy - zTDz) / abs(zTDz);
                cout << "k = " << k << endl;
                cout << "yTy is " <<  setprecision(15) << yTy << endl;
                cout << "Err_r = " << setprecision(15)  << Err << endl;
                cout << "Err_a = " << setprecision(15)  << abs(yTy - zTDz) << endl;
                cout << "difference between two iteration: r " << vec_norm(diff, n) / vec_norm(u_p, n) << endl;
                cout << "difference between two iteration: a " << vec_norm(diff, n) << endl;
                cout << "difference between two iteration is less than" << eps << endl;

                ofstream outdata;
                outdata.open("Rand_PVF0P12_Glambda_ModifyLambdaE2.txt");

                for (int i = 0; i < n; i++)
                    outdata << setprecision(25) << Q_v_sq_lambda_vT_e[i] << endl;

                outdata.close();

                //========== Measure error with the one computed through DS ======
                char point_data_file[64] = {0};
                sprintf(point_data_file, "./Rand_PVF012_k1000_Aug11_Glambda_2.txt");

                FILE *fp = fopen(point_data_file, "r");

                int counter = 0;
                double x1;
                while (fscanf(fp, "%lf", &x1) == 1) {
                    SLDM_DS_q[counter] = x1;

                    counter++;
                    if (counter == n)
                        break;
                }

                fclose(fp);

                for (int i = 0; i < n; i++)
                    diff[i] = Q_v_sq_lambda_vT_e[i] - SLDM_DS_q[i];

                cout << "difference with SLDM_DS: r " << vec_norm(diff, n) / vec_norm(SLDM_DS_q, n) << endl;
                cout << "difference with SLDM_DS: a " << vec_norm(diff, n) << endl;
                break;
            }

            // ========== Update and recover ===================
            for (int i = 0; i < n; i++)
                u_p[i] = Q_v_sq_lambda_vT_e[i]; // Update u_p
            
            for (int i = 0; i < kmax; i++) {
                d[i] = d_back[i]; // Recover d, e
                e[i] = e_back[i];
            }
        } // End if (k > 1)
    } // End Lanczos iteration

    delete[] d;
    delete[] e;
    delete[] d_back;
    delete[] e_back;
    delete[] diff;
    delete[] Q_v_sq_lambda_vT_e;
    delete[] u_p;
    delete[] SLDM_DS_q;
}

//================ main ========================
/* g++ -o testlapack SLDM.cpp liblapack.a libblas.a libf2c.a -fsanitize=address -fno-omit-frame-pointer -O0 -g */
int main()
{
    cout << " CASEIII number of beads is " << np << endl;
    cout << "theta = " << sqrt(sq_theta) << endl;
    cout << "P = " << P << endl;
    cout << "N_0 = "<< N0 << endl;
    cout << "eps = " << eps << endl;
    cout << "a = " << a << endl;
    torder = P;

    //============ The particles =============

    //=========  Read from a data file =========
    char point_data_file[64] = {0};
    sprintf(point_data_file, "./rand_%d.txt", np);

    FILE *fp = fopen(point_data_file, "r");
    if (fp == NULL) {
        cerr << "Cannot open point data file " << point_data_file << "!" << endl;
        return 1;
    }

    int counter = 0;
    double x1, x2, x3;
    while (fscanf(fp, "%lf %lf %lf", &x1, &x2, &x3) == 3) {
        particles.x[counter] = x1;
        particles.y[counter] = x2;
        particles.z[counter] = x3;
        particles.index[counter] = -1;
        particles.old_index[counter] = counter;

        counter++;
        if (counter == np)
            break;
    }

    fclose(fp);

    if (counter < np) {
        cerr << "Less lines of point data were read! counter = " << counter << endl;
        return 1;
    }

    //============= The vector =======================
    for (int i = 0; i < 3; i++)
        lambda[i] = new double[np];

    //=========== Read v from a data file ============
    sprintf(point_data_file, "./G_lambda_%d.txt", np);
    fp = fopen(point_data_file, "r");
    if (fp == NULL) {
        cerr << "Cannot open vector data file " << point_data_file << "!" << endl;
        return 1;
    }

    counter = 0;
    while (fscanf(fp, "%lf %lf %lf", &x1, &x2, &x3) == 3) {
        lambda[0][counter] = x1;
        lambda[1][counter] = x2;
        lambda[2][counter] = x3;

        counter++;
        if (counter == np)
            break;
    }

    fclose(fp);

    if (counter < np) {
        cerr << "Less lines of vector data were read! counter = " << counter << endl;
        return 1;
    }

    //========== Build tree: note that datas and lambda are permuted ==========================
    xyzminmax[0] = minval(particles.x);
    xyzminmax[1] = maxval(particles.x);
    xyzminmax[2] = minval(particles.y);
    xyzminmax[3] = maxval(particles.y);
    xyzminmax[4] = minval(particles.z);
    xyzminmax[5] = maxval(particles.z);

    long total_time, SLDM_cpu_time, treecode_cpu_time, ds_cpu_time;

    total_time = getTickCount();

    build_tree_init();
    build_tree_3D_Recursive(0, 0); // After this step, both particles and lambda has been permutted.
    size_t tree_size = tree.size();
    Cluster_Chev_Points(tree_size);

    velo = new vec_3d[np];
    for (size_t i = 1; i < tree_size; i++) // Skip roots
        Panel_Moment_B(i, tree[i].moments);

    for (int particle_index = 0; particle_index < np; particle_index++)
        velo[particle_index] = Compute_SUM(particle_index, 0);

    treecode_cpu_time = getTickCount() - total_time;

    cout << "treecode_cpu_time " << treecode_cpu_time << endl;

    //============= Compute z^TDz ===================
    ds_cpu_time = getTickCount();

    compute_Ds();

    ds_cpu_time = getTickCount() - ds_cpu_time;

    cout << "ds_cpu_time " << ds_cpu_time << endl;

    zTDz = 0.0;
    zTDz_a = 0.0;

    for (int i = 0; i < np; i++) {
        for (int j = 0; j < 3; j++)
            zTDz += velo_true[i].val[j] * lambda[j][i];
    }

    cout << "zTDz = " << setprecision(15) << zTDz << endl;

    double err2_ex = 0.0;
    double sum_d_ex = 0.0;
    double sum_n_ex = 0.0;
    for (int i = 0; i < np; i++) {
        for (int d = 0; d < 3; d++) {
            sum_n_ex += (velo[i].val[d] - velo_true[i].val[d]) * (velo[i].val[d] - velo_true[i].val[d]);
            sum_d_ex += velo_true[i].val[d] * velo_true[i].val[d];
        }
    }

    err2_ex = sqrt(sum_n_ex / sum_d_ex);

    cout << endl;
    cout << setprecision(25) << "treecode E_2_ex is " << err2_ex << endl;

    //============= Compute \sqrt{D}vec using SLDM with treecode =======
    SLDM_cpu_time = getTickCount();

    Compute_SLDM_treecode();

    SLDM_cpu_time = getTickCount() - SLDM_cpu_time;

    cout << "SLDM dse PVF 012 is " << SLDM_cpu_time << endl;

    //========= Delete pointers =======
    for (int i = 0; i < 3; i++)
        delete[] lambda[i];

    delete[] velo_true;
    delete[] velo;

    return 0;
}
