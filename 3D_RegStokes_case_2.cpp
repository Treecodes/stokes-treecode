
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>

// vector researve

using namespace std;

static const int P = 1; // order of Taylor approximation
static const double L = 21.3333 ; //size of box
static const int Pflat = (P + 1)*(P + 1)*(P + 1);
static const int N_cube = 60400; // N_cube points total
static const int N = 60400;
static const int N0 = 2000;
static const double sq_theta = 0.25; // theta = 0.5
static const double DEL = 0.3;
const bool UseSleep = false; // for testing memory usage purpose
int max_level = 0;


//**********************************************************//

struct vec_3d
{
    double val[6];
};

//**************//

struct xyz // particle coordinates (physical)
{
	double* x;
	double* y;
	double* z;
	size_t* index;
	size_t* old_index;
	size_t size;
	xyz(size_t N_cube_in)
	{
		size = N_cube_in;
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

//**************//

struct panel
{
	size_t members[2];
	double xinterval[2];
	double yinterval[2];
	double zinterval[2];
	double xc; // panel center x coordinate
	double yc; // panel center y coordinate
	double zc; // panel center z coordinate
	vector<size_t> children;
	double MAC; // r^2 / theta^2
	double moments[6][Pflat];
	int moment_flag;
    double t1[P + 1]; // interpolation points in x direction
    double t2[P + 1];
    double t3[P + 1];
	panel() // initialization
	{
		moment_flag = 0;
		members[0] = 0;
		members[1] = -1;
		for(size_t index = 0; index < 6; index ++)
		{
			for (size_t kk = 0; kk < Pflat + 1; kk++)
                moments[index][kk] = 0;
		}
        for (int i = 0; i < P + 1; i++)
        {
            t1[i] = 0.0;
            t2[i] = 0.0;
            t3[i] = 0.0;
        }
	}
};

//**************//
vector<panel> tree;
vector<size_t> leaf;

static size_t node_count = 0;

//*****************************************************************************//
long getTickCount()
{
    tms tm;
    return times(&tm);
}

//*****************************************************************************//
void build_tree_init()
{
	panel temp_panel;

	// indices of particles belonging to panel
	temp_panel.members[0] = 0;
	temp_panel.members[1] = N_cube - 1;

    temp_panel.xinterval[0] = -L/2 - 0.1234;//-8.1234; // interval defining the panel
    temp_panel.xinterval[1] = L/2 + 0.1122;//8.1111;
    temp_panel.yinterval[0] = -L/2 - 0.1234;//-8.1234;
    temp_panel.yinterval[1] = L/2 + 0.1122; //8.1111;
	temp_panel.zinterval[0] = -0.1234;
	temp_panel.zinterval[1] = 10;
    temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
    temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
    temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);
    temp_panel.MAC = (3 * 16 * 16/ 4) / sq_theta; // MAC = r^2 / theta^2

	tree.push_back(temp_panel);
	node_count = 1;
}

//*****************************************************************************//
void Swap(size_t i, size_t j, struct xyz &s)
{
	if (i == j)
		return;

	double x = s.x[i];
	double y = s.y[i];
	double z = s.z[i];
	size_t index = s.index[i];
	size_t old_index = s.old_index[i];

	s.x[i] = s.x[j];
	s.y[i] = s.y[j];
	s.z[i] = s.z[j];
	s.index[i] = s.index[j];
	s.old_index[i] = s.old_index[j];

	s.x[j] = x;
	s.y[j] = y;
	s.z[j] = z;
	s.index[j] = index;
	s.old_index[j] = old_index;
}
//*****************************************************************************//
void split_tree_node(size_t panel_index, struct xyz &particles)
{
	panel child[8];

	double tp_x0 = tree[panel_index].xinterval[0];
	double tp_x1 = tree[panel_index].xinterval[1];
	double tp_y0 = tree[panel_index].yinterval[0];
	double tp_y1 = tree[panel_index].yinterval[1];
	double tp_z0 = tree[panel_index].zinterval[0];
	double tp_z1 = tree[panel_index].zinterval[1];

	double midpointx = (tp_x0 + tp_x1) / 2.0;
	double midpointy = (tp_y0 + tp_y1) / 2.0;
    double midpointz = (tp_z0 + tp_z1) / 2.0;

	double xc0 = (tp_x0 + midpointx) / 2.0;
	double xc1 = (tp_x1 + midpointx) / 2.0;
	double yc0 = (tp_y0 + midpointy) / 2.0;
	double yc1 = (tp_y1 + midpointy) / 2.0;
	double zc0 = (tp_z0 + midpointz) / 2.0;
	double zc1 = (tp_z1 + midpointz) / 2.0;

	child[0].xinterval[0] = tp_x0;
	child[0].xinterval[1] = midpointx;
	child[0].yinterval[0] = tp_y0;
	child[0].yinterval[1] = midpointy;
	child[0].zinterval[0] = tp_z0;
	child[0].zinterval[1] = midpointz;
	child[0].xc = xc0;
	child[0].yc = yc0;
	child[0].zc = zc0;
	child[0].MAC = ((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[1].xinterval[0] = midpointx;
	child[1].xinterval[1] = tp_x1;
	child[1].yinterval[0] = tp_y0;
    child[1].yinterval[1] = midpointy;
	child[1].zinterval[0] = tp_z0;
	child[1].zinterval[1] = midpointz;
	child[1].xc = xc1;
	child[1].yc = yc0;
	child[1].zc = zc0;
	child[1].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[2].xinterval[0] = tp_x0;
	child[2].xinterval[1] = midpointx;
	child[2].yinterval[0] = midpointy;
	child[2].yinterval[1] = tp_y1;
	child[2].zinterval[0] = tp_z0;
	child[2].zinterval[1] = midpointz;
	child[2].xc = xc0;
	child[2].yc = yc1;
	child[2].zc = zc0;
    child[2].MAC = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[3].xinterval[0] = midpointx;
	child[3].xinterval[1] = tp_x1;
	child[3].yinterval[0] = midpointy;
	child[3].yinterval[1] = tp_y1;
	child[3].zinterval[0] = tp_z0;
	child[3].zinterval[1] = midpointz;
	child[3].xc = xc1;
	child[3].yc = yc1;
	child[3].zc = zc0;
    child[3].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[4].xinterval[0] = tp_x0;
	child[4].xinterval[1] = midpointx;
	child[4].yinterval[0] = tp_y0;
	child[4].yinterval[1] = midpointy;
    child[4].zinterval[0] = midpointz;
	child[4].zinterval[1] = tp_z1;
	child[4].xc = xc0;
	child[4].yc = yc0;
	child[4].zc = zc1;
	child[4].MAC = ((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	child[5].xinterval[0] = midpointx;
	child[5].xinterval[1] = tp_x1;
	child[5].yinterval[0] = tp_y0;
    child[5].yinterval[1] = midpointy;
	child[5].zinterval[0] = midpointz;
	child[5].zinterval[1] = tp_z1;
	child[5].xc = xc1;
	child[5].yc = yc0;
	child[5].zc = zc1;
	child[5].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	child[6].xinterval[0] = tp_x0;
	child[6].xinterval[1] = midpointx;
	child[6].yinterval[0] = midpointy;
	child[6].yinterval[1] = tp_y1;
	child[6].zinterval[0] = midpointz;
	child[6].zinterval[1] = tp_z1;
	child[6].xc = xc0;
	child[6].yc = yc1;
	child[6].zc = zc1;
	child[6].MAC = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	child[7].xinterval[0] = midpointx;
	child[7].xinterval[1] = tp_x1;
	child[7].yinterval[0] = midpointy;
	child[7].yinterval[1] = tp_y1;
	child[7].zinterval[0] = midpointz;
	child[7].zinterval[1] = tp_z1;
	child[7].xc = xc1;
	child[7].yc = yc1;
	child[7].zc = zc1;
	child[7].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	vector<size_t> v[8];
	size_t start = tree[panel_index].members[0];
	size_t end = tree[panel_index].members[1];
	size_t* addr_table = new size_t[end - start + 1];

	size_t index;
	for (index = start; index <= end; index++)
	{
		particles.index[index] = index;
		addr_table[index - start] = index;

		if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
			particles.z[index] <= midpointz)
			v[0].push_back(index);
		else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
				particles.z[index] <= midpointz )
			v[1].push_back(index);
		else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
				particles.z[index]<= midpointz)
			v[2].push_back(index);
		else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
				particles.z[index] <= midpointz)
			v[3].push_back(index);
		else if(particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
				particles.z[index] > midpointz )
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
	for (size_t j = 0; j < 8; j++)
	{
		size_t size = v[j].size();

		if (size >= 1)
		{
			for (size_t k = 0; k < size; k++)
			{
				if (k == 0)
					child[j].members[0] = seq;
				if (k == size - 1)
					child[j].members[1] = seq;

				index = v[j][k];
				/*
				// This is very slow
				size_t pos;
				for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++)
				{
					if (particles.index[pos] == index)
						break;
				}
				Swap(pos, seq, particles);
				*/
				// This uses an address table
				size_t pos = addr_table[index - start];
				size_t out = particles.index[seq];
				Swap(pos, seq, particles);
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

//*****************************************************************************//
void build_tree_3D_Recursive(size_t panel_index, struct xyz &particles, int level)
{
	if (level > max_level)
		max_level = level;
	
	size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;

	if (n >= (size_t)N0)
   // if (level < fix_level)
	{
		split_tree_node(panel_index, particles);

		for (size_t i = 0; i < tree[panel_index].children.size(); i++)
		{
			size_t panel_index_new = tree[panel_index].children[i];
			build_tree_3D_Recursive(panel_index_new, particles, level + 1);
		}
	}
	else
		leaf.push_back(panel_index);
}

//*****************************************************************************//

double mypow(double x, int n)
{
    double result = 1.0;
    while (n > 0) {
        if (n & 1)
            result = result * x;
        n = n >> 1;
        x = x * x;
    }
    
    return result;
}

//*****************************************************************************//
void Panel_Moment_B(size_t panel_index, double *lambda[6], struct xyz &particles,
                    double m[][Pflat])
{
    // intput : lambda : the RBF coeff
    //          particles : all particles' coordinate
    //          panel_index;
    // output :  m: moments for panel_index^th panel
    
    
    double t1[P + 1];
    double t2[P + 1];
    double t3[P + 1];
    
    int i,j,k;
    
    for (i = 0; i < P + 1; i++) {
        t1[i] = tree[panel_index].t1[i];
        t2[i] = tree[panel_index].t2[i];
        t3[i] = tree[panel_index].t3[i];
    }
    
    double w1i[P + 1];
    double w2j[P + 1];
    double w3k[P + 1];
    double dj[P + 1];
    dj[0] = 0.5;
    dj[P] = 0.5;
    for (j = 1; j<P; j++)
        dj[j] = 1;
    
    for (j = 0; j < P + 1; j++)
        w3k[j] = w2j[j] = w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    
    int tp0 = tree[panel_index].members[0];
    int tp1 = tree[panel_index].members[1];
    int cluster_size;
    cluster_size = tp1 - tp0 + 1;
    
    size_t tp_j;
    
    double **a1i;
    double **a2j;
    double **a3k;
    double *D;
    
    a1i = (double**)calloc(P + 1, sizeof(double*));
    a2j = (double**)calloc(P + 1, sizeof(double*));
    a3k = (double**)calloc(P + 1, sizeof(double*));
    D = (double*)calloc(cluster_size + 1, sizeof(double));
    for (i = 0; i < (P + 1); i++)
    {
        a1i[i] = (double*)calloc(cluster_size + 1, sizeof(double));
        a2j[i] = (double*)calloc(cluster_size + 1, sizeof(double));
        a3k[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    }
    
    double x, y, z;
    double SumA1 = 0.0;
    double SumA2 = 0.0;
    double SumA3 = 0.0;
    double temp1,temp2,temp3;
    
    
    double temp_D;
    for (tp_j = tp0; tp_j <= tp1 ; tp_j++) {
        x = particles.x[tp_j];
        y = particles.y[tp_j];
        z = particles.z[tp_j];
        for (j = 0; j < P + 1; j++) {
            temp1 = w1i[j] / (x - t1[j]);
            temp2 = w2j[j] / (y - t2[j]);
            temp3 = w3k[j] / (z - t3[j]);
            a1i[j][tp_j - tp0] = temp1;
            a2j[j][tp_j - tp0] = temp2;
            a3k[j][tp_j - tp0] = temp3;
            
            SumA1 += temp1;
            SumA2 += temp2;
            SumA3 += temp3;
        }
        D[tp_j - tp0] = 1.0 / (SumA1 * SumA2 * SumA3);
        SumA1 = 0.0;
        SumA2 = 0.0;
        SumA3 = 0.0;
    }
    
    double s;
    int kk = -1;
    for (i = 0; i < P + 1; i++) {
        for (j = 0; j < P + 1; j++) {
            for (k = 0; k < P + 1; k++) {
                kk = kk + 1;
                for (tp_j = tp0; tp_j <= tp1; tp_j++)
                {
                    s = a1i[i][tp_j - tp0] * a2j[j][tp_j - tp0] * a3k[k][tp_j - tp0] * D[tp_j - tp0];
                    m[0][kk] += s * lambda[0][tp_j];
                    m[1][kk]  += s * lambda[1][tp_j];
                    m[2][kk]  += s * lambda[2][tp_j];
                    m[3][kk] += s * lambda[3][tp_j];
                    m[4][kk]  += s * lambda[4][tp_j];
                    m[5][kk] += s * lambda[5][tp_j];
                }
            }
        }
    }
    
    for (i = 0; i < P + 1; i++) {
        free(a1i[i]);
        free(a2j[i]);
        free(a3k[i]);
    }
    
    free(a1i);
    free(a2j);
    free(a3k);
    free(D);
    
}

//*****************************************************************************//
vec_3d Call_Treecode(double x, double y, double z, int panel_index)
{
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;
    velocity.val[3] = 0.0;
    velocity.val[4] = 0.0;
    velocity.val[5] = 0.0;
    
    // Far_Expan_Taylor(a_t, px, py, pz, t1,t2,t3);
    
    double temp_i;
    double temp_j;
    double temp_k;
    double R2, Rinv, Rinv3, Rinv5, Rinv7;
    
    double dx[P + 1];
    double dy[P + 1];
    double dz[P + 1];
    
    double s;
    double temp_moments[6];
    double G,Q,D;
    double DD,CC;
    
    for (int i = 0; i < P + 1; i++)
    {
        dx[i] = x - tree[panel_index].t1[i];
        dy[i] = y - tree[panel_index].t2[i];
        dz[i] = z - tree[panel_index].t3[i];
    }
    
    
    int kk = -1;
    double xx,yy,zz;
    
    for (int i = 0; i < P + 1; i++)
    {
        xx = dx[i];
        temp_i = xx * xx;
        for (int j = 0; j < P + 1; j++)
        {
            yy = dy[j];
            temp_j =  yy * yy;
            for (int k = 0; k < P + 1; k++)
            {
                zz = dz[k];
                temp_k = zz * zz;
                kk = kk + 1;
                
                temp_moments[0] = tree[panel_index].moments[0][kk];
                temp_moments[1] = tree[panel_index].moments[1][kk];
                temp_moments[2] = tree[panel_index].moments[2][kk];
                temp_moments[3] = tree[panel_index].moments[3][kk];
                temp_moments[4] = tree[panel_index].moments[4][kk];
                temp_moments[5] = tree[panel_index].moments[5][kk];
                
                
                R2 = temp_i + temp_j + temp_k;
                Rinv = 1/ sqrt(R2 + DEL * DEL);
                Rinv3 = Rinv * Rinv * Rinv;
                Rinv5 = Rinv * Rinv * Rinv3;
                Rinv7 = Rinv * Rinv * Rinv5 * 0.25;
                
                Q = (5 * DEL * DEL + 2 * R2) * Rinv5 * 0.5;
                D = 21 * DEL * DEL + 6 * R2;
                G = (2 * R2 + 3 * DEL * DEL) * Rinv3 - Rinv;
                
                velocity.val[0] += temp_moments[0] * G;
                velocity.val[1] += temp_moments[1] * G;
                velocity.val[2] += temp_moments[2] * G;
                
                DD = temp_moments[0] * xx + temp_moments[1] * yy + temp_moments[2] * zz;
                DD = DD * Rinv3;
                 velocity.val[0] +=  DD * xx;
                 velocity.val[1] +=  DD * yy;
                 velocity.val[2] +=  DD * zz;
                
                velocity.val[0] +=  (-temp_moments[4] * zz + temp_moments[5] * yy) * Q;
                velocity.val[1] +=  (temp_moments[3] * zz - temp_moments[5] * xx) * Q;
                velocity.val[2] +=  (-temp_moments[3] * yy + temp_moments[4] * xx) * Q;
                
                velocity.val[3] += (-temp_moments[1] * zz + temp_moments[2] * yy) * Q;
                velocity.val[4] += (temp_moments[0] * zz - temp_moments[2] * xx) * Q;
                velocity.val[5] += (-temp_moments[0] * yy + temp_moments[1] * xx) * Q;
                
                CC = (10 * DEL * DEL * DEL * DEL - 7 * DEL * DEL * R2 - 2 * R2 * R2) * Rinv7;
              
                velocity.val[3] += temp_moments[3] * CC;
                velocity.val[4] += temp_moments[4] * CC;
                velocity.val[5] += temp_moments[5] * CC;
                
                DD = temp_moments[3] * xx + temp_moments[4] * yy + temp_moments[5] * zz;
                DD = DD * Rinv7 * D;
                velocity.val[3] +=  DD * xx;
                velocity.val[4] +=  DD * yy;
                velocity.val[5] +=  DD * zz;
            }
        }
    }
    
    return velocity;
}

//*****************************************************************************//
vec_3d Call_Ds(int limit_1, int limit_2, int particle_index, double p_x, double p_y, double p_z, struct xyz &particles, double *lambda[3])
{
    
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    
    double ff[6];
    double x1,y1,z1;
    double R,R2,Rinv,Rinv3, Rinv5, Rinv7;
    
    
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;
    velocity.val[3] = 0.0;
    velocity.val[4] = 0.0;
    velocity.val[5] = 0.0;

    double G,Q,D;
    double DD,CC;
    for (size_t jj = limit_1; jj <= limit_2; jj++)
    {
        if (jj == particle_index)
            continue;
        
        x = p_x - particles.x[jj];
        y = p_y - particles.y[jj];
        z = p_z - particles.z[jj];
        ff[0] = lambda[0][jj];
        ff[1] = lambda[1][jj];
        ff[2] = lambda[2][jj];
        ff[3] = lambda[3][jj];
        ff[4] = lambda[4][jj];
        ff[5] = lambda[5][jj];
        x1 = x * x;
        y1 = y * y;
        z1 = z * z;
        
        R2 = x1 + y1 + z1;
        R = sqrt(R2 + DEL * DEL);
        Rinv = 1.0/R;
        Rinv3 = Rinv * Rinv * Rinv;
        Rinv5 = Rinv * Rinv * Rinv3;
        Rinv7 = Rinv * Rinv * Rinv5 * 0.25;
        
        Q = (5 * DEL * DEL + 2 * R2) * Rinv5 * 0.5;
        D = 21 * DEL * DEL + 6 * R2;
        G = (2 * R2 + 3 * DEL * DEL) * Rinv3 - Rinv;
        
        velocity.val[0] += ff[0] * G;
        velocity.val[1] += ff[1] * G;
        velocity.val[2] += ff[2] * G;
        
        DD = ff[0] * x + ff[1] * y + ff[2] * z;
        DD = DD * Rinv3;
        
        velocity.val[0] +=  DD * x;
        velocity.val[1] +=  DD * y;
        velocity.val[2] +=  DD * z;
        
        velocity.val[0] += (-ff[4] * z + ff[5] * y) * Q;
        velocity.val[1] += (ff[3] * z - ff[5] * x) * Q;
        velocity.val[2] += (-ff[3] * y + ff[4] * x) * Q;
        
        velocity.val[3] += (-ff[1] * z + ff[2] * y) * Q;
        velocity.val[4] += (ff[0] * z - ff[2] * x) * Q;
        velocity.val[5] += (-ff[0] * y + ff[1] * x) * Q;
        
        CC = (10 * DEL * DEL * DEL * DEL - 7 * DEL * DEL * R2 - 2 * R2 * R2) * Rinv7;
        
        velocity.val[3] += ff[3] * CC;
        velocity.val[4] += ff[4] * CC;
        velocity.val[5] += ff[5] * CC;
        
        DD = ff[3] * x + ff[4] * y + ff[5] * z;
        DD = DD * Rinv7 * D;
        
        velocity.val[3] +=  DD * x;
        velocity.val[4] +=  DD * y;
        velocity.val[5] +=  DD * z;
        
    }
     return velocity;

}

//*****************************************************************************//

vec_3d Comput_RBF(double *lambda[3], struct xyz &particles,size_t particle_index, size_t panel_index)
{
    // input :
	//         lambda : RBF coefficients
	//         particles : all particles coordinates
    //         particle_index
    //         panel_index
	// output : 
	//          velocity in 3D
	
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;
    velocity.val[3] = 0.0;
    velocity.val[4] = 0.0;
    velocity.val[5] = 0.0;
	
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
	{
        vec_3d tree_result = Call_Treecode(p_x, p_y, p_z, panel_index);
        velocity.val[0] = velocity.val[0] + tree_result.val[0];
        velocity.val[1] = velocity.val[1] + tree_result.val[1];
        velocity.val[2] = velocity.val[2] + tree_result.val[2];
        velocity.val[3] = velocity.val[3] + tree_result.val[3];
        velocity.val[4] = velocity.val[4] + tree_result.val[4];
        velocity.val[5] = velocity.val[5] + tree_result.val[5];
    }

    else
	{
		if (limit_2 - limit_1 < N0) //  otherwise, if cluster is a leaf, use direct sum
		{
            vec_3d DS_result = Call_Ds(limit_1, limit_2, particle_index,p_x,  p_y, p_z, particles, lambda);
            velocity.val[0] = velocity.val[0] + DS_result.val[0] ;
            velocity.val[1] = velocity.val[1] + DS_result.val[1] ;
            velocity.val[2] = velocity.val[2] + DS_result.val[2] ;
            velocity.val[3] = velocity.val[3] + DS_result.val[3] ;
            velocity.val[4] = velocity.val[4] + DS_result.val[4] ;
            velocity.val[5] = velocity.val[5] + DS_result.val[5] ;
		}//
		else // othervise, if cluster is not a leaf, look at children
		{
			velocity.val[0] = 0.0;
            velocity.val[1] = 0.0;
            velocity.val[2] = 0.0;
            velocity.val[3] = 0.0;
            velocity.val[4] = 0.0;
            velocity.val[5] = 0.0;
            
			size_t length = tree[panel_index].children.size();
			for (size_t i = 0; i < length; i++)
			{
				size_t index = tree[panel_index].children[i];
                vec_3d temp_result = Comput_RBF(lambda, particles, particle_index, index);
                velocity.val[0] = velocity.val[0] + temp_result.val[0];
                velocity.val[1] = velocity.val[1] + temp_result.val[1];
                velocity.val[2] = velocity.val[2] + temp_result.val[2];
                velocity.val[3] = velocity.val[3] + temp_result.val[3];
                velocity.val[4] = velocity.val[4] + temp_result.val[4];
                velocity.val[5] = velocity.val[5] + temp_result.val[5];
			}
            
		}
        
	}
    return velocity;
}

//*****************************************************************************//
void Cluster_Chev_Points(size_t tree_size)
{
    double h;
    h = 3.14159265358979323846/P;
    double t[P + 1] = {0.0};
    for (int i = 0; i < P + 1; i++)
        t[i] = cos(i * h);  //Chebyshev interpolation points [-1,1]
    
    double x1,x2,y1,y2,z1,z2;
    size_t tree_index;
    
    for (tree_index = 0; tree_index < tree_size ; tree_index++)
    {
        x1 = tree[tree_index].xinterval[0];
        x2 = tree[tree_index].xinterval[1];
        y1 = tree[tree_index].yinterval[0];
        y2 = tree[tree_index].yinterval[1];
        z1 = tree[tree_index].zinterval[0];
        z2 = tree[tree_index].zinterval[1];
        
        for (int i = 0; i < P + 1; i++) // map to the cluster
        {
            tree[tree_index].t1[i] =  x1 + (t[i] + 1)/2 * (x2 - x1);
            tree[tree_index].t2[i] =  y1 + (t[i] + 1)/2 * (y2 - y1);
            tree[tree_index].t3[i] =  z1 + (t[i] + 1)/2 * (z2 - z1);
        }
    }
}


//*****************************************************************************//
int main()
{
    tree.reserve(5000);
    leaf.reserve(5000);
	struct xyz particles(N_cube);
	cout << " ===== No box shrink ===========" << endl;
	cout << "P is " << P << endl;
	cout << "N_cube is " << N_cube << endl;
	cout << "theta is " << sqrt(sq_theta) << endl;
    cout << "N0 is " << N0 << endl;
	

	FILE * fp;
    char S_data_file[64] = {0};
    sprintf(S_data_file, "./rand_%d.txt", N);
    fp = fopen(S_data_file, "r");
    
	double x1, x2, x3, x4, x5, x6;
	int count = -1;
	if (fp == NULL)
	{
		cout << "Cannot open random points file" << endl;
		getchar();
		exit(0);
	}

	while (true)
	{
		count++;
		fscanf(fp,"%lf%lf%lf", &x1, &x2, &x3);
		if (feof(fp))
			break;
		if (count >= N_cube)
		{
			cout << "Out of range" << endl;
			exit(0);
		}

		particles.x[count] = x1;
		particles.y[count] = x2;
		particles.z[count] = x3;
		particles.index[count] = -1;
		particles.old_index[count] = count;
	}

    
    // *****************  get data *******************************
    /*srand(time(NULL));
	for (size_t count = 0; count < N_cube; count++)
	{
		particles.x[count] = L * ((double)rand() / (double)RAND_MAX) ; // [0,L]
		particles.y[count] = L * ((double)rand() / (double)RAND_MAX) ; // [0,L]
		particles.z[count] = L * ((double)rand() / (double)RAND_MAX) ; // [0,L]
		particles.index[count] = -1;
		particles.old_index[count] = count;
	}*/

	double *lambda[6];
	lambda[0] = new double[N_cube];
	lambda[1] = new double[N_cube];
	lambda[2] = new double[N_cube];
    lambda[3] = new double[N_cube];
    lambda[4] = new double[N_cube];
    lambda[5] = new double[N_cube];
    
    char lambda_Str_data_file[64] = {0};
    sprintf(lambda_Str_data_file, "./lambda_%d.txt", N);
    fp = fopen(lambda_Str_data_file, "r");
	
	count = -1;
	if (fp == NULL)
	{
		cout << "Cannot open lambda file" << endl;
		getchar();
		exit(0);
	}
    
	while (true)
	{
		count++;
		fscanf(fp,"%lf%lf%lf%lf%lf%lf", &x1, &x2, &x3 , &x4, &x5, &x6);
		if (feof(fp))
			break;
		if (count >= N_cube)
		{
			cout << "Out of range" << endl;
			exit(0);
		}
        
		lambda[0][count] = x1;
        lambda[1][count] = x2;
        lambda[2][count] = x3;
        lambda[3][count] = x4;
        lambda[4][count] = x5;
        lambda[5][count] = x6;
	}
    
    //***************** Set up tree *******************************
    long Start_total, Start_btree;
    long End_total, End_btree;
    
    Start_total = getTickCount(); // Get currenct CPU time
    Start_btree = getTickCount();
    
    build_tree_init();
    build_tree_3D_Recursive(0, particles, 0);
    
    End_btree = getTickCount();
    
   	//***************** Compute moment for each panel **************
	size_t size = tree.size();
    Cluster_Chev_Points(size);
    
	for (size_t i = 1; i < size; i++) // skip root
		Panel_Moment_B(i, lambda, particles, tree[i].moments);

	//***************** Compute Velocity ***************************
    
    vec_3d *velo = new vec_3d[N_cube];
    vec_3d *velo_old = new vec_3d[N_cube];
	
	for (int i = 0; i < N_cube; i++)
	{
		velo[i] = Comput_RBF(lambda, particles, i, 0);
	}
    End_total = getTickCount(); // Time for all treecode computing
    long treecode_cpu_time;
    treecode_cpu_time = End_total - Start_total;
    
    // output data to a file
    // P12_N40_theta0.2_N0200_time
    // P = 12, N = 40, theta = 0.2, N0 = 200, current time
    time_t raw_time;
    struct tm* time_info;
    char time_buffer[80];
    
    time(&raw_time);
    time_info = localtime(&raw_time);
    strftime(time_buffer, 80, "%Y-%m-%d-%H-%M-%S", time_info);
    int N_all = N_cube;
    double theta = sqrt(sq_theta);
    char file_name[256];
    sprintf(file_name, "BR_P%d_N%d_theta%.2f_N0%d_%s", P, N_all, theta, N0, time_buffer);
    ofstream output_file(file_name);
    
    output_file << " ===== No box shrink ===========" << endl;
    output_file << "P is " << P << endl;
    output_file << "N_cube is " << N_cube << endl;
    output_file << "theta is " << sqrt(sq_theta) << endl;
    output_file << "N0 is " << N0 << endl;
    
    cout << "N_cube is "<< N_cube << endl;
    
    cout << "treecode_cpu_time " << treecode_cpu_time << endl;
    cout << "build tree time is " << End_btree - Start_btree << endl;
    
    output_file << "treecode_cpu_time " << treecode_cpu_time << endl;
    output_file << "build tree time is " << End_btree - Start_btree << endl;
    
    //***************** Director summation and L_2 Error *****************
    vec_3d *v_true = new vec_3d[N_cube];
    
    //***************** compute v_true here ******************************
    
    long Start_ds, End_ds;
    Start_ds = getTickCount(); // Get currenct CPU time
    
    // new direct sum
    double temp_x, temp_y, temp_z;
    double px,py,pz;
    double xx,yy,zz;
    double xsq,ysq,zsq;
    double ff[6];
    double sum[6];
    double R,R2;
    double Rinv, Rinv3, Rinv5,Rinv7;
    double DD,CC;
    double G,Q,D;
    
    for (size_t i = 0; i < N_cube; i++)
    {
        temp_x = particles.x[i];
        temp_y = particles.y[i];
        temp_z = particles.z[i];
        for (int k = 0; k < 6; k++)
            sum[k] = 0.0;
        for (size_t j = 0; j < N_cube; j++)
        {
            px = particles.x[j];
            py = particles.y[j];
            pz = particles.z[j];
            xx = temp_x - px;
            yy = temp_y - py;
            zz = temp_z - pz;
            ff[0] = lambda[0][j];
            ff[1] = lambda[1][j];
            ff[2] = lambda[2][j];
            ff[3] = lambda[3][j];
            ff[4] = lambda[4][j];
            ff[5] = lambda[5][j];
            xsq = xx * xx;
            ysq = yy * yy;
            zsq = zz * zz;
            
            R2 = xsq + ysq + zsq ;
            R = sqrt(R2 + DEL * DEL);
            
            if (i != j)
            {
                Rinv = 1.0/R;
                Rinv3 = Rinv * Rinv * Rinv;
                Rinv5 = Rinv * Rinv * Rinv3;
                Rinv7 = Rinv * Rinv * Rinv5 * 0.25;
                
                Q = (5 * DEL * DEL + 2 * R2) * Rinv5 * 0.5;
                D = 21 * DEL * DEL + 6 * R2;
                G = (2 * R2 + 3 * DEL * DEL) * Rinv3 - Rinv;
                
                sum[0] += ff[0] * G ;
                sum[1] += ff[1] * G;
                sum[2] += ff[2] * G;
                
                DD = ff[0] * xx + ff[1] * yy + ff[2] * zz;
                DD = DD * Rinv3;
                sum[0] +=  DD * xx;
                sum[1] +=  DD * yy;
                sum[2] +=  DD * zz;
                
                sum[0] += (-ff[4] * zz + ff[5] * yy) * Q;
                sum[1] += (ff[3] * zz - ff[5] * xx) * Q;
                sum[2] += (-ff[3] * yy + ff[4] * xx) * Q;
                
                sum[3] += (-ff[1] * zz + ff[2] * yy) * Q;
                sum[4] += (ff[0] * zz - ff[2] * xx) * Q;
                sum[5] += (-ff[0] * yy + ff[1] * xx) * Q;
                
                CC = (10 * DEL * DEL * DEL * DEL - 7 * DEL * DEL * R2 - 2 * R2 * R2) * Rinv7;
                
                sum[3] += ff[3] * CC;
                sum[4] += ff[4] * CC;
                sum[5] += ff[5] * CC;
                
                DD = ff[3] * xx + ff[4] * yy + ff[5] * zz;
                DD = DD * Rinv7 * D;
                sum[3] +=  DD * xx;
                sum[4] +=  DD * yy;
                sum[5] +=  DD * zz;
            }
        }
        for (int k = 0; k < 6; k++)
            v_true[i].val[k] = sum[k];
    }
    
    End_ds = getTickCount(); // Get currenct CPU time
    long ds_cpu_time;
    ds_cpu_time = End_ds - Start_ds;
    
    cout << "ds time is " << ds_cpu_time << endl;
    
    output_file << "ds_cpu_time " << ds_cpu_time << endl;
	   
    //======= Err ===========================================================
    
    for (size_t i = 0; i < N_cube; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            //velo_old[particles.old_index[i]].val[0] = velo[i].val[0];
            //velo_old[particles.old_index[i]].val[1] = velo[i].val[1];
            //velo_old[particles.old_index[i]].val[2] = velo[i].val[2];
            velo_old[i].val[j] = velo[i].val[j];
        }
    }
    
    
    //====== L_infty Err======================================================
    
    double E = 0.0;
    double temp_d = 0.0;
    double temp_n = 0.0;
    double max_d = 0.0;
    double max_n = 0.0;
    
    for (size_t i = 0; i < N_cube; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            temp_n += (v_true[i].val[j] - velo_old[i].val[j]) * (v_true[i].val[j] - velo_old[i].val[j]) ;
            temp_d += velo_old[i].val[j] * velo_old[i].val[j];
        }
        temp_n = sqrt(temp_n);
        temp_d = sqrt(temp_d);
        if (temp_n > max_n)
            max_n = temp_n;
        if (temp_d > max_d)
            max_d = temp_d;
        temp_d = 0.0;
        temp_n = 0.0;
    }
    
    E = max_n/max_d;
    
    cout << "E is " << E << endl;
    cout << "max_n is " << max_n << endl;
    cout << "max_d is " << max_d << endl;
    
    output_file << "E is " << E << endl;
    
    
    //******** Error : my definition *************************
    double err2[6] = {0.0};
    double sum_d[6] = {0.0};
    double sum_n[6]= {0.0};
    
    for (int j = 0; j < 6 ; j++)
    {
        for (size_t i = 0; i < N_cube; i++)
        {
            sum_n[j] += (v_true[i].val[j] - velo_old[i].val[j]) * (v_true[i].val[j] - velo_old[i].val[j]);
            sum_d[j] += v_true[i].val[j] * v_true[i].val[j];
        }
    }
    
    for (int j = 0; j < 6; j++)
        err2[j] = sqrt(sum_n[j]/sum_d[j]);
    
    cout << "L2 err[0] is " << err2[0] << endl;
    cout << "L2 err[1] is " << err2[1] << endl;
    cout << "L2 err[2] is " << err2[2] << endl;
    cout << "Tree depth is " << max_level << endl;
    
    output_file << "L2 err[0] is " << err2[0] << endl;
    output_file << "L2 err[1] is " << err2[1] << endl;
    output_file << "L2 err[2] is " << err2[2] << endl;
    output_file << "Tree depth is " << max_level << endl;
    
    //******** Error : extend from RBF paper *******************
    double err2_ex = 0.0;
    double sum_d_ex = 0.0;
    double sum_n_ex = 0.0;
    for (size_t i = 0; i < N_cube; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            sum_n_ex += (v_true[i].val[j] - velo_old[i].val[j]) * (v_true[i].val[j] - velo_old[i].val[j]);
            sum_d_ex += v_true[i].val[j] * v_true[i].val[j];
        }
    }
    
    err2_ex = sqrt(sum_n_ex/sum_d_ex);
    cout<<"sum_d_ex = "<< sum_d_ex << endl;
    cout<<"sum_n_ex = "<< sum_n_ex << endl;
    
    cout << "E_2_ex is " << err2_ex << endl;
    
    output_file << "E_2_ex is " << err2_ex << endl;

    //*********************************************************

		
	delete [] lambda[0];
	delete [] lambda[1];
	delete [] lambda[2];
    delete [] lambda[3];
    delete [] lambda[4];
    delete [] lambda[5];
    
    delete [] velo_old;
    delete [] v_true;
    
    output_file.close();
	
	delete [] velo;
	
    cout << "done"<< endl;

	return 0;
}
