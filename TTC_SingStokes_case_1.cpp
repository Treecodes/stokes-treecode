
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>


using namespace std;

static const double L = 7.3681 ; // box size, L = 3.684, (N = 125k),L = 5.4288 (N = 400k), L = 7.3681 (N = 1000k)
static const int P = 0; // order of Taylor approximation
static const int Pflat = (P + 1)*(P + 2)*(P + 3)/6;
static const int N_cube = 1000000; // N points in one dimension, N_cube points total
static const int N0 = 2000;
static const double sq_theta = 0.04; // theta = 0.5
const bool UseSleep = false; // for testing memory usage purpose
int max_level = 0;

//**********************************************************//

struct vec_3d
{
    double val[3];
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
	double moments[3][Pflat];
	int moment_flag;
	panel() // initialization
	{
		moment_flag = 0;
		members[0] = 0;
		members[1] = -1;
		for(size_t index = 0; index < 3; index ++)
		{
			for (size_t kk = 0; kk < Pflat + 1; kk++)
                moments[index][kk] = 0;
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

	temp_panel.xinterval[0] = 0.0; // interval defining the panel
	temp_panel.xinterval[1] = L;
	temp_panel.yinterval[0] = 0.0;
	temp_panel.yinterval[1] = L;
	temp_panel.zinterval[0] = 0.0;
	temp_panel.zinterval[1] = L; // r = sqrt(3) * L / 2, r^2 = 3 * L*L/ 4 ;
	temp_panel.xc = 0.5 * L;
	temp_panel.yc = 0.5 * L;
	temp_panel.zc = 0.5 * L;
	temp_panel.MAC = (3 * L * L/ 4) / sq_theta; // MAC = r^2 / theta^2

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
// Calculate the far field expansions associated with each panel by recurrence relation
void Far_Expan_Taylor(double a[][P + 3][P + 3], double x, double y, double z,
					  size_t panel_index)
{
	// intput :
	//          x,y,z is the particle's three coordinates
	//          panel_index
	// output : a[][][] is the coefficients for panel_inde^th panel

	double xc = tree[panel_index].xc;   // coordernate of the center of panel
	double yc = tree[panel_index].yc;
	double zc = tree[panel_index].zc;

	double tp_x = x - xc;
	double tp_y = y - yc;
	double tp_z = z - zc;

    double R2 = tp_x * tp_x + tp_y * tp_y + tp_z * tp_z; // R2 = R^2
	double R = sqrt(R2);
	double s, ijk1, ijk2, sum_ijk;
	
	const int Q = P + 1;

	// base case
	a[1][1][1] = 1.0 / R; // first coefficient is the Greens function itself

	// two of the indeces are zero
	sum_ijk = 1;
	s = sum_ijk * R2;
	ijk2 = 2.0 * sum_ijk - 1.0;

	a[2][1][1] = (ijk2 * tp_x * a[1][1][1]) / s;
	a[1][2][1] = (ijk2 * tp_y * a[1][1][1]) / s;
	a[1][1][2] = (ijk2 * tp_z * a[1][1][1]) / s;
	if (P == 1)
		return;

	for (int i = 2; i < Q + 1; i++)
	{
		sum_ijk = i;
		s = sum_ijk * R2;
		ijk2 = 2.0 * sum_ijk - 1.0;
		ijk1 = sum_ijk - 1.0;

		a[i + 1][1][1] = (ijk2 * (tp_x * a[i][1][1]) - ijk1 * (a[i - 1][1][1])) / s;
		a[1][i + 1][1] = (ijk2 * (tp_y * a[1][i][1]) - ijk1 * (a[1][i - 1][1])) / s;
		a[1][1][i + 1] = (ijk2 * (tp_z * a[1][1][i]) - ijk1 * (a[1][1][i - 1])) / s;
	}

	// one index = 0; one index = 1; other index >= 1
	sum_ijk = 2;
	s = sum_ijk * R2;
	ijk2 = 2.0 * sum_ijk - 1.0;
	ijk1 = sum_ijk - 1.0;

	a[2][2][1] = (ijk2 * (tp_x * a[1][2][1] + tp_y * a[2][1][1])) / s ;
	a[2][1][2] = (ijk2 * (tp_x * a[1][1][2] + tp_z * a[2][1][1])) / s;
    a[1][2][2] = (ijk2 * (tp_y * a[1][1][2] + tp_z * a[1][2][1])) / s;

	for (int i = 2; i + 1 <= Q; i++)
	{
		sum_ijk = i + 1;
		s = sum_ijk * R2;
		ijk2 = 2.0 * sum_ijk - 1.0;
		ijk1 = sum_ijk - 1.0;

		a[i + 1][2][1] = (ijk2 * (tp_x * a[i][2][1] + tp_y * a[i + 1][1][1])
					  - ijk1 * (a[i - 1][2][1])) / s;

		a[i + 1][1][2] = (ijk2 * (tp_x * a[i][1][2] + tp_z * a[i + 1][1][1])
					  - ijk1 * (a[i - 1][1][2])) / s;

		a[2][i + 1][1] = (ijk2 * (tp_x * a[1][i + 1][1] + tp_y * a[2][i][1])
					  - ijk1 * (a[2][i - 1][1]))/s;

		a[1][i + 1][2] = (ijk2 * (tp_y * a[1][i][2] + tp_z * a[1][i + 1][1])
					  - ijk1 * (a[1][i - 1][2])) / s;

		a[2][1][i + 1] = (ijk2 * (tp_x * a[1][1][i + 1] + tp_z * a[2][1][i])
					  - ijk1 * (a[2][1][i - 1])) / s;

		a[1][2][i + 1] = (ijk2 * (tp_y * a[1][1][i + 1] + tp_z * a[1][2][i])
					  - ijk1 * (a[1][2][i - 1])) / s;
	}

	// one index = 0; other indices >= 2
	for (int i = 2; i <= Q - 2; i++)
	{
		for (int j = 2; i + j <= Q; j++)
		{
			sum_ijk = i + j;
			s = sum_ijk * R2;
			ijk2 = 2.0 * sum_ijk - 1.0;
			ijk1 = sum_ijk - 1.0;

			a[i + 1][j + 1][1] = (ijk2 * (tp_x * a[i][j + 1][1] + tp_y * a[i + 1][j][1]) - ijk1 * (a[i - 1][j + 1][1] + a[i + 1][j - 1][1])) / s;

			a[i + 1][1][j + 1] = (ijk2 * (tp_x * a[i][1][j + 1] + tp_z * a[i + 1][1][j]) - ijk1 * (a[i - 1][1][j + 1] + a[i + 1][1][j - 1])) / s;

			a[1][i + 1][j + 1] = (ijk2 * (tp_y * a[1][i][j + 1] + tp_z * a[1][i + 1][j]) - ijk1 * (a[1][i - 1][j + 1] + a[1][i + 1][j - 1])) / s;
			}
	}

	// two indices = 1, other index >= 1
	sum_ijk = 3;
	s = sum_ijk * R2;
	ijk2 = 2.0 * sum_ijk - 1.0;
	ijk1 = sum_ijk - 1.0;
	a[2][2][2] = (ijk2 * (tp_x * a[1][2][2] + tp_y * a[2][1][2] + tp_z * a[2][2][1])) / s;

	for (int i = 2; i + 2 <= Q; i++)
	{
		sum_ijk = i + 2;
		s = sum_ijk * R2;
		ijk2 = 2.0 * sum_ijk - 1.0;
		ijk1 = sum_ijk - 1.0;

		a[i + 1][2][2] = (ijk2 * (tp_x * a[i][2][2] + tp_y * a[i + 1][1][2] + tp_z * a[i + 1][2][1]) - ijk1 * (a[i - 1][2][2])) / s;

		a[2][i + 1][2] = (ijk2 * (tp_x * a[1][i + 1][2] + tp_y * a[2][i][2] + tp_z * a[2][i + 1][1]) - ijk1 * (a[2][i - 1][2])) / s;

		a[2][2][i + 1] = (ijk2 * (tp_x * a[1][2][i + 1] + tp_y * a[2][1][i + 1] + tp_z * a[2][2][i]) - ijk1 * (a[2][2][i - 1])) / s;
	}

	// one index = 1; other indeces >= 2
	for (int i = 2; i <= Q - 3; i++)
	{
		for (int j = 2; i + j < Q; j++)
		{
			sum_ijk = i + j + 1;
			s = sum_ijk * R2;
			ijk2 = 2.0 * sum_ijk - 1.0;
			ijk1 = sum_ijk - 1.0;

			a[i + 1][j + 1][2] = (ijk2 * (tp_x * a[i][j + 1][2] + tp_y * a[i + 1][j][2] + tp_z * a[i + 1][j + 1][1]) - ijk1 * (a[i - 1][j + 1][2] + a[i + 1][j - 1][2])) / s;
            a[i + 1][2][j + 1] = (ijk2 * (tp_x * a[i][2][j + 1] + tp_y * a[i + 1][1][j + 1] + tp_z * a[i + 1][2][j]) - ijk1 * (a[i - 1][2][j + 1] + a[i + 1][2][j - 1])) / s;
            a[2][i + 1][j + 1] = (ijk2 * (tp_x * a[1][i + 1][j + 1] + tp_y * a[2][i][j + 1] + tp_z * a[2][i + 1][j]) - ijk1 * (a[2][i - 1][j + 1] + a[2][i + 1][j - 1])) / s;
		}
	}

	// all indices >= 2
	for (int i = 2; i <= Q - 4; i++)
	{
		for (int j = 2; i + j <= Q - 2; j++)
		{
			for (int k = 2; i + j + k <= Q; k++)
			{
				sum_ijk = i + j + k;
				s = sum_ijk * R2;
				ijk2 = 2.0 * sum_ijk - 1.0;
				ijk1 = sum_ijk - 1.0;

				a[i + 1][j + 1][k + 1] = (ijk2 * (tp_x * a[i][j + 1][k + 1] + tp_y * a[i + 1][j][k + 1] + tp_z * a[i + 1][j + 1][k])
				- ijk1 * (a[i - 1][j + 1][k + 1] + a[i + 1][j - 1][k + 1] + a[i + 1][j + 1][k - 1])) / s;
			}
		}
	}
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
void Panel_Moment_Taylor(size_t panel_index, double *lambda[3], struct xyz &particles,
						 double m[][Pflat])
{
    // intput : lambda : the RBF coeff
    //          particles : all particles' coordinate
    //          panel_index;
    // output :  m: moments for panel_index^th panel
    double xc = tree[panel_index].xc;
    double yc = tree[panel_index].yc;
    double zc = tree[panel_index].zc;
    
    double tp0 = tree[panel_index].members[0];
    double tp1 = tree[panel_index].members[1];
    
    double x, y, z;
    double s;
    double sum = 0.0;
    size_t tp_j;
    
    int kk = 0;
 
    for (int k1 = 0; k1 < P + 1; k1++)
    {
        for (int k2 = 0; k1 + k2 < P + 1; k2++)
        {
            for (int k3 = 0; k1 + k2 + k3 < P + 1; k3++)
            {
                kk = kk + 1;
                for (tp_j = tp0; tp_j <= tp1; tp_j++)
                {
                    x = particles.x[tp_j];
                    y = particles.y[tp_j];
                    z = particles.z[tp_j];
                    s = mypow(x - xc, k1) * mypow(y - yc,k2) * mypow(z - zc, k3);
                    
                    m[0][kk] += s*lambda[0][tp_j];
                    m[1][kk] += s*lambda[1][tp_j];
                    m[2][kk] += s*lambda[2][tp_j];
                }
            }
        }
    }
}
//*****************************************************************************//
vec_3d Call_Treecode(double px, double py, double pz, int panel_index)
{
    double a_t[P + 3][P + 3][P + 3] = {{{0.0}}};
    
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;
    
    Far_Expan_Taylor(a_t, px, py, pz, panel_index);
    
    double xc = tree[panel_index].xc;   // coordernate of the center of panel
    double yc = tree[panel_index].yc;
    double zc = tree[panel_index].zc;
    
    double tp_x = px - xc;
    double tp_y = py - yc;
    double tp_z = pz - zc;
    double s;
    double moments[3];
    double term2;
    int k1p1,k2p1,k3p1;
    
    int kk = 0;
    for (int k1 = 0; k1 < P + 1; k1++)
    {
        for (int k2 = 0; k1 + k2 < P + 1; k2++)
        {
            for (int k3 = 0; k1 + k2 + k3 < P + 1; k3++)
            {
                kk = kk + 1;
                moments[0] = tree[panel_index].moments[0][kk];
                moments[1] = tree[panel_index].moments[1][kk];
                moments[2] = tree[panel_index].moments[2][kk];
                
                term2 = tp_x * moments[0] + tp_y * moments[1] + tp_z * moments[2];
                
                k1p1 = k1 + 1;
                k2p1 = k2 + 1;
                k3p1 = k3 + 1;
                
                s = a_t[k1p1][k2p1][k3p1];
                
                velocity.val[0] += (1 - k1) * s * moments[0] + k1p1 * (a_t[k1 + 2][k2p1][k3p1] * term2 - a_t[k1 + 2][k2][k3p1] * moments[1] - a_t[k1 + 2][k2p1][k3] * moments[2]);
                
                velocity.val[1] += (1 - k2) * s * moments[1] + k2p1 * (a_t[k1p1][k2 + 2][k3p1] * term2 - a_t[k1][k2 + 2][k3p1] * moments[0] - a_t[k1p1][k2 + 2][k3] * moments[2]);
                
                velocity.val[2] += (1 - k3) * s * moments[2] + k3p1 * (a_t[k1p1][k2p1][k3 + 2] * term2 - a_t[k1][k2p1][k3 + 2] * moments[0] - a_t[k1p1][k2][k3 + 2] * moments[1]);
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
    
    double ff[3];
    double R, Rinv, Rinv3;
    double DD;
    
    vec_3d velocity;
    velocity.val[0] = 0.0;
    velocity.val[1] = 0.0;
    velocity.val[2] = 0.0;

    
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
        R = sqrt(x*x + y*y + z*z);
        Rinv = 1.0/R;
        Rinv3 = Rinv * Rinv * Rinv;
        velocity.val[0] += ff[0] * Rinv;
        velocity.val[1] += ff[1] * Rinv;
        velocity.val[2] += ff[2] * Rinv;
        
        DD = ff[0] * x + ff[1] * y + ff[2] * z;
        DD = DD * Rinv3;
        velocity.val[0] += DD * x;
        velocity.val[1] += DD * y;
        velocity.val[2] += DD * z;
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
    }

    else
	{
		if (limit_2 - limit_1 < N0) //  otherwise, if cluster is a leaf, use direct sum
		{
            vec_3d DS_result = Call_Ds(limit_1, limit_2, particle_index,p_x,  p_y, p_z, particles, lambda);
            velocity.val[0] = velocity.val[0] + DS_result.val[0] ;
            velocity.val[1] = velocity.val[1] + DS_result.val[1] ;
            velocity.val[2] = velocity.val[2] + DS_result.val[2] ;
		}//
		else // othervise, if cluster is not a leaf, look at children
		{
			velocity.val[0] = 0.0;
            velocity.val[1] = 0.0;
            velocity.val[2] = 0.0;
			size_t length = tree[panel_index].children.size();
			for (size_t i = 0; i < length; i++)
			{
				size_t index = tree[panel_index].children[i];
                vec_3d temp_result = Comput_RBF(lambda, particles, particle_index, index);
                velocity.val[0] = velocity.val[0] + temp_result.val[0];
                velocity.val[1] = velocity.val[1] + temp_result.val[1];
                velocity.val[2] = velocity.val[2] + temp_result.val[2];
			}
            
		}
        
	}
    return velocity;
}

//*****************************************************************************//
int main()
{
	struct xyz particles(N_cube);
	cout << " ===== No box shrink ===========" << endl;
	cout << "P is " << P << endl;
	cout << "N_cube is " << N_cube << endl;
	cout << "theta is " << sqrt(sq_theta) << endl;
    cout << "N0 is " << N0 << endl;
	

	FILE * fp;
    char S_data_file[64] = {0};
    sprintf(S_data_file, "./data/rand_%d.txt", N_cube);
    fp = fopen(S_data_file, "r");
    
   // fp = fopen("rand_125000.txt", "r");
	double x1, x2, x3;
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

	double *lambda[3];
	lambda[0] = new double[N_cube];
	lambda[1] = new double[N_cube];
	lambda[2] = new double[N_cube];
    
    char lambda_Str_data_file[64] = {0};
    sprintf(lambda_Str_data_file, "./data/lambda_%d.txt", N_cube);
    fp = fopen(lambda_Str_data_file, "r");
    
   // fp = fopen("lambda_125000.txt", "r");
	
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
		fscanf(fp,"%lf%lf%lf", &x1, &x2, &x3);
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
	for (size_t i = 1; i < size; i++) // skip root
		Panel_Moment_Taylor(i, lambda, particles, tree[i].moments);

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
    double ff[3];
    double sum[3];
    double R;
    double Rinv, Rinv3;
    double DD;
    for (size_t i = 0; i < N_cube; i++)
    {
        temp_x = particles.x[i];
        temp_y = particles.y[i];
        temp_z = particles.z[i];
        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;
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
            R = sqrt(xx*xx + yy*yy + zz*zz);
            if (i != j)
            {
                Rinv = 1.0/R;
                sum[0] += ff[0] * Rinv;
                sum[1] += ff[1] * Rinv;
                sum[2] += ff[2] * Rinv;
                Rinv3 = Rinv * Rinv * Rinv;
                DD = ff[0] * xx + ff[1] * yy + ff[2] * zz;
                DD = DD * Rinv3;
                sum[0] = sum[0] + DD * xx;
                sum[1] = sum[1] + DD * yy;
                sum[2] = sum[2] + DD * zz;
            }
        }
        v_true[i].val[0] = sum[0];
        v_true[i].val[1] = sum[1];
        v_true[i].val[2] = sum[2];
    }
    
    End_ds = getTickCount(); // Get currenct CPU time
    long ds_cpu_time;
    ds_cpu_time = End_ds - Start_ds;
    
    cout << "ds time is " << ds_cpu_time << endl;
    
    output_file << "ds_cpu_time " << ds_cpu_time << endl;
	   
    //======= Err ===========================================================
    
    for (size_t i = 0; i < N_cube; i++)
    {
        //velo_old[particles.old_index[i]].val[0] = velo[i].val[0];
        //velo_old[particles.old_index[i]].val[1] = velo[i].val[1];
        //velo_old[particles.old_index[i]].val[2] = velo[i].val[2];
        velo_old[i].val[0] = velo[i].val[0];
        velo_old[i].val[1] = velo[i].val[1];
        velo_old[i].val[2] = velo[i].val[2];
    }
    
    
    //====== L_infty Err======================================================
    
    double E = 0.0;
    double temp_d = 0.0;
    double temp_n = 0.0;
    double max_d = 0.0;
    double max_n = 0.0;
    
    for (size_t i = 0; i < N_cube; i++)
    {
        for (int j = 0; j < 3; j++)
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
    double err2[3] = {0.0};
    double sum_d[3] = {0.0};
    double sum_n [3]= {0.0};
    
    for (int j = 0; j < 3; j++)
    {
        for (size_t i = 0; i < N_cube; i++)
        {
            sum_n[j] += (v_true[i].val[j] - velo_old[i].val[j]) * (v_true[i].val[j] - velo_old[i].val[j]);
            sum_d[j] += v_true[i].val[j] * v_true[i].val[j];
        }
    }
    
    for (int j = 0; j < 3; j++)
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
        sum_n_ex += (v_true[i].val[0] - velo_old[i].val[0]) * (v_true[i].val[0] - velo_old[i].val[0])
        + (v_true[i].val[1] - velo_old[i].val[1]) * (v_true[i].val[1] - velo_old[i].val[1])
        + (v_true[i].val[2] - velo_old[i].val[2]) * (v_true[i].val[2] - velo_old[i].val[2]);
        
        sum_d_ex += v_true[i].val[0] * v_true[i].val[0] + v_true[i].val[1] * v_true[i].val[1] +
        v_true[i].val[2] * v_true[i].val[2];
    }
    
    err2_ex = sqrt(sum_n_ex/sum_d_ex);
    
    cout << "E_2_ex is " << err2_ex << endl;
    
    output_file << "E_2_ex is " << err2_ex << endl;

    //*********************************************************

		
	delete [] lambda[0];
	delete [] lambda[1];
	delete [] lambda[2];
    delete [] velo_old;
    delete [] v_true;
    
    output_file.close();
	
	delete [] velo;
	
    cout << "done"<< endl;

	return 0;
}
