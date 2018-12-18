
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>

using namespace std;

static const int N = 81920; // N points in each direction
static const int P = 8; // order of Taylor approximation
static const int Pflat = (P + 1)*(P + 2)*(P + 3)/6;
static const size_t N_cube = N ; // N_cube points total
static const int N0 = 2000;
static const double sq_theta = 0.64; // theta = 0.1
const bool UseSleep = false; // for testing memory usage purpose
int max_level = 0;

//**********************************************************//

struct vec_3d
{
    double val_Str[3];
    double val_Sto[3];
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
	double moments_Sto[3][Pflat];
    double moments_Str[3][3][Pflat];
	int moment_flag;
	panel() // initialization
	{
		moment_flag = 0;
		members[0] = 0;
		members[1] = -1;
		for(size_t index = 0; index < 3; index ++)
		{
            for (size_t kk = 0; kk < Pflat + 1; kk++)
                moments_Sto[index][kk] = 0.0;
		}
        for(size_t index1 = 0; index1 < 3; index1 ++)
		{
            for(size_t index2 = 0; index2 < 3; index2 ++)
            {
                for (size_t kk = 0; kk < Pflat + 1; kk++)
                moments_Str[index1][index2][kk] = 0.0;
			}
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

double delta(int i, int j)
{
    double d = 0.0;
    if (fabs(i - j) < 1e-5)
        d = 1.0;
    return d;
}

//*****************************************************************************//
void min_max(struct xyz &s, double min[], double max[])
{
    for(int i = 0; i< N_cube; i++)
    {
        if (s.x[i] <  min[0])
            min[0] = s.x[i];
        if (s.y[i] <  min[1])
            min[1] = s.y[i];
        if (s.z[i] <  min[2])
            min[2] = s.z[i];
        if (s.x[i] >  max[0])
            max[0] = s.x[i];
        if (s.y[i] >  max[1])
            max[1] = s.y[i];
        if (s.z[i] >  max[2])
            max[2] = s.z[i];
    }
}

//*****************************************************************************//
void build_tree_init(struct xyz &s)
{
	panel temp_panel;
    double min[3] = {10};
    double max[3] = {-10};
    
    min_max(s, min,max);

	// indices of particles belonging to panel
	temp_panel.members[0] = 0;
	temp_panel.members[1] = N_cube - 1;
    double r2 = (max[0] - min[0])*(max[0] - min[0]) + (max[1] - min[1])*(max[1] - min[1]) + (max[2] - min[2])*(max[2] - min[2]);

	temp_panel.xinterval[0] = min[0]; // interval defining the panel
	temp_panel.xinterval[1] = max[0];
	temp_panel.yinterval[0] = min[1];
	temp_panel.yinterval[1] = max[1];
	temp_panel.zinterval[0] = min[2];
	temp_panel.zinterval[1] = max[2]; // r = sqrt(3) / 2, r^2 = 3 / 4 = 0.75;
	temp_panel.xc = 1;
	temp_panel.yc = 0;
	temp_panel.zc = -2;
	temp_panel.MAC = r2 / sq_theta; // MAC = r^2 / theta^2

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
void ShrinkCluster(struct xyz &particles)
{
    double min_x = 0.0;;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    double min_z = 0.0;
    double max_z = 0.0;
    double r;
    size_t size = tree.size();
    int limit_1;
    int limit_2;
	for (size_t i = 0; i < size; i++)
    {
        
		limit_1 = tree[i].members[0];
		limit_2 = tree[i].members[1];
        
        min_x = particles.x[limit_1]; // the first point in the current cluster
        max_x = particles.x[limit_1];
        min_y = particles.y[limit_1];
        max_y = particles.y[limit_1];
        min_z = particles.z[limit_1];
        max_z = particles.z[limit_1];
        
        for (size_t jj = limit_1; jj <= limit_2; jj++)
        {
            if (particles.x[jj] < min_x)
                min_x = particles.x[jj];
            if (particles.x[jj] > max_x)
                max_x = particles.x[jj];
            if (particles.y[jj] < min_y)
                min_y = particles.y[jj];
            if (particles.y[jj] > max_y)
                max_y = particles.y[jj];
            if (particles.z[jj] < min_z)
                min_z = particles.z[jj];
            if (particles.z[jj] > max_z)
                max_z = particles.z[jj];
        }
        tree[i].xc = (min_x + max_x)/2.0;
        tree[i].yc = (min_y + max_y)/2.0;
        tree[i].zc = (min_z + max_z)/2.0;
        r = sqrt ((max_x - min_x) * (max_x - min_x) + (max_y - min_y) * (max_y - min_y) + (max_z - min_z) * (max_z - min_z))/2;
        tree[i].MAC = r*r/ sq_theta;
    }
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
void Far_Expan_Taylor(double a[][P + 4][P + 4], double x, double y, double z,
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
	
	const int Q = P + 2;

	// base case
	a[0 + 1][0 + 1][0 + 1] = 1.0 / R; // first coefficient is the Greens function itself

	// two of the indeces are zero
	sum_ijk = 1;
	s = sum_ijk * R2;
	ijk2 = 2.0 * sum_ijk - 1.0;

	a[1 + 1][0 + 1][0 + 1] = (ijk2 * tp_x * a[0 + 1][0 + 1][0 + 1]) / s;
	a[0 + 1][1 + 1][0 + 1] = (ijk2 * tp_y * a[0 + 1][0 + 1][0 + 1]) / s;
	a[0 + 1][0 + 1][1 + 1] = (ijk2 * tp_z * a[0 + 1][0 + 1][0 + 1]) / s;
	if (P == 1)
		return;

	for (int i = 2; i < Q + 1; i++)
	{
		sum_ijk = i;
		s = sum_ijk * R2;
		ijk2 = 2.0 * sum_ijk - 1.0;
		ijk1 = sum_ijk - 1.0;

		a[i + 1][0 + 1][0 + 1] = (ijk2 * (tp_x * a[i - 1 + 1][0 + 1][0 + 1]) - ijk1 * (a[i - 2 + 1][0 + 1][0 + 1])) / s;
		a[0 + 1][i + 1][0 + 1] = (ijk2 * (tp_y * a[0 + 1][i - 1 + 1][0 + 1]) - ijk1 * (a[0 + 1][i - 2 + 1][0 + 1])) / s;
		a[0 + 1][0 + 1][i + 1] = (ijk2 * (tp_z * a[0 + 1][0 + 1][i - 1 + 1]) - ijk1 * (a[0 + 1][0 + 1][i - 2 + 1])) / s;
	}

	// one index = 0; one index = 1; other index >= 1
	sum_ijk = 2;
	s = sum_ijk * R2;
	ijk2 = 2.0 * sum_ijk - 1.0;
	ijk1 = sum_ijk - 1.0;

	a[1 + 1][1 + 1][0 + 1] = (ijk2 * (tp_x * a[0 + 1][1 + 1][0 + 1] + tp_y * a[1 + 1][0 + 1][0 + 1])) / s ;
	a[1 + 1][0 + 1][1 + 1] = (ijk2 * (tp_x * a[0 + 1][0 + 1][1 + 1] + tp_z * a[1 + 1][0 + 1][0 + 1])) / s;
    a[0 + 1][1 + 1][1 + 1] = (ijk2 * (tp_y * a[0 + 1][0 + 1][1 + 1] + tp_z * a[0 + 1][1 + 1][0 + 1])) / s;

	for (int i = 2; i + 1 <= Q; i++)
	{
		sum_ijk = i + 1;
		s = sum_ijk * R2;
		ijk2 = 2.0 * sum_ijk - 1.0;
		ijk1 = sum_ijk - 1.0;

		a[i + 1][1 + 1][0 + 1] = (ijk2 * (tp_x * a[i - 1 + 1][1 + 1][0 + 1] + tp_y * a[i + 1][0 + 1][0 + 1])
					  - ijk1 * (a[i - 2 + 1][1 + 1][0 + 1])) / s;

		a[i + 1][0 + 1][1 + 1] = (ijk2 * (tp_x * a[i - 1 + 1][0 + 1][1 + 1] + tp_z * a[i + 1][0 + 1][0 + 1])
					  - ijk1 * (a[i - 2 + 1][0 + 1][1 + 1])) / s;

		a[1 + 1][i + 1][0 + 1] = (ijk2 * (tp_x * a[0 + 1][i + 1][0 + 1] + tp_y * a[1 + 1][i - 1 + 1][0 + 1])
					  - ijk1 * (a[1 + 1][i - 2 + 1][0 + 1]))/s;

		a[0 + 1][i + 1][1 + 1] = (ijk2 * (tp_y * a[0 + 1][i - 1 + 1][1 + 1] + tp_z * a[0 + 1][i + 1][0 + 1])
					  - ijk1 * (a[0 + 1][i - 2 + 1][1 + 1])) / s;

		a[1 + 1][0 + 1][i + 1] = (ijk2 * (tp_x * a[0 + 1][0 + 1][i + 1] + tp_z * a[1 + 1][0 + 1][i - 1 + 1])
					  - ijk1 * (a[1 + 1][0 + 1][i - 2 + 1])) / s;

		a[0 + 1][1 + 1][i + 1] = (ijk2 * (tp_y * a[0 + 1][0 + 1][i + 1] + tp_z * a[0 + 1][1 + 1][i - 1 + 1])
					  - ijk1 * (a[0 + 1][1 + 1][i - 2 + 1])) / s;
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

			a[i + 1][j + 1][0 + 1] = (ijk2 * (tp_x * a[i - 1 + 1][j + 1][0 + 1] + tp_y * a[i + 1][j - 1 + 1][0 + 1])
						  - ijk1 * (a[i - 2 + 1][j + 1][0 + 1] + a[i + 1][j - 2 + 1][0 + 1])) / s;

			a[i + 1][0 + 1][j + 1] = (ijk2 * (tp_x * a[i - 1 + 1][0 + 1][j + 1] + tp_z * a[i + 1][0 + 1][j - 1 + 1])
						  - ijk1 * (a[i - 2 + 1][0 + 1][j + 1] + a[i + 1][0 + 1][j - 2 + 1])) / s;

			a[0 + 1][i + 1][j + 1] = (ijk2 * (tp_y * a[0 + 1][i - 1 + 1][j + 1] + tp_z * a[0 + 1][i + 1][j - 1 + 1])
						  - ijk1 * (a[0 + 1][i - 2 + 1][j + 1] + a[0 + 1][i + 1][j - 2 + 1])) / s;
			}
	}

	// two indices = 1, other index >= 1
	sum_ijk = 3;
	s = sum_ijk * R2;
	ijk2 = 2.0 * sum_ijk - 1.0;
	ijk1 = sum_ijk - 1.0;
	a[1 + 1][1 + 1][1 + 1] = (ijk2 * (tp_x * a[0 + 1][1 + 1][1 + 1] + tp_y * a[1 + 1][0 + 1][1 + 1] + tp_z * a[1 + 1][1 + 1][0 + 1])) / s;

	for (int i = 2; i + 2 <= Q; i++)
	{
		sum_ijk = i + 2;
		s = sum_ijk * R2;
		ijk2 = 2.0 * sum_ijk - 1.0;
		ijk1 = sum_ijk - 1.0;

		a[i + 1][1 + 1][1 + 1] = (ijk2 * (tp_x * a[i - 1 + 1][1 + 1][1 + 1] + tp_y * a[i + 1][0 + 1][1 + 1] + tp_z * a[i + 1][1 + 1][0 + 1])
					  - ijk1 * (a[i - 2 + 1][1 + 1][1 + 1])) / s;

		a[1 + 1][i + 1][1 + 1] = (ijk2 * (tp_x * a[0 + 1][i + 1][1 + 1] + tp_y * a[1 + 1][i - 1 + 1][1 + 1] + tp_z * a[1 + 1][i + 1][0 + 1])
					  - ijk1 * (a[1 + 1][i - 2 + 1][1 + 1])) / s;

		a[1 + 1][1 + 1][i + 1] = (ijk2 * (tp_x * a[0 + 1][1 + 1][i + 1] + tp_y * a[1 + 1][0 + 1][i + 1] + tp_z * a[1 + 1][1 + 1][i - 1 + 1])
					  - ijk1 * (a[1 + 1][1 + 1][i - 2 + 1])) / s;
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

			a[i + 1][j + 1][1 + 1] = (ijk2 * (tp_x * a[i - 1 + 1][j + 1][1 + 1] + tp_y * a[i + 1][j - 1 + 1][1 + 1] + tp_z * a[i + 1][j + 1][0 + 1])
						  - ijk1 * (a[i - 2 + 1][j + 1][1 + 1] + a[i + 1][j - 2 + 1][1 + 1])) / s;
            a[i + 1][1 + 1][j + 1] = (ijk2 * (tp_x * a[i - 1 + 1][1 + 1][j + 1] + tp_y * a[i + 1][0 + 1][j + 1] + tp_z * a[i + 1][1 + 1][j - 1 + 1])
						  - ijk1 * (a[i - 2 + 1][1 + 1][j + 1] + a[i + 1][1 + 1][j - 2 + 1])) / s;
            a[1 + 1][i + 1][j + 1] = (ijk2 * (tp_x * a[0 + 1][i + 1][j + 1] + tp_y * a[1 + 1][i - 1 + 1][j + 1] + tp_z * a[1 + 1][i + 1][j - 1 + 1])
						  - ijk1 * (a[1 + 1][i - 2 + 1][j + 1] + a[1 + 1][i + 1][j - 2 + 1])) / s;
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

				a[i + 1][j + 1][k + 1] = (ijk2 * (tp_x * a[i - 1 + 1][j + 1][k + 1] + tp_y * a[i + 1][j - 1 + 1][k + 1] + tp_z * a[i + 1][j + 1][k - 1 + 1])
				- ijk1 * (a[i - 2 + 1][j + 1][k + 1] + a[i + 1][j - 2 + 1][k + 1] + a[i + 1][j + 1][k - 2 + 1])) / s;
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
void Panel_Moment_Taylor(size_t panel_index,double *lambda_Sto[3], double *lambda_Str[3], double *normal[3],struct xyz &particles, double m_Sto[][Pflat], double m_Str[][3][Pflat])
{
	// intput : lambda : the RBF coeff
	//          particles : all particles' coordinate
	//          panel_index;
	// output :  m: moments for panel_index^th panel
	double xc = tree[panel_index].xc;
	double yc = tree[panel_index].yc;
	double zc = tree[panel_index].zc;

    double x, y, z;
	size_t tp_j;
    double tp0 = tree[panel_index].members[0];
    double tp1 = tree[panel_index].members[1];

    int kk = -1;

    double s;
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
                    
                    m_Sto[0][kk] += s * lambda_Sto[0][tp_j];
                    m_Sto[1][kk] += s * lambda_Sto[1][tp_j];
                    m_Sto[2][kk] += s * lambda_Sto[2][tp_j];
                    
                    for (int index1 = 0; index1 < 3; index1++)
                    {
                        for (int index2 = 0; index2 < 3; index2++)
                        {
                            m_Str[index1][index2][kk] += s * lambda_Str[index1][tp_j] * normal[index2][tp_j];
                        }
                    }
					
				}
			}
		}
	}
}

//*****************************************************************************//
vec_3d Comput_RBF(double *lambda_Sto[3], double *lambda_Str[3], double *normal[3],  struct xyz &particles, size_t particle_index, size_t panel_index, int dep = 0 )
{
    // input :
	//         lambda : RBF coefficients
	//         particles : all particles coordinates
    //         particle_index
    //         panel_index
	// output : 
	//          velocity in 3D
	
	double a_t[P + 4][P + 4][P + 4] = {{{0.0}}};
    double C[3][3][3] = {{{0.0}}};
    vec_3d velocity;
    vec_3d temp_result;
    
    for (int index1 = 0; index1 < 3; index1++)
    {
        velocity.val_Str[index1] = 0.0;
        temp_result.val_Str[index1] = 0.0;
        velocity.val_Sto[index1] = 0.0;
        temp_result.val_Sto[index1] = 0.0;
    }
	
	
	double dx = 0.0;
 	double dy = 0.0;
    double dz = 0.0;
    double p_x = 0.0;
    double p_y = 0.0;
    double p_z = 0.0;
	double xc = 0.0;
	double yc = 0.0;
    double zc = 0.0;
	
	double R_sq = 0.0;
    double R;
    double ff_sto[3];
    double ff_str[3];
    double rnorm[3];
    double Rinv, Rinv3, Rinv5;
    double DD_sto;
    double DD_str;
    
	size_t limit_1;
	size_t limit_2;
	
	p_x = particles.x[particle_index];
	p_y = particles.y[particle_index];
    p_z = particles.z[particle_index];
	
	xc = tree[panel_index].xc;
	yc = tree[panel_index].yc;
    zc = tree[panel_index].zc;
	
    double tp_x,tp_y,tp_z;
    tp_x = p_x - xc;
    tp_y = p_y - yc;
    tp_z = p_z - zc;
	
	R_sq = tp_x * tp_x + tp_y * tp_y + tp_z * tp_z;
    double term2;
    int k1p1,k2p1,k3p1;
    double moments_sto[3];
    double nv1,nv2,nv3;
    double MM;
    double M1p2,M1p3,M2p3;
    double moments_str[3][3];
    
    
    if (tree[panel_index].MAC < R_sq)
	{
		// if MAC is satisfied, use Taylor approximation
    
        Far_Expan_Taylor(a_t, p_x, p_y, p_z, panel_index);
        double s;
  
        int kk = -1;
        
        for (int k1 = 0; k1 < P + 1; k1++)
        {
            for (int k2 = 0; k1 + k2 < P + 1; k2++)
            {
                for (int k3 = 0; k1 + k2 + k3 < P + 1; k3++)
                {
                    kk = kk + 1;
                    moments_str[0][0] = tree[panel_index].moments_Str[0][0][kk];
                    moments_str[0][1] = tree[panel_index].moments_Str[0][1][kk];
                    moments_str[0][2] = tree[panel_index].moments_Str[0][2][kk];
                    
                    moments_str[1][0] = tree[panel_index].moments_Str[1][0][kk];
                    moments_str[1][1] = tree[panel_index].moments_Str[1][1][kk];
                    moments_str[1][2] = tree[panel_index].moments_Str[1][2][kk];
                    
                    moments_str[2][0] = tree[panel_index].moments_Str[2][0][kk];
                    moments_str[2][1] = tree[panel_index].moments_Str[2][1][kk];
                    moments_str[2][2] = tree[panel_index].moments_Str[2][2][kk];
                    
                    nv1 = tp_x * moments_str[0][0] + tp_y * moments_str[0][1] + tp_z * moments_str[0][2];
                    nv2 = tp_x * moments_str[1][0] + tp_y * moments_str[1][1] + tp_z * moments_str[1][2];
                    nv3 = tp_x * moments_str[2][0] + tp_y * moments_str[2][1] + tp_z * moments_str[2][2];
                    
                    MM = moments_str[0][0] + moments_str[1][1] + moments_str[2][2];
                    
                    M1p2 = moments_str[0][1] + moments_str[1][0];
                    M1p3 = moments_str[0][2] + moments_str[2][0];
                    M2p3 = moments_str[1][2] + moments_str[2][1];
                    
                    k1p1 = k1 + 1;
                    k2p1 = k2 + 1;
                    k3p1 = k3 + 1;
                    
                    velocity.val_Str[0] += (k1p1 * ((k1 + 2) * a_t[k1p1 + 2][k2p1][k3p1] * nv1 +
                                                       k2p1 * a_t[k1p1 + 1][k2p1 + 1][k3p1] * nv2 +
                                                       k3p1 * a_t[k1p1 + 1][k2p1][k3p1 + 1] * nv3)
                           -k1p1 * ((k1 + 2) * (a_t[k1p1 + 1][k2p1][k3p1] * moments_str[0][0]
                                            + a_t[k1p1 + 2][k2p1 - 1][k3p1] * moments_str[0][1]
                                            + a_t[k1p1 + 2][k2p1][k3p1 - 1] * moments_str[0][2])
                                      + k2p1 * (a_t[k1p1][k2p1 + 1][k3p1] * moments_str[1][0]
                                              + a_t[k1p1 + 1][k2p1][k3p1] * moments_str[1][1]
                                          + a_t[k1p1 + 1][k2p1 + 1][k3p1 - 1] * moments_str[1][2])
                                    + k3p1 * (a_t[k1p1][k2p1][k3p1 + 1] * moments_str[2][0]
                                              + a_t[k1p1 + 1][k2p1 - 1][k3p1 + 1] * moments_str[2][1]
                                          + a_t[k1p1 + 1][k2p1][k3p1] * moments_str[2][2]))
                                          + k1p1 * a_t[k1p1 + 1][k2p1][k3p1] * MM +
                                          k1p1 * a_t[k1p1 + 1][k2p1][k3p1] * 2 * moments_str[0][0]
                                        + k2p1 * a_t[k1p1][k2p1 + 1][k3p1] * M1p2
                                        + k3p1 * a_t[k1p1][k2p1][k3p1 + 1] * M1p3)/3.0;
                    
                    velocity.val_Str[1] += (k2p1 * (k1p1 * a_t[k1p1 + 1][k2p1 + 1][k3p1] * nv1 +
                                                    (k2 + 2) * a_t[k1p1][k2p1 + 2][k3p1] * nv2 +
                                                     k3p1 * a_t[k1p1][k2p1 + 1][k3p1 + 1] * nv3)
                          -k2p1 * (k1p1 * (a_t[k1p1][k2p1 + 1][k3p1] * moments_str[0][0]
                                         + a_t[k1p1 + 1][k2p1][k3p1] * moments_str[0][1]
                                         + a_t[k1p1 + 1][k2p1 + 1][k3p1 - 1] * moments_str[0][2])
                                  +(k2 + 2) * (a_t[k1p1 - 1][k2p1 + 2][k3p1] * moments_str[1][0]
                                             + a_t[k1p1][k2p1 + 1][k3p1] * moments_str[1][1]
                                             + a_t[k1p1][k2p1 + 2][k3p1 - 1] * moments_str[1][2])
                                   +k3p1 * (a_t[k1p1 - 1][k2p1 + 1][k3p1 + 1] * moments_str[2][0]
                                          + a_t[k1p1][k2p1][k3p1 + 1] * moments_str[2][1]
                                          + a_t[k1p1][k2p1 + 1][k3p1] * moments_str[2][2]))
                                  +k2p1 * a_t[k1p1][k2p1 + 1][k3p1] * MM +
                                   k1p1 * a_t[k1p1 + 1][k2p1][k3p1] * M1p2 +
                                   k2p1 * a_t[k1p1][k2p1 + 1][k3p1] * 2 * moments_str[1][1]
                                  +k3p1 * a_t[k1p1][k2p1][k3p1 + 1] * M2p3)/3.0;
                    
                    velocity.val_Str[2] += (k3p1 * (k1p1 * a_t[k1p1 + 1][k2p1][k3p1 + 1] * nv1 +
                                                    k2p1 * a_t[k1p1][k2p1 + 1][k3p1 + 1] * nv2 +
                                                    (k3 + 2) * a_t[k1p1][k2p1][k3p1 + 2] * nv3)
                            - k3p1 * ( k1p1 * (a_t[k1p1][k2p1][k3p1 + 1] * moments_str[0][0]
                                     + a_t[k1p1 + 1][k2p1 - 1][k3p1 + 1] * moments_str[0][1]
                                             + a_t[k1p1 + 1][k2p1][k3p1] * moments_str[0][2])
                                      + k2p1 * (a_t[k1p1 - 1][k2p1 + 1][k3p1 + 1] * moments_str[1][0]
                                                + a_t[k1p1][k2p1][k3p1 + 1] * moments_str[1][1]
                                                + a_t[k1p1][k2p1 + 1][k3p1] * moments_str[1][2])
                                  + (k3 + 2) * (a_t[k1p1 - 1][k2p1][k3p1 + 2] * moments_str[2][0]
                                                    + a_t[k1p1][k2p1 - 1][k3p1 + 2] * moments_str[2][1]
                                                    + a_t[k1p1][k2p1][k3p1 + 1] * moments_str[2][2]))
                                        + k3p1 * a_t[k1p1][k2p1][k3p1 + 1] * MM
                                        +k1p1 * a_t[k1p1 + 1][k2p1][k3p1] * M1p3
                                        + k2p1 * a_t[k1p1][k2p1 + 1][k3p1] * M2p3
                                        +k3p1 * a_t[k1p1][k2p1][k3p1 + 1] * 2 * moments_str[2][2])/3.0;
                    // stokeslet
                    moments_sto[0] = tree[panel_index].moments_Sto[0][kk] ;
                    moments_sto[1] = tree[panel_index].moments_Sto[1][kk] ;
                    moments_sto[2] = tree[panel_index].moments_Sto[2][kk] ;
                    
                    term2 = tp_x * moments_sto[0] + tp_y * moments_sto[1] + tp_z * moments_sto[2];
                    
                    s = a_t[k1p1][k2p1][k3p1];
                    
                    velocity.val_Sto[0] += (1 - k1) * s * moments_sto[0] + k1p1 * (a_t[k1 + 2][k2p1][k3p1] * term2 - a_t[k1 + 2][k2][k3p1] * moments_sto[1] - a_t[k1 + 2][k2p1][k3] * moments_sto[2]);
                    
                    velocity.val_Sto[1] += (1 - k2) * s * moments_sto[1] + k2p1 * (a_t[k1p1][k2 + 2][k3p1] * term2 - a_t[k1][k2 + 2][k3p1] * moments_sto[0] - a_t[k1p1][k2 + 2][k3] * moments_sto[2]);
                    
                    velocity.val_Sto[2] += (1 - k3) * s * moments_sto[2] + k3p1 * (a_t[k1p1][k2p1][k3 + 2] * term2 - a_t[k1][k2p1][k3 + 2] * moments_sto[0] - a_t[k1p1][k2][k3 + 2] * moments_sto[1]);
                    
                    
                    
                    
                    
                }
            }
        }
}

    else
	{
		limit_1 = tree[panel_index].members[0];
		limit_2 = tree[panel_index].members[1];
		if (limit_2 - limit_1 < N0) //  otherwise, if cluster is a leaf, use direct sum
		{
			for (size_t jj = limit_1; jj <= limit_2; jj++)
			{
                if (jj == particle_index)
                    continue;
				dx = p_x - particles.x[jj];
                dy = p_y - particles.y[jj];
                dz = p_z - particles.z[jj];
                
                ff_sto[0] = lambda_Sto[0][jj];
                ff_sto[1] = lambda_Sto[1][jj];
                ff_sto[2] = lambda_Sto[2][jj];
               
                ff_str[0] = lambda_Str[0][jj];
                ff_str[1] = lambda_Str[1][jj];
                ff_str[2] = lambda_Str[2][jj];
                
                rnorm[0] = normal[0][jj];
                rnorm[1] = normal[1][jj];
                rnorm[2] = normal[2][jj];
                
                R = sqrt(dx * dx + dy * dy + dz * dz);
                Rinv = 1.0/R;
                velocity.val_Sto[0] += ff_sto[0] * Rinv;
                velocity.val_Sto[1] += ff_sto[1] * Rinv;
                velocity.val_Sto[2] += ff_sto[2] * Rinv;
                Rinv3 = Rinv * Rinv * Rinv;
                DD_sto = ff_sto[0] * dx + ff_sto[1] * dy + ff_sto[2] * dz;
                DD_sto = DD_sto * Rinv3;
                velocity.val_Sto[0] +=  DD_sto * dx;
                velocity.val_Sto[1] +=  DD_sto * dy;
                velocity.val_Sto[2] +=  DD_sto * dz;
               
                DD_str = (ff_str[0] * dx + ff_str[1] * dy + ff_str[2] * dz) *
                (rnorm[0] * dx + rnorm[1] * dy + rnorm[2] * dz);
                Rinv5 = Rinv3 * Rinv * Rinv;
                DD_str = DD_str * Rinv5;
                velocity.val_Str[0] += DD_str * dx;
                velocity.val_Str[1] += DD_str * dy;
                velocity.val_Str[2] += DD_str * dz;
               
                
			}
		}// end if
		else // othervise, if cluster is not a leaf, look at children
		{
			for (int index1 = 0; index1 < 3; index1++)
            {
                velocity.val_Str[index1] = 0.0;
                velocity.val_Sto[index1] = 0.0;
            }
            
			size_t length = tree[panel_index].children.size();
			for (size_t i = 0; i < length; i++)
			{
				size_t index = tree[panel_index].children[i];
                temp_result = Comput_RBF(lambda_Sto,lambda_Str,normal, particles, particle_index, index, dep + 1);
                
                for(int index_1 = 0; index_1 < 3; index_1++)
                {
                    velocity.val_Str[index_1] = velocity.val_Str[index_1] + temp_result.val_Str[index_1];
                    velocity.val_Sto[index_1] = velocity.val_Sto[index_1] + temp_result.val_Sto[index_1];
                }
			}
            
		}
        
	}
   /* cout << "velocity.val_Str[0] = "<< velocity.val_Str[0] << endl;
    cout << "velocity.val_Str[1] = "<< velocity.val_Str[1] << endl;
    cout << "velocity.val_Str[2] = "<< velocity.val_Str[2] << endl;
    
    cout << "velocity.val_Sto[0] = "<< velocity.val_Sto[0] << endl;
    cout << "velocity.val_Sto[1] = "<< velocity.val_Sto[1] << endl;
    cout << "velocity.val_Sto[2] = "<< velocity.val_Sto[2] << endl;*/
    
    return velocity;
}

//*****************************************************************************//
int main()
{
	struct xyz particles(N_cube);
   
	cout << "P is " << P << endl;
	cout << "N is " << N << endl;
	cout << "theta is " << sqrt(sq_theta) << endl;
    cout << "N0 is " << N0 << endl;
	   
    // *****************  Generate Random data *******************************
    
    FILE * fp;
    double x1, x2, x3,x4,x5,x6,x7,x8,x9;
    
    int count;
    char S_data_file[64] = {0};
    sprintf(S_data_file, "./data/Points_N%d", N);
    fp = fopen(S_data_file, "r");
	count = -1;
	if (fp == NULL)
	{
		cout << "Cannot open random points file" << endl;
		getchar();
		exit(0);
	}
    
	while (true)
	{
		count++;
        fscanf(fp,"%lf%lf%lf%lf%lf", &x1, &x2, &x3,&x4,&x5); // x4: area, x5 Farea
		if (feof(fp))
			break;
		if (count >= N_cube)
		{
			cout << "S Out of range" << endl;
			exit(0);
		}
        
		particles.x[count] = x1;
		particles.y[count] = x2;
		particles.z[count] = x3;
		particles.index[count] = -1;
		particles.old_index[count] = count;
	}


    
    double *normal[3]; // normal direction at each point
	normal[0] = new double[N_cube];
	normal[1] = new double[N_cube];
	normal[2] = new double[N_cube];
    
	for (size_t i = 0; i < N_cube; i++)
    {
        normal[0][i] = particles.x[i];
        normal[1][i] = particles.y[i];
        normal[2][i] = particles.z[i];
	}
    
    // weights

	double *lambda_Sto[3];
	lambda_Sto[0] = new double[N_cube];
	lambda_Sto[1] = new double[N_cube];
	lambda_Sto[2] = new double[N_cube];
    
    char lambda_Sto_data_file[64] = {0};
    sprintf(lambda_Sto_data_file, "./data/lambda_sto_%d.txt", N);
    fp = fopen(lambda_Sto_data_file, "r");
    
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
			cout << "lambda Sto Out of range" << endl;
			exit(0);
		}
        
		lambda_Sto[0][count] = x1;
        lambda_Sto[1][count] = x2;
        lambda_Sto[2][count] = x3;
		
	}
    
    char lambda_Str_data_file[64] = {0};
    sprintf(lambda_Str_data_file, "./data/lambda_str_%d.txt", N);
    fp = fopen(lambda_Str_data_file, "r");
    
	count = -1;
	if (fp == NULL)
	{
		cout << "Cannot open lambda file" << endl;
		getchar();
		exit(0);
	}
    
    double *lambda_Str[3];
	lambda_Str[0] = new double[N_cube];
	lambda_Str[1] = new double[N_cube];
	lambda_Str[2] = new double[N_cube];
    
	while (true)
	{
		count++;
		fscanf(fp,"%lf%lf%lf", &x1, &x2, &x3);
		if (feof(fp))
			break;
		if (count >= N_cube)
		{
			cout << "lambda Str Out of range" << endl;
			exit(0);
		}
        
		lambda_Str[0][count] = x1;
        lambda_Str[1][count] = x2;
        lambda_Str[2][count] = x3;
	}
    
    


	//***************** Set up tree *******************************
	long Start_total, Start_btree;
    long End_total, End_btree;

	Start_total = getTickCount(); // Get currenct CPU time
    Start_btree = getTickCount();
    
    build_tree_init(particles);
	build_tree_3D_Recursive(0, particles, 0);
    ShrinkCluster(particles);
    
    End_btree = getTickCount();
	
	//***************** Compute moment for each panel **************
	size_t size = tree.size();
	for (size_t i = 0; i < size; i++)
		Panel_Moment_Taylor(i, lambda_Sto, lambda_Str,normal, particles,  tree[i].moments_Sto,tree[i].moments_Str);

	//***************** Compute Velocity ***************************
    
    vec_3d *tree_app = new vec_3d[N_cube];
	
	for (int i = 0; i < N_cube; i++)
	{
		tree_app[i] = Comput_RBF(lambda_Sto,lambda_Str,normal, particles, i, 0);
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
	
	double theta = sqrt(sq_theta);
	char file_name[256];
	sprintf(file_name, "BR_P%d_N%d_theta%.2f_N0%d_%s", P, N, theta, N0, time_buffer);
	ofstream output_file(file_name);
	
	output_file << " ===== No box shrink ===========" << endl;
	output_file << "P is " << P << endl;
	output_file << "N is " << N << endl;
	output_file << "theta is " << sqrt(sq_theta) << endl;
    output_file << "N0 is " << N0 << endl;
	
	cout << "N_cube is "<< N_cube << endl;

    cout << "treecode_cpu_time " << treecode_cpu_time << endl;
	cout << "build tree time is " << End_btree - Start_btree << endl;

    output_file << "treecode_cpu_time " << treecode_cpu_time << endl;
	output_file << "build tree time is " << End_btree - Start_btree << endl;

	//***************** Director summation and L_2 Error *****************
    int N_sub = N_cube;
    
	double *sum_ex_Str[3];
    double *sum_ex_Sto[3];
    
    for (int i = 0; i < 3; i++)
    {
        sum_ex_Sto[i] = new double[N_sub];
        sum_ex_Str[i] = new double[N_sub];
    }
    
	//***************** compute v_true first N_sub points here ******************************
    
    long Start_ds, End_ds;
    Start_ds = getTickCount(); // Get currenct CPU time
    
	double x,y,z;
	double px,py,pz;
	double dx,dy,dz;
	double sum_p_Str[3] = {0.0};
	double sum_p_Sto[3] = {0.0};
	double R,Rinv,Rinv3,Rinv5;
    double ff_sto[3];
    double ff_str[3];
    double rnorm[3];
    double DD_sto;
    double DD_str;
    
    
	for (size_t i = 0; i < N_sub; i++)
	{
		x = particles.x[i];
		y = particles.y[i];
		z = particles.z[i];
		for (size_t j = 0; j < N_cube; j++)
		{
            if (i!=j)
            {
                px = particles.x[j];
                py = particles.y[j];
                pz = particles.z[j];
                dx = x - px;
                dy = y - py;
                dz = z - pz;
                R = sqrt(dx * dx + dy * dy + dz * dz);
        
                ff_sto[0] = lambda_Sto[0][j];
                ff_sto[1] = lambda_Sto[1][j];
                ff_sto[2] = lambda_Sto[2][j];
            
                ff_str[0] = lambda_Str[0][j];
                ff_str[1] = lambda_Str[1][j];
                ff_str[2] = lambda_Str[2][j];
            
                rnorm[0] = normal[0][j];
                rnorm[1] = normal[1][j];
                rnorm[2] = normal[2][j];
                
                Rinv = 1.0/R;
                sum_p_Sto[0] += ff_sto[0] * Rinv;
                sum_p_Sto[1] += ff_sto[1] * Rinv;
                sum_p_Sto[2] += ff_sto[2] * Rinv;
                Rinv3 = Rinv * Rinv * Rinv;
                DD_sto = ff_sto[0] * dx + ff_sto[1] * dy + ff_sto[2] * dz;
                DD_sto = DD_sto * Rinv3;
                sum_p_Sto[0] +=  DD_sto * dx;
                sum_p_Sto[1] +=  DD_sto * dy;
                sum_p_Sto[2] +=  DD_sto * dz;
                
                DD_str = (ff_str[0] * dx + ff_str[1] * dy + ff_str[2] * dz) *
                (rnorm[0] * dx + rnorm[1] * dy + rnorm[2] * dz);
                Rinv5 = Rinv3 * Rinv * Rinv;
                DD_str = DD_str * Rinv5;
                sum_p_Str[0] += DD_str * dx;
                sum_p_Str[1] += DD_str * dy;
                sum_p_Str[2] += DD_str * dz;
            }
        }

        sum_ex_Str[0][i] = sum_p_Str[0];
        sum_ex_Str[1][i] = sum_p_Str[1];
        sum_ex_Str[2][i] = sum_p_Str[2];
    
        sum_ex_Sto[0][i] = sum_p_Sto[0];
        sum_ex_Sto[1][i] = sum_p_Sto[1];
        sum_ex_Sto[2][i] = sum_p_Sto[2];
        
       
    
        for (int j = 0; j < 3; j++)
        {
            sum_p_Sto[j] = 0;
            sum_p_Str[j] = 0;
        }

}

    End_ds = getTickCount(); // Get currenct CPU time
    long ds_cpu_time;
    ds_cpu_time = End_ds - Start_ds;
    
    cout << "ds time is " << ds_cpu_time << endl;
    
    output_file << "ds_cpu_time " << ds_cpu_time << endl;
    
    //====== L_infty Err======================================================

	double E8_Str = 0.0;
    double E8_Sto = 0.0;
    
    double temp_d_Str = 0.0;
    double temp_n_Str = 0.0;
    double max_d_Str = 0.0;
    double max_n_Str = 0.0;
    
    double temp_d_Sto = 0.0;
    double temp_n_Sto = 0.0;
    double max_d_Sto = 0.0;
    double max_n_Sto = 0.0;
    
	for (size_t i = 0; i < N_sub; i++)
	{
		for (int index1 = 0; index1<3; index1++)
        {
            temp_n_Str += (tree_app[i].val_Str[index1] - sum_ex_Str[index1][i]) * (tree_app[i].val_Str[index1] - sum_ex_Str[index1][i]);
            temp_d_Str += sum_ex_Str[index1][i] * sum_ex_Str[index1][i];
            
        }
        temp_n_Str= sqrt(temp_n_Str);
        temp_d_Str = sqrt(temp_d_Str);
        if (temp_n_Str > max_n_Str)
            max_n_Str = temp_n_Str;
        if (temp_d_Str > max_d_Str)
            max_d_Str = temp_d_Str;
        temp_d_Str = 0.0;
        temp_n_Str = 0.0;
    }
    
    E8_Str = max_n_Str/max_d_Str;
    
	cout << "Str E8 is " << E8_Str << endl;
    output_file << "Str E8 is " << E8_Str << endl;
    
    for (size_t i = 0; i < N_sub; i++)
	{
		for (int index1 = 0; index1<3; index1++)
        {
            temp_n_Sto += (tree_app[i].val_Sto[index1] - sum_ex_Sto[index1][i]) * (tree_app[i].val_Sto[index1] - sum_ex_Sto[index1][i]);
            temp_d_Sto += sum_ex_Sto[index1][i] * sum_ex_Sto[index1][i];
        }
        temp_n_Sto= sqrt(temp_n_Sto);
        temp_d_Sto = sqrt(temp_d_Sto);
        if (temp_n_Sto > max_n_Sto)
            max_n_Sto = temp_n_Sto;
        if (temp_d_Sto > max_d_Sto)
            max_d_Sto = temp_d_Sto;
        temp_d_Sto = 0.0;
        temp_n_Sto = 0.0;
    }
    
    E8_Sto = max_n_Sto/max_d_Sto;
    
	cout << "Sto E8 is " << E8_Sto << endl;
    output_file << "Sto E8 is " << E8_Sto << endl;
	
    //******** L2 Error : extend from RBF paper *******************
    double E2_Str = 0.0;
    double sum_d_ex_Str = 0.0;
    double sum_n_ex_Str = 0.0;
    
    double E2_Sto = 0.0;
    double sum_d_ex_Sto = 0.0;
    double sum_n_ex_Sto = 0.0;
    
    for (size_t i = 0; i < N_sub; i++)
	{
        for (int index1 = 0; index1<3; index1++)
        {
            sum_n_ex_Str += (tree_app[i].val_Str[index1] - sum_ex_Str[index1][i]) * (tree_app[i].val_Str[index1] - sum_ex_Str[index1][i]);
            sum_d_ex_Str += sum_ex_Str[index1][i] * sum_ex_Str[index1][i];
        }
    }
    
    for (size_t i = 0; i < N_sub; i++)
	{
        for (int index1 = 0; index1<3; index1++)
        {
            sum_n_ex_Sto += (tree_app[i].val_Sto[index1] - sum_ex_Sto[index1][i]) * (tree_app[i].val_Sto[index1] - sum_ex_Sto[index1][i]);
            sum_d_ex_Sto += sum_ex_Sto[index1][i] * sum_ex_Sto[index1][i];
        }
    }
    
    E2_Str = sqrt(sum_n_ex_Str/sum_d_ex_Str);
    E2_Sto = sqrt(sum_n_ex_Sto/sum_d_ex_Sto);
    
    //********** total error ***************************
    double total_error;
    double sum_d = 0.0;
    double sum_n = 0.0;
    for (size_t i = 0; i < N_sub; i++)
    {
        for (int index1 = 0; index1<3; index1++)
        {
            sum_d += (sum_ex_Str[index1][i] + sum_ex_Sto[index1][i]) * (sum_ex_Str[index1][i] + sum_ex_Sto[index1][i]);
            sum_n += (sum_ex_Sto[index1][i] - tree_app[i].val_Sto[index1]  + sum_ex_Str[index1][i] - tree_app[i].val_Str[index1]) * (sum_ex_Sto[index1][i] - tree_app[i].val_Sto[index1]  + sum_ex_Str[index1][i] - tree_app[i].val_Str[index1]);
        }
    }
    total_error = sqrt(sum_n/sum_d);
    
    cout << "E2_Str is " << E2_Str << endl;
    output_file << "E2_Str is " << E2_Str << endl;
    cout << "E2_Sto is " << E2_Sto << endl;
    output_file << "E2_Sto is " << E2_Sto << endl;
    cout << "total error is " << total_error << endl;
    output_file << "total error is " << total_error << endl;
    
    
    //*********************************************************
    
    delete [] tree_app;
    
    for (int i = 0; i < 3; i++)
    {
        delete [] sum_ex_Sto[i];
        delete [] sum_ex_Str[i];
    }
    
    
    delete [] normal[0];
	delete [] normal[1];
	delete [] normal[2];
    
	delete [] lambda_Sto[0];
	delete [] lambda_Sto[1];
	delete [] lambda_Sto[2];
    
    delete [] lambda_Str[0];
	delete [] lambda_Str[1];
	delete [] lambda_Str[2];
    

    
	output_file.close();

	return 0;
}
