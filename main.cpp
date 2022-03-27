#include <iostream>
#include "gauss_quad.h"
using namespace std;

int main()
{
    double tetra1_x[4] = { 0, 1,-1, 0 };
	double tetra1_y[4] = { 1, 0, 0, 0 };
	double tetra1_z[4] = { 0, 0, 0, 1 };
	double tetra2_x[4] = { 0,-1, 1, 0 };
	double tetra2_y[4] = { 0, 0, 0,-1 };
	double tetra2_z[4] = { 1, 0, 0, 0 };

	int order = 15;
	int num_quad_2d = 49;
	double areas[4];
	const double* quad_2D_wts = QUAD_2D_P15_wts;
	double** quad_pts_face = new double* [12];
	for (int i = 0; i < 12; ++i) {
		quad_pts_face[i] = new double[num_quad_2d];
	}
	quad_points_2d(quad_pts_face, tetra1_x, tetra1_y, tetra1_z);
	triangle_areas_on_tetra(areas, tetra1_x, tetra1_y, tetra1_z);

	double ret1 = 0, ret2 = 0;
#pragma region method 1
	int idxf = 2; // the interface between two tetrahedrons is face[2]
	double* f1 = new double[num_quad_2d];
	double* f2 = new double[num_quad_2d];
	double* dyf2 = new double[num_quad_2d];
	for (int i = 0; i < num_quad_2d; ++i) {
		double x = quad_pts_face[idxf * 3][i];
		double y = quad_pts_face[idxf * 3 + 1][i];
		double z = quad_pts_face[idxf * 3 + 2][i];
		f1[i] = 2 * x + 4 * y + z;
		f2[i] = 3 * x + 2 * y * y - z * y;
		dyf2[i] = 4 * y - z;
	}
	for (int i = 0; i < num_quad_2d; ++i) {
		ret1 += f1[i] * dyf2[i] * quad_2D_wts[i] * areas[idxf];
	}
#pragma endregion 

#pragma region method 2
    double ref_tetra_x[4] = { -1, 1,-1,-1 };
    double ref_tetra_y[4] = { -1,-1, 1,-1 };
    double ref_tetra_z[4] = { -1,-1,-1, 1 };
    double** ref_face_quad_pts = new double* [12];
    for (int i = 0; i < 12; ++i) {
        ref_face_quad_pts[i] = new double[num_quad_2d];
    }
    double ref_areas[4];
    quad_points_2d(ref_face_quad_pts, ref_tetra_x, ref_tetra_y, ref_tetra_z);
    triangle_areas_on_tetra(ref_areas, ref_tetra_x, ref_tetra_y, ref_tetra_z);
    double* ref_f1 = new double[num_quad_2d];
    double* drf2 = new double[num_quad_2d];
    double* dsf2 = new double[num_quad_2d];
    double* dtf2 = new double[num_quad_2d];
    const double* ref_quad_2D_wts = QUAD_2D_P15_wts;

    // r s t to x y z on gauss integral point
    int idxf1 = 2; int idxf2 = 0;
	double* r1 = new double[num_quad_2d];
	double* s1 = new double[num_quad_2d];
	double* t1 = new double[num_quad_2d];
	double* r2 = new double[num_quad_2d];
	double* s2 = new double[num_quad_2d];
	double* t2 = new double[num_quad_2d];

    for (int i = 0; i < num_quad_2d; ++i) {
        r1[i] = ref_face_quad_pts[idxf1 * 3][i];
        s1[i] = ref_face_quad_pts[idxf1 * 3 + 1][i];
        t1[i] = ref_face_quad_pts[idxf1 * 3 + 2][i];
    }

    for (int i = 0; i < num_quad_2d; ++i) {
        r2[i] = ref_face_quad_pts[idxf2 * 3][i];
        s2[i] = ref_face_quad_pts[idxf2 * 3 + 1][i];
        t2[i] = ref_face_quad_pts[idxf2 * 3 + 2][i];
    }

	double* x1 = new double[num_quad_2d];
	double* y1 = new double[num_quad_2d];
	double* z1 = new double[num_quad_2d];
	double* x2 = new double[num_quad_2d];
	double* y2 = new double[num_quad_2d];
	double* z2 = new double[num_quad_2d];

	for (int i = 0; i != num_quad_2d; ++i)
    {
		x1[i] = 0.5 * (-(1 + r1[i] + s1[i] + t1[i]) * tetra1_x[0] 
            + (1 + r1[i]) * tetra1_x[1] + (1 + s1[i]) * tetra1_x[2] + (1 + t1[i]) * tetra1_x[3]);
		y1[i] = 0.5 * (-(1 + r1[i] + s1[i] + t1[i]) * tetra1_y[0] 
            + (1 + r1[i]) * tetra1_y[1] + (1 + s1[i]) * tetra1_y[2] + (1 + t1[i]) * tetra1_y[3]);
		z1[i] = 0.5 * (-(1 + r1[i] + s1[i] + t1[i]) * tetra1_z[0] 
            + (1 + r1[i]) * tetra1_z[1] + (1 + s1[i]) * tetra1_z[2] + (1 + t1[i]) * tetra1_z[3]);
		x2[i] = 0.5 * (-(1 + r2[i] + s2[i] + t2[i]) * tetra2_x[0] 
            + (1 + r2[i]) * tetra2_x[1] + (1 + s2[i]) * tetra2_x[2] + (1 + t2[i]) * tetra2_x[3]);
		y2[i] = 0.5 * (-(1 + r2[i] + s2[i] + t2[i]) * tetra2_y[0] 
            + (1 + r2[i]) * tetra2_y[1] + (1 + s2[i]) * tetra2_y[2] + (1 + t2[i]) * tetra2_y[3]);
		z2[i] = 0.5 * (-(1 + r2[i] + s2[i] + t2[i]) * tetra2_z[0] 
            + (1 + r2[i]) * tetra2_z[1] + (1 + s2[i]) * tetra2_z[2] + (1 + t2[i]) * tetra2_z[3]);
    }

    double tol = 1e-8;
    int* quad_pts_map = new int[num_quad_2d];
    for (int i = 0; i < num_quad_2d; ++i) {
        quad_pts_map[i] = i;
        for (int j = 0; j < num_quad_2d; ++j) {
            double distmp = sqrt((x1[i] - x2[j]) * (x1[i] - x2[j])
                + (y1[i] - y2[j]) * (y1[i] - y2[j])
                + (z1[i] - z2[j]) * (z1[i] - z2[j]));
            if (distmp < tol) {
                quad_pts_map[i] = j; break;
            }
        }
    }

    double sJ = areas[idxf1] / ref_areas[idxf1];
    double ry = 1, sy = 1, ty = -2;
    for (int i = 0; i < num_quad_2d; ++i) {
        double r_1 = ref_face_quad_pts[idxf1 * 3][i];
        double s_1 = ref_face_quad_pts[idxf1 * 3 + 1][i];
        double t_1 = ref_face_quad_pts[idxf1 * 3 + 2][i];
        double r_2 = ref_face_quad_pts[idxf2 * 3][quad_pts_map[i]];
        double s_2 = ref_face_quad_pts[idxf2 * 3 + 1][quad_pts_map[i]];
        double t_2 = ref_face_quad_pts[idxf2 * 3 + 2][quad_pts_map[i]];
        drf2[i] = -1.5 - 0.25 * (t_2 + 1);
        dsf2[i] = 1.5 - 0.25 * (t_2 + 1);
        dtf2[i] = 0.5 * (t_2 + 1) - 0.25 * (r_2 + s_2);
        ref_f1[i] = (-r_1 - 3 * s_1 - 1.5 * t_1 - 1.5);
    }

    for (int i = 0; i < num_quad_2d; ++i) {
        ret2 += ref_f1[i] * (drf2[i] * ry + dsf2[i] * sy
            + dtf2[i] * ty) * quad_2D_wts[i] * ref_areas[idxf1];
    }

#pragma endregion 
    cout << " result1 = " << ret1 << endl;
    cout << " result2 = " << ret2 * sJ << endl;
	return 0;
}
