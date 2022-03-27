#include "gauss_quad.h"
void quad_points_2d(double** p, double x[], double y[], double z[])
{
	int num_quad = 49, i;
    const double* points = 0;
    points = QUAD_2D_P15_pts;
    double x_f[3], y_f[3], z_f[3];

    //face 1:1 2 3
    x_f[0] = x[0];  x_f[1] = x[1];  x_f[2] = x[2];
    y_f[0] = y[0];  y_f[1] = y[1];  y_f[2] = y[2];
    z_f[0] = z[0];  z_f[1] = z[1];  z_f[2] = z[2];
    for (i = 0; i < num_quad; ++i)
    {
        p[0][i] = x_f[0] * points[3 * i] + x_f[1] * points[3 * i + 1] + x_f[2] * points[3 * i + 2];
        p[1][i] = y_f[0] * points[3 * i] + y_f[1] * points[3 * i + 1] + y_f[2] * points[3 * i + 2];
        p[2][i] = z_f[0] * points[3 * i] + z_f[1] * points[3 * i + 1] + z_f[2] * points[3 * i + 2];
    }

    //face 2:1 2 4
    x_f[0] = x[0]; x_f[1] = x[1]; x_f[2] = x[3]; 
    y_f[0] = y[0]; y_f[1] = y[1]; y_f[2] = y[3];
    z_f[0] = z[0]; z_f[1] = z[1]; z_f[2] = z[3];
    for (i = 0; i < num_quad; ++i)
    {
        p[3][i] = x_f[0] * points[3 * i] + x_f[1] * points[3 * i + 1] + x_f[2] * points[3 * i + 2];
        p[4][i] = y_f[0] * points[3 * i] + y_f[1] * points[3 * i + 1] + y_f[2] * points[3 * i + 2];
        p[5][i] = z_f[0] * points[3 * i] + z_f[1] * points[3 * i + 1] + z_f[2] * points[3 * i + 2];
    }

    //face 3:2 3 4
    x_f[0] = x[1]; x_f[1] = x[2]; x_f[2] = x[3];
    y_f[0] = y[1]; y_f[1] = y[2]; y_f[2] = y[3]; 
    z_f[0] = z[1]; z_f[1] = z[2]; z_f[2] = z[3];
    for (i = 0; i < num_quad; ++i)
    {
        p[6][i] = x_f[0] * points[3 * i] + x_f[1] * points[3 * i + 1] + x_f[2] * points[3 * i + 2];
        p[7][i] = y_f[0] * points[3 * i] + y_f[1] * points[3 * i + 1] + y_f[2] * points[3 * i + 2];
        p[8][i] = z_f[0] * points[3 * i] + z_f[1] * points[3 * i + 1] + z_f[2] * points[3 * i + 2];
    }

    //face 4:1 3 4
    x_f[0] = x[0]; x_f[1] = x[2]; x_f[2] = x[3];  
    y_f[0] = y[0]; y_f[1] = y[2]; y_f[2] = y[3]; 
    z_f[0] = z[0]; z_f[1] = z[2]; z_f[2] = z[3];
    for (i = 0; i < num_quad; ++i)
    {
        p[9][i]  = x_f[0] * points[3 * i] + x_f[1] * points[3 * i + 1] + x_f[2] * points[3 * i + 2];
        p[10][i] = y_f[0] * points[3 * i] + y_f[1] * points[3 * i + 1] + y_f[2] * points[3 * i + 2];
        p[11][i] = z_f[0] * points[3 * i] + z_f[1] * points[3 * i + 1] + z_f[2] * points[3 * i + 2];
    }
}

void triangle_areas_on_tetra(double* areas, double x[], double y[], double z[])
{
    double len[3], x_f[3], y_f[3], z_f[3];

    //face 1:1 2 3
    x_f[0] = x[0]; x_f[1] = x[1]; x_f[2] = x[2];
    y_f[0] = y[0]; y_f[1] = y[1]; y_f[2] = y[2];
    z_f[0] = z[0]; z_f[1] = z[1]; z_f[2] = z[2];

    len[0] = sqrt(pow((x_f[1] - x_f[0]), 2.0) + pow((y_f[1] - y_f[0]), 2.0) + pow((z_f[1] - z_f[0]), 2.0));
    len[1] = sqrt(pow((x_f[2] - x_f[1]), 2.0) + pow((y_f[2] - y_f[1]), 2.0) + pow((z_f[2] - z_f[1]), 2.0));
    len[2] = sqrt(pow((x_f[2] - x_f[0]), 2.0) + pow((y_f[2] - y_f[0]), 2.0) + pow((z_f[2] - z_f[0]), 2.0));
	areas[0] = sqrt((len[0] + len[1] + len[2]) * (len[0] + len[1] - len[2]) * (len[0] + len[2] - len[1]) * (len[1] + len[2] - len[0])) / 4;

    //face 2:1 2 4
    x_f[0] = x[0]; x_f[1] = x[1]; x_f[2] = x[3];
    y_f[0] = y[0]; y_f[1] = y[1]; y_f[2] = y[3];
    z_f[0] = z[0]; z_f[1] = z[1]; z_f[2] = z[3];

    len[0] = sqrt(pow((x_f[1] - x_f[0]), 2.0) + pow((y_f[1] - y_f[0]), 2.0) + pow((z_f[1] - z_f[0]), 2.0));
    len[1] = sqrt(pow((x_f[2] - x_f[1]), 2.0) + pow((y_f[2] - y_f[1]), 2.0) + pow((z_f[2] - z_f[1]), 2.0));
    len[2] = sqrt(pow((x_f[2] - x_f[0]), 2.0) + pow((y_f[2] - y_f[0]), 2.0) + pow((z_f[2] - z_f[0]), 2.0));
	areas[1] = sqrt((len[0] + len[1] + len[2]) * (len[0] + len[1] - len[2]) * (len[0] + len[2] - len[1]) * (len[1] + len[2] - len[0])) / 4;

    //face 3:2 3 4
    x_f[0] = x[1]; x_f[1] = x[2]; x_f[2] = x[3];
    y_f[0] = y[1]; y_f[1] = y[2]; y_f[2] = y[3];
    z_f[0] = z[1]; z_f[1] = z[2]; z_f[2] = z[3];

    len[0] = sqrt(pow((x_f[1] - x_f[0]), 2.0) + pow((y_f[1] - y_f[0]), 2.0) + pow((z_f[1] - z_f[0]), 2.0));
    len[1] = sqrt(pow((x_f[2] - x_f[1]), 2.0) + pow((y_f[2] - y_f[1]), 2.0) + pow((z_f[2] - z_f[1]), 2.0));
    len[2] = sqrt(pow((x_f[2] - x_f[0]), 2.0) + pow((y_f[2] - y_f[0]), 2.0) + pow((z_f[2] - z_f[0]), 2.0));
	areas[2] = sqrt((len[0] + len[1] + len[2]) * (len[0] + len[1] - len[2]) * (len[0] + len[2] - len[1]) * (len[1] + len[2] - len[0])) / 4;

    //face 4:1 3 4
    x_f[0] = x[0]; x_f[1] = x[2]; x_f[2] = x[3];
    y_f[0] = y[0]; y_f[1] = y[2]; y_f[2] = y[3];
    z_f[0] = z[0]; z_f[1] = z[2]; z_f[2] = z[3];

    len[0] = sqrt(pow((x_f[1] - x_f[0]), 2.0) + pow((y_f[1] - y_f[0]), 2.0) + pow((z_f[1] - z_f[0]), 2.0));
    len[1] = sqrt(pow((x_f[2] - x_f[1]), 2.0) + pow((y_f[2] - y_f[1]), 2.0) + pow((z_f[2] - z_f[1]), 2.0));
    len[2] = sqrt(pow((x_f[2] - x_f[0]), 2.0) + pow((y_f[2] - y_f[0]), 2.0) + pow((z_f[2] - z_f[0]), 2.0));
	areas[3] = sqrt((len[0] + len[1] + len[2]) * (len[0] + len[1] - len[2]) * (len[0] + len[2] - len[1]) * (len[1] + len[2] - len[0])) / 4;
}