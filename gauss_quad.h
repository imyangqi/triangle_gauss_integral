#pragma once
#include <math.h>

/* 2D */
#define Length(a)      (sizeof(a)/sizeof(a[0]))
#define Perm3(a)       double(a),double(a),double(a)
#define Dup3(w)        w
#define Perm21(a)      double(1.)-a-a,a,a,a,double(1.)-a-a,a,a,a,double(1.)-a-a
#define Dup21(w)       w,w,w
#define Perm111(a,b)   a,b,double(1.)-a-b,a,double(1.)-a-b,b,b,a,double(1.)-a-b,b,double(1.)-a-b,a,double(1.)-a-b,a,b,double(1.)-a-b,b,a
#define Dup111(w)      w,w,w,w,w,w

const double QUAD_2D_P15_wts[] = {
	Dup3(.02357126703190634206659321140821418),
	Dup21(.01517314955721170450311858877690239),
	Dup21(.01297600128392884154979521077280757),
	Dup21(.01706629596800615670942600046160914),
	Dup21(.04576001946273760698482638108892258),
	Dup111(.00222757447282223154006065426298478),
	Dup111(.02701014165986947101315702212247500),
	Dup111(.02608377963958756403057720483642768),
	Dup111(.01211015327702828337230795926322736),
	Dup111(.01564785059680444573399007149035058),
	Dup111(.03417088937929479242522512890637806)
};

const double QUAD_2D_P15_pts[Length(QUAD_2D_P15_wts) * 3] = {
   Perm3(.33333333333333333333333333333333333),
   Perm21(.11022229622834687297855264132259850),
   Perm21(.05197643301003435047003197947889073),
   Perm21(.49114565807532554119014945122395425),
   Perm21(.39315718888435884048226809785071794),
   Perm111(.03737440487572919066543605209836625,
		.96251835223001214880811969560396873),
   Perm111(.24824877798467321198263980694374938,
		.19316669854521416819773100288721521),
   Perm111(.20699402274830217740486528153682148,
		.08689590883549962551575259619781217),
   Perm111(.14854110526954708137688902238435510,
		.01743682539845430796259020511767948),
   Perm111(.30674237923596382376588728350286621,
		.01749251095825766163254977051260599),
   Perm111(.36703198754220473278855469116984882,
		.09034802175864556044634095119222305)
};

void quad_points_2d(double** p, double x[], double y[], double z[]);
void triangle_areas_on_tetra(double* Area, double x[], double y[], double z[]);