#ifndef __MASTERCELL_C__
#define __MASTERCELL_C__

template<const int Ndim>
struct masterCell
{
  const int dim;                // dimension 1, 2, or 3 ...
 
  int nsubs[Ndim+1];            // number of sub-cells of dimension 0,1,..,dim

  int subSize[Ndim+1];          // subSize[d] number of vertices in a d-cell

  int (*subs)(int,int,int);     // subs[d][i][.] = local node numbers of
                                // the i^th sub-cell of dimension d
  double volume;
};

int Point_Subs(int i, int j, int k)
{
  int  ss[1][1][1] =   
    {
      {{0}}
    };
 
  return( ss[i][j][k] );
}

masterCell<0> Point = 
  {0, 
   {1}, 
   {1},
   Point_Subs,
   1.0              // Point Mass
   };

int Interval_Subs(int i, int j, int k)
{
  int  ss[2][2][2] =   
    {
      {{0},{1}},
      {{0,1}}
    };
 
  return( ss[i][j][k] );
}

masterCell<1> Interval = 
  {1, 
   {2,1}, 
   {1,2},
   Interval_Subs,
   2.0              // lenght of master interval [-1,1]
   };

int Square_Subs(int i, int j, int k)
{
  int  ss[3][4][4] =   
    {
      {{0},{1},{2},{3}}, 
      {{0,1},{1,2},{2,3},{3,0}}, 
      {{0,1,2,3}}
    };
 
  return( ss[i][j][k] );
}

masterCell<2> Square = 
  {2, 
   {4,4,1}, 
   {1,2,4}, 
   Square_Subs,
   4.0              // area of master square [-1,1]^2
   };

int Triangle_Subs(int i, int j, int k)
{
  int  ss[3][3][3] =   
    {
      {{0},{1},{2}}, 
      {{0,1},{1,2},{2,0}},    // better is edge opposite 0, 1, 2
      {{0,1,2}}
    };
 
  return( ss[i][j][k] );
}

masterCell<2> Triangle = 
  {2, 
   {3,3,1}, 
   {1,2,3}, 
   Triangle_Subs,
   0.5              // area of unit triangle
   };

int Cube_Subs(int i, int j, int k)
{
  int  ss[4][12][8] =   
    {
      {{0},{1},{2},{3},{4},{5},{6},{7}}, 
      {{0,1}, {3,2}, {5,4}, {7,6}, {2,0}, {1,3}, 
       {6,4}, {5,7}, {0,4}, {1,5}, {2,6}, {3,7}},
      {{0,2,3,1}, {4,5,7,6}, {0,1,5,4}, {2,6,7,3}, {0,4,6,2}, {1,3,7,5}},
      {{0,1,2,3,4,5,6,7}}
    };
 
  return( ss[i][j][k] );
}

masterCell<3> Cube = 
  {3, 
   {8,12,6,1}, 
   {1, 2,4,8}, 
   Cube_Subs,
   8.0              // volume of master cube [-1,1]^3
   };

int Tetrahedron_Subs(int i, int j, int k)
{
  int  ss[4][6][4] =   
    {
      {{0},{1},{2},{3}}, 
      {{0,1}, {1,2}, {2,0}, {3,0}, {3,1}, {3,2}},
      {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}},
      {{0,1,2,3}}
    };
 
  return( ss[i][j][k] );
}

masterCell<3> Tetrahedron = 
  {3, 
   {4,6,4,1},        // # vertices, edges, faces, cells
   {1,2,3,4},        // # vertices each has
   Tetrahedron_Subs,
   1.0/6.0           // volume of unit tet
   };
#endif                          // !defined(__MASTERCELL_C__)
