// Set up a mesh with each trignale affine equivalent to the standard 

#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

using namespace std;

int nex = 8, ney = 8;
double xmin = -1, xmax = 1, ymin = -1, ymax = 1;
int bcType[4] = {0,1,0,0};

int nelem;
int nodes;
int nbelem;

int main(int argc, char *argv[])
{
 int i,j,n;

 if(argc > 1)
   {
   istringstream ist(argv[1]);

   if( !(ist >> nex) )
     {
     cerr << "Invalid value for nex: " << argv[1] << ", terminating" << endl;
     throw("Invalid input");
     }

   if(nex <= 0) 
     {
     cerr << "Invalid value for nex: " << nex << ", terminating" << endl;
     throw("Invalid input argument");
     }
   }

 if(argc > 2)
   {
   istringstream ist(argv[2]);

   if( !(ist >> ney) )
     {
     cerr << "Invalid value for ney: " << argv[1] << ", terminating" << endl;
     throw("Invalid input");
     }

   if(ney <= 0) 
     {
     cerr << "Invalid value for ney: " << ney << ", terminating" << endl;
     throw("Invalid input argument");
     }
   }
 else
   ney = nex;

  nelem = 2*nex*ney;
  nodes = (nex+1)*(ney+1);
  nbelem = 2*(nex+ney);
  
  double hx = (xmax-xmin) / nex;
  double hy = (ymax-ymin) / ney;
  
  cout << "Boundary codes bcType[botom,rhs,top,lhs] = ["
       << bcType[0] << "," << bcType[1] << ","
       << bcType[2] << "," << bcType[3] << "]" 
       << endl;
  cout << "Number of elements in each x-direction = ......... " << nex 
       << endl;
  cout << "Number of elements in each y-direction = ......... " << ney 
       << endl;
  cout << "[xmin,xmax] = ............................... [" 
       << xmin << "," << xmax << "]" << endl;
  cout << "[ymin,ymax] = ............................... [" 
       << ymin << "," << ymax << "]" << endl;


  ostringstream fname;

  fname << "rectangle" << fixed << setprecision(3) << hx << ".tplc";

  ofstream out((fname.str()).c_str());

  if( !out.is_open() ) 
    {
      cerr << "dir3Mesh::input(): can't open .tplc file:" 
	   << fname.str() << endl;

      throw("Can't open file");
    }

  out << "TPLC" << endl 
      << "DIMENSION " << 2 << endl 
      << "POINTS " << nodes << endl;
   
  out.setf(ios::scientific,ios::floatfield);
  out.precision(16);

  for(j = 0; j < ney+1; j++)
  for(i = 0; i < nex+1; i++)
    {
      out << j * (nex+1) + i << " " << xmin+i*hx << " " << ymin+j*hy << endl;
    }

  out << "SEGMENTS " << nbelem << endl;   // next line is SEGMENTS nSegment

  int segId = 100;

  for(n = 0; n < nex; n++)
    {
      out << segId + bcType[0] << " 2 " << n << " " << n+1 << endl;

      segId += 2;
    }

  for(n = nex; n < nex+ney; n++)
    {
    out << segId + bcType[1] << " 2 " 
	<< ((n-nex)+1) * (nex+1) - 1 << " "
	<< ((n-nex)+2) * (nex+1) - 1 << endl;

    segId += 2;
    }

  for(n = nex+ney; n < 2*nex+ney; n++)
    {
    out << segId + bcType[2] << " 2 " 
	<< (nex+1)*ney + (n-nex-ney) + 1 << " "
	<< (nex+1)*ney + (n-nex-ney) << endl;

    segId += 2;
    }

  for(n = 2*nex+ney; n < 2*nex+2*ney; n++)
    {
    out << segId + bcType[3] << " 2 "
	<< (n-2*nex-ney + 1) * (nex+1) << " "
	<< (n-2*nex-ney    ) * (nex+1) << endl; 

    segId += 2;
    }

  out << "TRIANGLES " << nelem << endl;

  for(n = 0; n < nelem/2; n++)
    {
      int je = n / nex;
      int ie = n - je * nex;
      
      out << (nex+1) * je + ie << " " 
	  << (nex+1) * je + ie + 1 << " "
	  << (nex+1) * je + ie + (nex+1) << endl;

      out << (nex+1) * je + ie + 1 << " " 
	  << (nex+1) * je + ie + 1 + (nex+1) << " "
	  << (nex+1) * je + ie + 1 +  nex    << endl;
    }

  out.close();
 
  return(0); 
}
