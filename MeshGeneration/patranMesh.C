// Set up a mesh from patran.rpt for plasticity problem

#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<map>
#include<set>

#include"../LinearAlgebra/svector.C"
#include"gelement.C"
#include"../FiniteElement/basisFns.C"

//#define nDim 2 //change it, and gBasis(Square4) or Cube8
#if NDIM == 2
#define BASIS Square4
#elif NDIM ==3
#define BASIS Cube8
#endif

//#define MYPATRAN

//typedef GeometricElement<3, 8> gelement; //GeometricElement<Ndim, Nchi>
// typedef BoundaryGelement<2, 4> bgelement; //BoundaryGelement<Ndim, Nchi>

class patranMesh
{
public:
  patranMesh(std::string ss);

  patranMesh(const patranMesh &mm) 
  {
    std::cerr << "patranMesh(): don't copy a patranMesh\n";
    throw("copying a patranMesh");
  }
    
  patranMesh() 
  {
    std::cerr << "patranMesh(): Default constructor called, nyi\n";
    throw("bad patranMesh constructor");
  }

  ~patranMesh() {return;}

  int nNode()  const {return nNds;}
  int nElem()  const {return nElms;}
  int nBelem() const {return nbelem;}

  vector<double> params(); 
  void locate(const svector<NDIM>&, gelement&, svector<NDIM>&) const;
  void setBoundaryElems(int); 
  gelement intToGe(int);
  bgelement intToBge(int);
  void define_slip_planes(void);

  //information&utilities
  void info();
  void print_bdyf(std::ofstream&);
  std::string convertInt(int);
  void swap(vector<int>&, int,int);
  void construct_elm_lib();
  void check_numbering();
  void addto(set<vector<int> > s, vector<int> f);
  void print(set<int> s);
  vector<double> center_elem(int tag);
  int getNodeByElm(int);
  double setLayerElems(int, double, double, vector<double>& );
  int addtoLayerElms(int, double, double, vector<double>&);
  int getNumElmsInLayer() {return nex;}   //return number of 1 layer of elements
  int getCaseId(){return caseid;}
  int getLayer(){return layer;}
  int getSlipPlanes(){return slip_planes;}

  class geiterator
  {
  public:
    geiterator();
    geiterator(int l, patranMesh *g);

    gelement operator*();  
    void operator++();
    void operator++(int);
    geiterator& operator=(const geiterator & it);
    bool operator==(const geiterator & it);
    bool operator!=(const geiterator & it);
  private:
    int loc;
    patranMesh * gm;
  };
  
  geiterator gebegin();
  geiterator geend();
  
  class bgeiterator
  {
  public:
    bgeiterator();
    bgeiterator(int l, patranMesh * g);

    bgelement operator*();
    void operator++();
    void operator++(int);
    bgeiterator& operator=(const bgeiterator & it);
    bool operator==(const bgeiterator & it);
    bool operator!=(const bgeiterator & it);
  private:
    int loc;
    patranMesh * gm;
  };
  
  bgeiterator bgebegin();
  bgeiterator bgeend();

  LagrangeBasisFunction<NDIM> *gBasis;

private:
  int nElms;
  int nNds;
  int nbelem;
  
  int MPCs, Matrls, loads, Elprops, Groups, Points;
  
  int nex, ney;
  int bandth;
  int layer;
  int caseid;
  int slip_planes;  //each slip_Plane has layer layers of iso quad elems

  vector<double> max; 
  vector<double> min; 
  vector<double> cc;
  int bcType[4];

  map<int, vector<double> > nodes; //tag->coords
  map<int, vector<int> > elms; //tag->nodes
  map<int, vector<int> > belms; //tag->nodes
  map<int, std::string> topology; //tag->eleType
  map<int, std::string> btopology;//tag->beleType
  map<std::string, vector<int> > elm_lib;//eleTyp -> vector<number of nodes in element, number of nodes on face>
  map<std::string, std::set<vector<int> > > elmFaceLib; //eleType -> set of faces. face = set of node
  int input(std::string meshFile);

  map<int, pair<int, double> > layerElms; //<global eleID, pair<layerElmID, coordinate>
  std::vector<double> ys;
};


vector<double> patranMesh::params()
{
  double vals[] = {min[0],min[1],max[0],max[1],cc[0],cc[1]}; 
  return vector<double>(vals, vals + sizeof(vals) / sizeof(double) );
}

patranMesh::patranMesh(std::string ss) : gBasis(BASIS)//gBasis(Cube8) //gBasis(Square4) 
{ 
  construct_elm_lib();

  if(!ss.empty())
    {
      if(input(ss)!=0) std::cout<<"INPUT MESH FILE ERROR."<<std::endl;      
    }
  else
    {
      throw("patran mesh file cannot be empty.");
    }

  check_numbering();

  define_slip_planes();
  //setBoundaryElems();

  //info();
  
  //parameters for layer model 
}


void patranMesh::define_slip_planes()
{
  assert(nNds>0);
  //below is hard coding for each problem
  switch(nNds)
  {
  case 118546:
    caseid = 1000;
    break;
  case 235169:
    caseid = 2000;
    break;
  case 366278:
    caseid = 3061;
    break;
  case 357300:
    caseid = 3125;
    break;
  case 354338:
    caseid = 3186;
    break;
  case 352920:
    caseid = 3250;
    break;

  case 593162:
    caseid = 5000;
    break;

  case 118596:
    caseid = 1100;
    break;

    //not used below
  case 24666:
    caseid = 1;
    break;
  case 56595:
    caseid = 2;
    break;
  case 121278:
    caseid = 5;
    break;
  case 79572:
    caseid = 31;
    break;
  case 72722:
    caseid = 32;
    break;
  case 72830:
    caseid = 33;
    break;
  case 73020:
    caseid = 34;
    break;
  case 16366:
    caseid = -1;
    break;
  default:
    std::cerr<< "not implemented mesh in patranWallMesh.C."<<std::endl;
    throw("exception");
    break;
  }

  double ys_3_250[3] = {-250,0,250};   //case 34
  double ys_3_186[3] = {-186,0,186};   //case 33
  double ys_3_125[3] = {-125,0,125};   //case 32
  double ys_3_61[3] = {-61,0,61};      //case 31
  double ys_2[2] = {-166,166};             //case 2
  double ys_1[1] = {0};                   //case 1
  double ys_5[5] = {-332,-166,0,166,332}; //case 5
   
  //   assert(slip_planes<=10);
  switch(caseid)
    {
    case 1:
      slip_planes = 1;
      for(int i = 0;i<1;i++)
	ys.push_back(ys_1[i]);
      break;
    case 1000:
      slip_planes = 1;
      std::cerr<<"slip_planes = "<<slip_planes<<std::endl;
      for(int i = 0;i<1;i++)
	ys.push_back(ys_1[i]);
      break;
    case 1100:
      slip_planes = 1;
      for(int i = 0;i<1;i++)
	ys.push_back(ys_1[i]);
      break;

    case -1:
      slip_planes = 1;
      for(int i = 0;i<1;i++)
	ys.push_back(ys_1[i]);
      break;

    case 2: 
      slip_planes = 2;
      for(int i = 0;i<2;i++)
	ys.push_back(ys_2[i]);
      break;
    case 2000: 
      slip_planes = 2;
      for(int i = 0;i<2;i++)
	ys.push_back(ys_2[i]);
      break;

    case 5:
      slip_planes = 5;
      for(int i = 0;i<5;i++)
	ys.push_back(ys_5[i]);     
      break;
    case 5000:
      slip_planes = 5;
      for(int i = 0;i<5;i++)
	ys.push_back(ys_5[i]);     
      break;

    case 31:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_61[i]);
     break;
    case 3061:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_61[i]);
     break;
    case 32:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_125[i]);
     break;
    case 3125:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_125[i]);
     break;
    case 33:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_186[i]);
     break;
    case 3186:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_186[i]);
     break;
    case 34:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_250[i]);
     break;
    case 3250:
      slip_planes = 3;
      for(int i = 0;i<3;i++)
	ys.push_back(ys_3_250[i]);
     break;

    default:
      std::cerr<<"not implemented mesh"<<std::endl;
      throw("not implemented");
      break;
   }
 
  //  layer = slip_planes;

  std::cout<<"ys arr = ";
  for(int i = 0;i<slip_planes;i++)
    std::cout<<ys[i]<<"  "; 
  std::cout<<std::endl;
}


int patranMesh::getNodeByElm(int tag)
{
  if(layerElms.find(tag)==layerElms.end())
    return -1;
  else
    return layerElms[tag].first;
}


int cmpx(const std::pair<int,vector<double> > &x, const std::pair<int,vector<double> > &y)  
{    
  //vector<double> xv, yv;
  
  //xv = x.second; yv = y.second;
  assert(x.second.size()>0); assert(y.second.size()>0);
  //std::cout<<"xv0  = "<<xv[0]<<"yv0 = "<<yv[0]<<std::endl;
  
  //return (x.second).front < (y.second).front;  
  return (x.second)[0]<(y.second)[0];
}  

int cmpy(const std::pair<int,std::vector<double> > &x, const std::pair<int,std::vector<double> > &y )  
{  

  //vector<double> xv, yv;

  //xv = x.second; yv = y.second;
  assert(x.second.size()>0); assert(y.second.size()>0);
  //std::cout<<"xv1  = "<<xv[1]<<"yv1 = "<<yv[1]<<std::endl;
  //return (x.second).front < (y.second).front;  
  //return xv[1]<yv[1];
  return (x.second)[1] < (y.second)[1];  
}  



double patranMesh::setLayerElems(int layer, double layer_top, double layer_bot, vector<double>& imagNodes)
{
  assert(layer%2==1);
  double dx_min = 1e50;
  assert(imagNodes.size() == 0);
  double halfbandwidth = (fabs(layer_top)+fabs(layer_bot))/2.0;
  map<int,vector<int> > :: iterator it;

  vector<pair<int,vector<double> > > sortedLayerElm; //pair<tag, cc>
  //map<int,vector<double> > lymts; //tag->cc
  //add the middle layer of elems, and set the 1d node`s coordinates.
  for(int ind_ys = 0;ind_ys<slip_planes;ind_ys++)
    {
      std::cerr<<"ind_ys = "<<ind_ys<<"; layer = "<<layer<<std::endl;
      for(it = elms.begin();it!=elms.end();it++)
	{
	  vector<double> v = center_elem(it->first);
	  
	  //std::cerr<<"v[1] = "<<v[1]<<"; ys = "<<ys[ind_ys]+layer_top <<"; ys_bot = "<<ys[ind_ys]-layer_bot<<std::endl;
	  if((v[1]<ys[ind_ys]+halfbandwidth) && (v[1]>ys[ind_ys]-halfbandwidth))//0.001 is a tolerance
	    {	      
	      std::pair<int,std::vector<double> > p; 
	      p = std::make_pair(it->first, v);
	      sortedLayerElm.push_back(p);  
	      //lymts.insert(std::pair<int, vector<double> >(it->first, v));
	    }
	}
    }
  //assert(lymts.size()%layer == 0);  

  // for (map<int,vector<double> >::iterator curr = lymts.begin(); curr != lymts.end(); ++curr)  
  //   {  
  //   }  
  assert(sortedLayerElm.size()%(layer*slip_planes) == 0);
  std::cout<<"szImagNodes = " << sortedLayerElm.size() <<"; layer = "<<layer<<"; slip_planes = "<<slip_planes<<std::endl;
  int szImagNodes = sortedLayerElm.size()/(layer*slip_planes);

  std::sort(sortedLayerElm.begin(), sortedLayerElm.end(), cmpy); 
  
  //for(vector<pair<int,vector<double> > >::iterator vit = sortedLayerElm.begin(); vit!=sortedLayerElm.end();vit++)  std::cout<<vit->first<<"--->"<<vit->second[1]<<std::endl;

  //std::vector<pair<int,vector<double> > >subVec; 

  for(int i = 0;i<layer*slip_planes;i++)			  
    {
      //std::cout<<"loop of layer = "<<i<<"; szImagNodes*i = "<< szImagNodes*i <<"; sortedLayaerElm.size() = "<<sortedLayerElm.size()<<std::endl;
      //" sortedLayerElm.begin()+szImagNodes*(layer+1)-1 = "<<szImagNodes*(i+1)-1<<"; (sortedLayerElm.begin()+szImagNodes*i = "<<*(sortedLayerElm.begin()+szImagNodes*i).first <<"; sortedLayerElm.begin()+szImagNodes*(i+1)-1 = "<<*(sortedLayerElm.begin()+szImagNodes*(i+1)-1).first<<std::endl;
      // for(int j = i*szImagNodes;j<szImagNodes*(i+1);j++)
      // 	{
      // 	  subVec.push_back(sortedLayerElm[j]);
      // 	  //std::cout<<"sortedLayerElm["<<j<<"] = "<<sortedLayerElm[j].first<<std::endl;
      // 	}
      // std::sort(subVec.begin(), subVec.end(),cmpx);
      
      // for(int j = 0;j<szImagNodes;j++)
      // 	{
      // 	  sortedLayerElm[j] = subVec[j];
      // 	}
      // subVec.empty();

      std::sort(sortedLayerElm.begin()+i*szImagNodes, sortedLayerElm.begin()+(i+1)*szImagNodes, cmpx);

      this->nex = szImagNodes;
    }			  

	
  //  std::cout<<"after x sort"<<std::endl;		
  //for(vector<pair<int,vector<double> > >::iterator vit = sortedLayerElm.begin(); vit!=sortedLayerElm.end();vit++) std::cout<<vit->first<<"---> "<<vit->second[0]<<std::endl;

  dx_min = (max[0]-min[0])/(double)szImagNodes;

  std::cerr<<"dx_min first = "<<dx_min<<std::endl;

  for(int j = 0;j<slip_planes;j++)    
    for(int i = 0;i<szImagNodes;i++)	
      {	  	
	if(i == 0)
	  imagNodes.push_back(min[0]+dx_min/2.0);
	else
	  imagNodes.push_back(imagNodes[i-1]+dx_min);

	for(int k = 0;k<layer;k++)
	  layerElms.insert(std::pair<int,pair<int,double> >(sortedLayerElm[(j*layer+k)*szImagNodes+i].first,std::pair<int,double>(j*szImagNodes+i, imagNodes[i]) ));
      }    


  std::ofstream dump("layerMeshCheck.rpt");

  dump.setf(std::ios::fixed, std::ios::floatfield);
  for(int i = 0;i<layer;i++)
    {
      for(int j = i*szImagNodes;j<szImagNodes*(i+1);j++)
	dump<<sortedLayerElm[j].first+1<<"  ";
      for(int k = 0;k<10;k++) dump<<"tab next line"<<std::endl;
    }

  std::cerr<<"szImagNodes = "<<szImagNodes<<std::endl;
  for(int i = 1;i<szImagNodes;i++)
    if(dx_min>(imagNodes[i]-imagNodes[i-1]))  
      dx_min = imagNodes[i]-imagNodes[i-1];
     
  return dx_min;
}


vector<double> patranMesh::center_elem(int tag)
{
  assert(tag>=0 && (unsigned)tag<elms.size());
  vector<double> v;
  vector<int> ns = elms[tag];
  for(int i = 0;i<NDIM;i++)
    {
      double val = 0;
      for(int j = 0;(unsigned)j<ns.size();j++)
	{
	  val+= nodes[ns[j]][i];
	}
      val/=(double)ns.size();
      v.push_back( val);
    }
  return v;
}

//add more element categories here  
void patranMesh::construct_elm_lib()
{
#if NDIM == 2  
  vector<int> myvector(2);
  myvector[0] = 4; myvector[1] = 2; //nodes in elem, number of nodes on face
  elm_lib["Quad4"] = myvector;
  std::set< vector<int> > faces;
  int ints[4][2] = {{0,1},{1,2},{2,3},{3,0}};
  for(int i = 0;i<4;i++)
    {      
      vector<int> subVec (ints[i],ints[i]+2); 
      faces.insert(subVec);
    }
  elmFaceLib["Quad4"] = faces;
 
#elif NDIM == 3
  //Hex8
  vector<int> myvector(2);
  myvector[0] = 8; myvector[1] = 4;
  elm_lib["Hex8"] = myvector;  
  int ints[6][4] = {{6,5,1,2}, {6,7,4,5}, {2,6,7,3}, {7,4,0,3}, {1,2,3,0}, {5,1,0,4}};
  std::set< vector<int> > faces;
  for(int i = 0;i<6;i++)
    {
      vector<int> subVec (ints[i], ints[i]+4);
      faces.insert(subVec);
    }
  elmFaceLib["Hex8"] = faces;

#endif
}

//8624  8756 4873 7513  2431 6215

int patranMesh::input(std::string rpt)
{
  std::cout<<"INPUT MESH FILE "<<rpt<<"..."<<std::endl;
  
  std::ifstream in(rpt.c_str());
  if(!in.is_open())
    {
      std::cerr<<"can`t open files:"<< rpt<<std::endl;
      throw("can`t open file");
    }

  string line;
  std::string summary("Total number of entities:");
  std::string ElmAttr("Element Attributes:");
  std::string ElmConn("Element Connectivity:");
  std::string NdAttr("Node Attributes (Coordinates are in Reference CF)");
  std::string PtAttr("Point Attributes (Coordinates are in Global CF)");
  std::string BelmConn("Boundary Element Connectivity");

  int cnt_Sum = 0;
  int cnt_EAt = 0;
  int cnt_ECo = 0;
  int cnt_NAt = 0;
  int cnt_PAt = 0;
  int cnt_BCo = 0;

  bool found_Sum = false;
  bool found_EAt = false;
  bool found_ECo = false;
  bool found_NAt = false;
  bool found_PAt = false; 
  bool found_BCo = false;

  while(getline(in, line)) 
    {
      if(in.eof() || line.empty()) break;

      //std::cout<<line<<std::endl;
      std::stringstream ist(line);

      if(!found_Sum && line.find(summary)!=std::string::npos)
	{
	  cnt_Sum = 0;
	  found_Sum = true;
	  //std::cout<<"found sum"<<std::endl;
	}

      if(!found_EAt && line.find(ElmAttr)!=std::string::npos)
	{
	  cnt_EAt = 0;
	  cnt_Sum = 3; //stop read sum
	  found_EAt = true;
	  //std::cout<<"found eat"<<std::endl;
	}


      if(!found_ECo && line.find(ElmConn)!=std::string::npos)
	{
	  cnt_ECo = 0;
	  cnt_EAt = 3; //stop read EAt
	  found_ECo = true;
	  //std::cout<<"found eco"<<std::endl;
	}

      if(!found_NAt && line.find(NdAttr)!=std::string::npos)
	{
	  cnt_NAt = 0;
	  cnt_ECo = 3; //stop read ECo
	  found_NAt = true;
	  //std::cout<<"found nat"<<std::endl;
	}

      if(!found_PAt && line.find(PtAttr)!=std::string::npos)
	{
	  cnt_PAt = 0;
	  cnt_NAt = 3; //stop read NAt
	  found_PAt = true;
	  //std::cout<<"found pat"<<std::endl;
	}

      if(!found_BCo && line.find(BelmConn)!=std::string::npos)
	{
	  cnt_BCo = 0;
	  cnt_PAt = 3; //stop read NAt
	  found_BCo = true;
	  //std::cout<<"found pat"<<std::endl;
	}

      //if( found_Sum && found_EAt &&  found_ECo &&  found_NAt &&  found_PAt &&  found_BCo &&;)
      
      //read summary
      if(found_Sum && cnt_Sum==2)
	{
	  ist>>nNds>>nElms>>MPCs>>Matrls>>loads>>Elprops>>Groups>>Points;
	  nbelem = 0;
	  //std::cout<<"summary = "<<setw(7)<<nNds<<setw(7)<<nElms<<setw(7)<<MPCs<<setw(7)<<Matrls<<setw(7)<<loads<<setw(7)<<Elprops<<setw(7)<<Groups<<setw(7)<<Points<<setw(7)<<std::endl;
	}

      //read topology
      if(found_EAt && cnt_EAt==2)
	{
	  int tmpInt;
	  string tmpStr;
	  ist>>tmpInt;
	  ist>>tmpStr;
	  if(elm_lib.find(tmpStr)==elm_lib.end())
	    {
	      std::cout<<"ELEM TYPE: "<<tmpStr<<"NEED TO BE IMPLEMENTED"<<std::endl;
	      throw("not implemented.");
	    }
	  //std::cout<<"tmpInt = "<<tmpInt<<" tmpStr = "<<tmpStr<<std::endl;
	  int i = tmpInt-1;
	  topology[i] = tmpStr;
	  //NDIM = elm_lib[tmpStr][1]; 
	  continue;
	}
      //read Connectivity
      if(found_ECo && cnt_ECo==2)
	{
	  int tmpID;
	  int tmpInt;
	  ist>>tmpID;
	  int i = tmpID-1;
	  //assert(tmpID = (i+1));
	  int n = elm_lib[topology[i]][0]; 
	  //std::cout<<"n = "<<n<<std::endl;
	  vector<int> nodes_in_elm(n);
	  for(int j = 0;j<n;j++)
	    {		 
	      ist>>tmpInt;
	      //std::cout<<"tmpInt = "<<tmpInt<<"  ";
	      nodes_in_elm[j] = tmpInt-1; //start from 0.
	    }
	 
	  // if(topology[i].compare("Hex8")==0)
	  //   {
	  //     swap(nodes_in_elm,2,3);
	  //     swap(nodes_in_elm,6,7);
	  //   }
	  //std::cout<<std::endl;
	  elms[i] = nodes_in_elm;	    
	}

      //read Nodes
      if(found_NAt && cnt_NAt == 2)
	{
	  //std::cout<<"this line = "<<line<<std::endl;
	  vector<double> coords(NDIM);
	  int tmpInt;
	  double tmpDouble;
	  ist>>tmpInt;
	  int i = tmpInt - 1;
	  //assert(tmpInt == (i+1));
	  for(int j = 0;j<5;j++)
	    {		  
	      if(j<NDIM)
		{
		  ist>>tmpDouble;
		  coords[j] = tmpDouble;
		}		    
	    }
	      
	  nodes[i] = coords;	    
	  //std::cout<<"i = "<<i <<" coords[0] = "<<coords[0]<<" coords[1] = "<<coords[1]<< " coords[2] = "<<coords[2]<<std::endl;
	}

      if(found_BCo && cnt_BCo == 2)
	{
	  int tmpID;
	  ist>>tmpID;
	  int i = tmpID - 1;

	  string tmpStr;
	  ist>>tmpStr;
	 
	  int n = elm_lib[tmpStr][0];
	  vector<int> nodes_in_elm(n);
	  
	  //j = 0: bc type. j =1--4
	  for(int j = 0;j<n+1;j++)
	    {
	      double tmpInt;
	      ist>>tmpInt;		
	      nodes_in_elm[j] = tmpInt-1;
	    }			 
	  belms[i] = nodes_in_elm;
	  nbelem++;
	}
    
      cnt_Sum = (cnt_Sum==2)?cnt_Sum:cnt_Sum+1; 
      cnt_EAt = (cnt_EAt==2)?cnt_EAt:cnt_EAt+1;
      cnt_ECo = (cnt_ECo==2)?cnt_ECo:cnt_ECo+1;
      cnt_NAt = (cnt_NAt==2)?cnt_NAt:cnt_NAt+1;      
      cnt_PAt = (cnt_PAt==2)?cnt_PAt:cnt_PAt+1;
      cnt_BCo = (cnt_BCo==2)?cnt_BCo:cnt_BCo+1;
    }
  
  in.close();

  max = vector<double>(NDIM);
  min = vector<double>(NDIM);
  cc = vector<double>(NDIM);
  for(int i = 0;i<NDIM;i++)
    {
      //std::cout<<"i ="<<i<<std::endl;
      min[i] = 1e50;
      max[i] = -1e50;

      map<int,vector<double> >::iterator it;
      for(it=nodes.begin();it!=nodes.end();it++)
	{	  
	  min[i] = (it->second)[i]<min[i]?(it->second)[i]:min[i];
	  max[i] = (it->second)[i]>max[i]?(it->second)[i]:max[i];
	}
    }

  for(int i = 0;i<NDIM;i++) cc[i] = (max[i]+min[i])/2.0;

  // std::cout<<"nodes.size() = "<<nodes.size()<<" nNds = "<<nNds<<std::endl;
  // std::cout<<"elms.size() = "<<elms.size()<<" nElms = "<<nElms<<std::endl;
 
  //size
  assert (nodes.size() == (unsigned)nNds);
  assert (elms.size() == (unsigned)nElms);

  //topo
  assert( topology.size() == (unsigned)nElms);
  //std::cout<<"top.size = "<<topology.size()<<std::endl;
  return 0;
}

void patranMesh::check_numbering()
{
  std::cout<<"START TO CHECK NUMBERING OF ELEMENTS"<<std::endl;

  map<int, vector<int> > :: iterator it;
  
  double cross;
  vector<double> v0, v1, v2, v3;
  //double x0, x1, x2, x3;
  //double y0, y1, y2, y3;
  vector<double> w(2);
  vector<double> z(2); 
  int cnt = 0;
  for(it = elms.begin();it!=elms.end();it++)
    {
      int index = it->first;
      vector<int> nds = it->second;
      if(topology[index] == "Quad4"){
	
	v0 = nodes[nds[0]];
	v1 = nodes[nds[1]];
	v2 = nodes[nds[2]];
	v3 = nodes[nds[3]];
	//     v[j] = nodes[nds[j]];
	//     std::cout<<"start"<<std::endl;
	//     std::cout<<"v["<<j<<"]" <<v[j][0]<<"; "<<v[j][1]<<std::endl;
	//     for(int i = 0;i<2;i++)
	//       cc[i] += v[j][i];
	//   }
	// std::cout<<"v0 = "<<v0[0]<<"; "<<v0[1]<<std::endl;
	// std::cout<<"v0 = "<<v1[0]<<"; "<<v1[1]<<std::endl;
	// std::cout<<"v0 = "<<v2[0]<<"; "<<v2[1]<<std::endl;
	// std::cout<<"v0 = "<<v3[0]<<"; "<<v3[1]<<std::endl;
	for(int i = 0;i<2;i++)
	  {
	    w[i] = v2[i] - v1[i];
	    z[i] = v3[i] - v2[i];
	  }
	
	cross = w[0]*z[1] - w[1]*z[0];
	if(cross<0)
	  {
	    //std::cout<<"cross < 0"<<std::endl;
	    vector<int> new_nodes;
	    for(int i =nds.size()-1;i>=0;i--)
	      new_nodes.push_back(nds[i]);
	    elms[index] = new_nodes;
	    cnt++;
	  }
      }
      if(topology[index] == "Hex8"){
	
	v0 = nodes[nds[0]];
	v1 = nodes[nds[1]];
	v2 = nodes[nds[2]];
	v3 = nodes[nds[3]];
	// v4 = nodes[nds[4]];
	// v5 = nodes[nds[5]];
	// v6 = nodes[nds[6]];
	// v7 = nodes[nds[7]];

	//     v[j] = nodes[nds[j]];
	//     std::cout<<"start"<<std::endl;
	//     std::cout<<"v["<<j<<"]" <<v[j][0]<<"; "<<v[j][1]<<std::endl;
	//     for(int i = 0;i<2;i++)
	//       cc[i] += v[j][i];
	//   }
	// std::cout<<"v0 = "<<v0[0]<<"; "<<v0[1]<<std::endl;
	// std::cout<<"v0 = "<<v1[0]<<"; "<<v1[1]<<std::endl;
	// std::cout<<"v0 = "<<v2[0]<<"; "<<v2[1]<<std::endl;
	// std::cout<<"v0 = "<<v3[0]<<"; "<<v3[1]<<std::endl;
	for(int i = 0;i<2;i++)
	  {
	    w[i] = v2[i] - v1[i];
	    z[i] = v3[i] - v2[i];
	  }
	
	cross = w[0]*z[1] - w[1]*z[0];
	if(cross<0)
	  {
	    //std::cout<<"cross < 0"<<std::endl;
	    vector<int> new_nodes;
	    for(int i =nds.size()-1;i>=0;i--)
	      new_nodes.push_back(nds[i]);
	    elms[index] = new_nodes;
	    cnt++;
	  }
      }    
    }
  std::cout<<"NO. OF ELEMENTS THAT CHANGE NUMBERING = "<<cnt<<std::endl;
}

void patranMesh::setBoundaryElems(int rb)
{

  if(rb == 1)
    {
      belms.clear();
      std::cout<<"INPUT BDY FILE: bdy.out"<<std::endl;
      std::ifstream in("bdy.out");
      if(!in.is_open())
	{
	  std::cerr<<"cannot open file: bdy.out"<<std::endl;
	  throw("cannot open bdy.out");
	}
      else
	{
	  int k1,k2,k;
	  std::string line;
	  while(getline(in,line))
	    {
	      if(in.eof()||line.empty()) break;
	      std::stringstream ist(line);
	      ist>>k1; //be id
	      ist>>k2; //no of nodes on face
	      vector<int> nodes_in_elm;
	      for(int i = 0;i<k2;i++)
		{
		  ist>>k;
		  nodes_in_elm.push_back(k);
		}
	      belms[k1] =  nodes_in_elm;
	    }
	}
      nbelem = belms.size();
      std::cout<<"belms.siz() = "<<nbelem<<std::endl;
    }
  else
    {
      std::cout<<"SET BOUNDARY ELEMENTS"<<std::endl;

      set< vector<int> > face_pool;
  
      set< vector<int> > :: iterator sit_faces;
      map<int,std::string >:: iterator it;
      set< vector<int> > :: iterator sit;
      int elemID = 0;

      for(it=topology.begin();it!=topology.end();it++)
	{
	  if((++elemID)%500==0) 
	  std::cout<<"ELEM ID = "<<elemID<<std::endl;

	  for(sit_faces = elmFaceLib[it->second].begin();sit_faces!=elmFaceLib[it->second].end();sit_faces++)
	    {
	      vector<int> ref_face = *sit_faces; //raw index of i th face node
	      vector<int> face;  //real face 
	      vector<int> ::iterator it_faceIndex;
	      for(it_faceIndex = ref_face.begin();it_faceIndex!=ref_face.end();it_faceIndex++)
		{
		  //std::cout<<"elms["<<it->first<<"]["<<*it_faceIndex<<"] = "<<elms[it->first][*it_faceIndex]<<"; ";
		  face.push_back(elms[it->first][*it_faceIndex]);
		}
	      //std::cout<<std::endl;
	  
	      set<int> fs(face.begin(),face.end());
	      //std::cout<<"fs = ";print(fs);
	      set<vector<int> > ::iterator pit;
	      bool insertFlag = true;
	      for(pit=face_pool.begin();pit!=face_pool.end();pit++)
		{
		  set<int> s ( (*pit).begin(), (*pit).end());
		  //std::cout<<"s = ";print(s);
		  if(s == fs)
		    {
		      //std::cout<<"s == fs"<<std::endl;
		      face_pool.erase(pit);
		      insertFlag = false;
		    }
		  //std::cout<<"step"<<std::endl;
		}
	      if(insertFlag)  face_pool.insert(face);
	    }
	}

      int cnt = 0;
      for(sit=face_pool.begin();sit!=face_pool.end();sit++)
	{
	  //std::cout<<"sit->second = "<<(*sit)[0]<<(*sit)[1]<<std::endl;
	  vector<int> face = (*sit);
	  belms[cnt] = face;
	  cnt++;
	}
      nbelem = cnt;

      std::cout<<"boundary elements number = "<<belms.size()<<std::endl;
      
      // for(int i = 0;i<belms.size();i++)
      // 	{
      // 	  for(int j = 0;j<4;j++)
      // 	    std::cout<<"  belms["<<i<<"]["<<j<<"] = "<<belms[i][j] ;
      // 	  std::cout<<std::endl;
      // 	}
      // std::cout<<"end of out put b.elm"<<std::endl;
    }
}

void 
patranMesh::print_bdyf(std::ofstream &out)
{
  for(std::map<int,std::vector<int> >::iterator it = belms.begin();it!=belms.end();it++)
    {
      out<< it->first <<"  "<<(it->second).size()<<"  ";
      for(int j = 0;j<(int)(it->second).size();j++)
	out<<(it->second)[j]<<"  ";
      out<<std::endl;
    }
}


void patranMesh::print(set<int> s)
{
  set<int>::iterator it;
  for(it=s.begin();it!=s.end();it++)
    std::cout<<*it<<"  ";
  std::cout<<std::endl;
}

void patranMesh::swap(vector<int>& v, int i, int j)
{
  double swapV;
  swapV = v[i];
  v[i] = v[j];
  v[j] = swapV;
}

std::string patranMesh::convertInt(int number)
{
  stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

void patranMesh::info()
{
  //  for(int i = 0;i<nElms;i++) std::cout<<topology[i]<<std::endl;
 
  for(int i = 0; i<NDIM;i++)
    {
      std::cout<<"min["<<i<<"] = " <<min[i]<<"; max["<<i<<"] = " <<max[i]<<"; cc["<<i<<"] = "<<cc[i]<<std::endl;
    }

  // std::cout << "patranMesh::info(): patranrilateral mesh" << std::endl;
  // std::cout << "Boundary elements number = ["<<this->nbelem<<"]" << std::endl;

  std::ofstream dump("meshCheck.rpt");

  dump.setf(std::ios::fixed, std::ios::floatfield);
  dump<<"  eeID -> nodes"<<std::endl;
  for (int i = 0;i<nElms;i++ )
    {
      dump<<setw(6)<<i<<setw(8)<<topology[i];
      for(int j = 0;j<elm_lib[topology[i]][0];j++)
	dump<<setw(6)<<elms[i][j];
      dump<<std::endl;
    }
  dump<<"  nodeID -> coords"<<std::endl;

  //std::cout<<"nodes.size()" << nodes.size()<<std::endl;
  for(int i = 0;i<nNds;i++)
    {
      dump<<setw(6)<<i;
      for(int j = 0;j<NDIM;j++)
	dump<<setw(15)<<nodes[i][j];
	//std::cout<<"nodes["<<i<<"]["<<j<<"]= " <<nodes[i][j]<<"  ";
	dump<<std::endl;
	//std::cout<<std::endl;
    } 
  //std::cout<<"get here"<<std::endl;
  return;
}

void patranMesh::locate(const svector<NDIM>& xx, gelement& ee, svector<NDIM>& xi) const
{
  throw("not yet implemented");
}


patranMesh::geiterator::geiterator() {
  loc = 0;
  gm = NULL;
}

patranMesh::geiterator::geiterator(int l, patranMesh * g) {
  loc = l;
  gm = g;
}

void patranMesh::geiterator::operator++() {
  (*this)++;
}

void patranMesh::geiterator::operator++(int) {
  loc++;
}

patranMesh::geiterator& patranMesh::geiterator::operator=(const geiterator & it) {
  loc = it.loc;
  gm = it.gm;
  return *this;
}

bool patranMesh::geiterator::operator==(const geiterator & it) {
  if (gm != it.gm || loc != it.loc) {
    return false;
  }
  return true;
}

bool patranMesh::geiterator::operator!=(const geiterator & it) {
  return !(*this == it);
}

gelement patranMesh::geiterator::operator*() {

  gelement ee = gm->intToGe(loc);
    
  return ee;
}


patranMesh::geiterator patranMesh::gebegin() {
  return geiterator(0, this);
}

patranMesh::geiterator patranMesh::geend() {
  return geiterator(nElms, this);
}

patranMesh::bgeiterator::bgeiterator() {
  loc = 0;
  gm = NULL;
}

patranMesh::bgeiterator::bgeiterator(int l, patranMesh * g) {
  loc = l;
  gm = g;
}

void patranMesh::bgeiterator::operator++() {
  (*this)++;
}

void patranMesh::bgeiterator::operator++(int) {
  loc++;
}

patranMesh::bgeiterator& patranMesh::bgeiterator::operator=(const bgeiterator & it) {
  loc = it.loc;
  gm = it.gm;
  return *this;
}

bool patranMesh::bgeiterator::operator==(const bgeiterator & it) {
  if (loc != it.loc || gm != it.gm) {
    return false;
  }
  return true;
}

bool patranMesh::bgeiterator::operator!=(const bgeiterator & it) {
  return !(*this == it);
}

bgelement patranMesh::bgeiterator::operator*() {

  bgelement be = gm->intToBge(loc);

  return be;

}

patranMesh::bgeiterator patranMesh::bgebegin() {
  return bgeiterator(0,this);
}

patranMesh::bgeiterator patranMesh::bgeend() {
  return bgeiterator(nbelem,this);
}


gelement patranMesh::intToGe(int n)
{
  gelement ee;
  if((n < 0) || (n >= nElms))
    {
      std::cerr << "Element number " << n << " out of range, nElms = " 
		<< nElms << std::endl;

      throw("Element number out of range");
    }

  vector<int> v = elms[n];

  int num_elm_nodes = v.size();
  for(int i = 0;i<num_elm_nodes;i++)
    {
      ee.nodes[i] = elms[n][i];
      for(int j = 0;j<NDIM;j++)
	ee.xx[i][j] = nodes[elms[n][i]][j];
    }
 
  ee.master = gBasis->master;
  ee.basis  = gBasis->basis;
  ee.tag = n;

  return(ee);
}

bgelement patranMesh::intToBge(int n)
{
  if((n < 0) || (n >= nbelem))
    {
      std::cerr << "Boundary element number " << n << " out of range, nElms = " 
		<< nbelem << std::endl;

      throw("Element number out of range");
    }

  bgelement be;
  vector<int> v = belms[n];
  
  int num_elm_nodes  = v.size();

  for(int i = 0;i<num_elm_nodes;i++)
    {
      be.nodes[i] = belms[n][i];
      for(int j = 0;j<NDIM;j++)
	{
	  be.xx[i][j] = nodes[belms[n][i]][j];	  
	}
      //std::cout<<"intobge n = "<<n <<"be.nodes["<<i<<"]="<<be.nodes[i]<<std::endl;
    }

  if(fabs(be.xx[0][0]-max[0])<1e-4 && fabs(be.xx[1][0]-max[0])<1e-4)
    {          
      be.type = 0;
    }
  else    
    be.type = 1; //1: neumann. 0: dirichlet

  be.type = 1;

  be.master = gBasis->boundary->master;
  be.basis  = gBasis->boundary->basis;
  be.tag = n;

  return(be);
}



