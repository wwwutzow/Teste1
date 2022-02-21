#include "NOS.h"

CNo::CNo(const int& _NumberNo, const double& _x1, const double& _x2, const double& _y1, const double& _y2)
{
  NumberNo=_NumberNo;
  X[0][0]=_x1;
  X[0][1]=_x2;
  X[1][0]=_y1;
  X[1][1]=_y2;
  CondCountour[0]=0;
  CondCountour[1]=0;
  vCondCountour[0]=0.0;
  vCondCountour[1]=0.0;
  PConcentrated[0]=0.0;
  PConcentrated[1]=0.0;
  Adress[0]=0;
  Adress[1]=0;
  FP=0;
  no_E.resize(3*3,0.0);
  no_S.resize(3*3,0.0);
  no_ue=0.0;
}

CNo::~CNo()
{
  no_E.clear();
  no_S.clear();
}

void Le_Nos(tvNo& No, const std::string& NAr)
{
  int nNo;            ///< Numero de nos a ser lido
  char s[1000];       ///< Variavel auxiliar
  std::ifstream Ent;  ///< Arquivo de leirura
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  nNo; Ent.getline(s,1000);
  Ent.getline(s,1000);
  for (int i=0; i<nNo; i++) 
  {
    double x[2],y[2];
	int NumberNo;
	Ent >> NumberNo >> x[0] >> x[1]; Ent.getline(s,1000);
    y[0]=x[0];
	y[1]=x[1];
	No.push_back(CNo(NumberNo,x[0],x[1],y[0],y[1]));
  }
  Ent.close();
}

void Le_CC(tvNo& No, const std::string& NAr, int& ngl)
{
  int nCC;            ///< Numero de nos a ser lido
  char s[1000];       ///< Variavel auxiliar
  std::ifstream Ent;  ///< Arquivo de leirura
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  nCC; Ent.getline(s,1000);
  Ent.getline(s,1000);
  for (int i=0; i<nCC; i++) 
  {
    int no;
	int dir;
	double vCondCountour;
	Ent >> no >> dir >> vCondCountour; Ent.getline(s,1000);
	No[no].w_CondCountour(dir,1);
	No[no].w_vCondCountour(dir,vCondCountour);
	No[no].w_X(1,dir,vCondCountour);
  }
  int k,n;
  k=0;
  n=0;
  for (int i=0; i<No.size(); i++)
  {
	for (int j=0; j<2; j++)
    {
	  if (No[i].r_CondCountour(j)==0)
	  {
		k=k+1;
		No[i].w_Adress(j,k-1);
	  }
      else
	  {
		No[i].w_Adress(j,-1);
	  }
    }
  }
  ngl=k;
  Ent.close();
}


void Le_P(tvNo& No, const std::string& NAr)
{
  int nP;            ///< Numero de nos a ser lido
  char s[1000];       ///< Variavel auxiliar
  std::ifstream Ent;  ///< Arquivo de leirura
  Ent.open(NAr.c_str());
  Ent.getline(s,1000);
  Ent >>  nP; Ent.getline(s,1000);
  Ent.getline(s,1000);
  for (int i=0; i<nP; i++) 
  {
    int no;
    int dir;
	double PConcentrated;
	Ent >> no >> dir >> PConcentrated; Ent.getline(s,1000);
	No[no].w_PConcentrated(dir,PConcentrated);
  }
  Ent.close();
}


